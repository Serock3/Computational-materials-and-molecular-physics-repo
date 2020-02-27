#%%
from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
from ase.units import Hartree
from gpaw.lrtddft import LrTDDFT


def discrete_spectrum(lr, fpath):
    '''
    Returns the discrete spectrum from a GPAW LrTDDFT calculator.

    Inputs:
    lr:		GPAW LrTDDFT calculator
            If you saved the calculator to file (i.e. calcFile.dat)
            you can load it via calc=LrTDDFT('calcFile.dat')
    fpath:	Path and filename you wish to write the discrete spectrum to

    The file will contain the energy for each excitation, the average
    oscillator strength, and the (x,y,z) oscillator strengths.
    '''
    with open(fpath, 'w') as f:
        f.write('# %4s %12s %12s %12s %12s %12s\n' % ('No', 'Energy (eV)', 'osc', 'osc x', 'osc y', 'osc z'))
        for i, exc in enumerate(lr):
            energy = exc.get_energy() * Hartree
            osc_j = exc.get_oscillator_strength()
            data = [i, energy] + list(osc_j)
            f.write('%6d %12.5f %12.6e %12.6e %12.6e %12.6e\n' % tuple(data))


def fold(x_t, x_i, y_i, width):
    '''
        Convolutes each peak in the discrete spectrum
        with a Gaussian of chosen width.
        
        inputs:
        x_t:    vector with energies (i.e. linspace/np.arange)
        x_i:    stick spectrum energies
        y_i:    stick spectrum intensities
        width:  width of the Gaussian

        outputs:
        y_t:    convoluted spectrum, i.e. intensities
    '''
    def Gauss(x0):
        norm = 1.0 / (width * np.sqrt(2 * np.pi))
        return norm * np.exp(-0.5 * (x_t - x0)**2 / width**2)

    y_t = np.zeros_like(x_t)
    for x, y in zip(x_i, y_i):
        y_t += y * Gauss(x)
    return y_t


def dump_data(lr, fpath):
	'''
		NOTE: This dumps everything in atomic units

		Inputs:
		lr:		GPAW LrTDDFT calculator
					If you saved the calculator to file (i.e. calcFile.dat)
					you can load it via calc=LrTDDFT('calcFile.dat')
		fpath:	Path and filename you wish to dump to

		Keywords for dumped quantities:
		'i_p':		Kohn-Sham state that gets excited from
		'a_p':		Kohn-Sham state that gets excited to
		'ediff_p':	Kohn-Sham eigenvalue differences for each KS excitation
		'fdiff_p':	Kohn_sham occupation differe for each KS excitation
		'mux_p':	x-component of dipole matrix elements
		'muy_p':	y-component of dipole matrix elements
		'muz_p':	z-component of dipole matrix elements
		'K_pp':		K-matrix
	'''
    # Kohn-Sham electron-hole pairs
	kss = lr.kss
	Np = len(kss)

	# Read arrays
	i_p = np.zeros(Np, dtype=int)
	a_p = np.zeros(Np, dtype=int)
	w_p = np.zeros(Np)
	f_p = np.zeros(Np)
	mu_vp = np.zeros((3, Np))
	for p, ks in enumerate(kss):
		i_p[p] = ks.i
		a_p[p] = ks.j
		w_p[p] = ks.energy
		f_p[p] = ks.fij
		mu_vp[:, p] = ks.mur

	# Read K matrix
	Omega_pp = lr.Om.full
	sqfw_p = np.sqrt(f_p * w_p)
	K_pp = 0.5 * (Omega_pp - np.diag(w_p**2)) / np.outer(sqfw_p, sqfw_p)

	np.savez_compressed(fpath, i_p=i_p, a_p=a_p, ediff_p=w_p, fdiff_p=f_p,
						mux_p=mu_vp[0], muy_p=mu_vp[1], muz_p=mu_vp[2],
						K_pp=K_pp)




c = GPAW('Na_gs_110bands.gpw')

dE = 6  # maximal Kohn-Sham transition energy to consider in eV
print("loading lr")
lr = LrTDDFT(c, filename = 'LrTDDFTresults.dat')
print("discrete_spectrum")
discrete_spectrum(lr, 'discrete_spectrum_GPAW.dat')

print("dump_data")
dump_data(lr,'dumpFile.npz')
#%%
dump = np.load('dumpFile.npz')
print("getting stuff")
K_pp = dump['K_pp']
omega_p= dump['ediff_p']
fdiff_p = dump['fdiff_p']

Omega=np.zeros(np.shape(K_pp))

for i in range(np.shape(K_pp)[0]):
    for j in range(np.shape(K_pp)[1]):
        Omega[i,j]=4*np.sqrt(omega_p[i]*omega_p[j])*K_pp[i,j]
        if i==j:
            Omega[i,j]+=omega_p[i]**2

[Omega2,F_I]=sp.linalg.eig(Omega)

# sort_perm = Omega2.argsort()

# Omega2 = Omega2[sort_perm]     # <-- This sorts the list in place.
# eig_vecs = F_I[:, sort_perm]

# Omegacalc=np.zeros(np.shape(K_pp)[0])
# for i in range(np.shape(K_pp)[0]):
#     Omegacalc += omega_p[i]**2*F_I[i,:]**2
#     for j in range(np.shape(K_pp)[1]):
#         Omegacalc+=F_I[i,:]*F_I[j,:]*K_pp[i,j]*4*np.sqrt(omega_p[i]*omega_p[j])

mu_x=dump['mux_p']
mu_y=dump['muy_p']
mu_z=dump['muz_p']

print("making f")
fx=2*np.abs(np.sum(mu_x*np.sqrt(2*omega_p)*F_I.T,axis=1))**2
fy=2*np.abs(np.sum(mu_y*np.sqrt(2*omega_p)*F_I.T,axis=1))**2
fz=2*np.abs(np.sum(mu_z*np.sqrt(2*omega_p)*F_I.T,axis=1))**2

f=(fx+fy+fz)/3
print("F ", f.shape)

x_t = np.linspace(np.sqrt(np.min(Omega2)*2*13.6),np.sqrt(np.max(Omega2))*2*13.6,200)
y_t = fold(x_t, np.sqrt(Omega2)*2*13.6, f, 0.06)

plt.plot(x_t,y_t)
plt.savefig('our_spectrum.pdf')

discrete_GPAW = np.loadtxt('discrete_spectrum_GPAW.dat')
x_t = np.linspace(np.min(discrete_GPAW[:,1]),np.max(discrete_GPAW[:,1]),200)
y_t = fold(x_t, discrete_GPAW[:,1], discrete_GPAW[:,2], 0.06)
plt.plot(x_t,y_t)
