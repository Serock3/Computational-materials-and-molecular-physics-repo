from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import helper
c = GPAW('Na_gs_110bands.gpw')

dE = 6  # maximal Kohn-Sham transition energy to consider in eV
lr = LrTDDFT(c, filename = 'LrTDDFTresults.dat')
helper.discrete_spectrum(lr, 'discrete_spectrum_GPAW.dat')

helper.dump_data(lr,'dumpFile.npz')
dump = numpy.load('dumpFile.npz')
K_pp = dump['K_pp']
omega_p= dump['ediff_p']

Omega=np.zeros(np.shape(K_pp))

for i in range(np.shape(K_pp)[0]):
    for j in rnage(np.shape(K_pp)[1]):
        Omega[i,j]=4*np.sqrt(omega_p[i]*omega_p[j])*K_pp[i,j]
        if i==j:
            Omega[i,j]+=omega_p[i]**2

[Omega2,F_I]=sp.linalg.eig(Omega)

mu_x=dump['mux_p']
mu_y=dump['muy_p']
mu_z=dump['muz_p']

fx=2*np.abs(np.sum(mu_x*np.sqrt(2*omega_p)*F_I))**2

x_t = np.linspace(np.sqrt(np.min(Omega2))-1,np.sqrt(np.max(Omega2))+1)
y_t = helper.fold(x_t, np.sqrt(Omega2), fx, 0.06)