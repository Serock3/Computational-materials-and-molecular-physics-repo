#%%
from ase.build import bulk
from ase.phonons import Phonons
from gpaw import GPAW, FermiDirac, PW
atoms = bulk('Si', 'diamond', a=5.4)
calc = GPAW(mode=PW(300), kpts=(7, 7, 7),h=0.2,occupations=FermiDirac(0.),symmetry='off')
ph = Phonons(atoms, calc, supercell=(2, 2, 2))
#%%
ph.run()
#%%
ph.read(method='Frederiksen', acoustic=True)
#ph.clean()

path = atoms.cell.bandpath('GXKGLXWL', npoints=100)
bs = ph.get_band_structure(path)

dos = ph.get_dos(kpts=(7, 7, 7)).sample_grid(npts=100, width=1e-3)
# %%
import matplotlib.pyplot as plt
fig = plt.figure(1, figsize=(8, 6))
ax = fig.add_axes([.12, .07, .67, .85])

emax = 0.07
PLOT=bs.plot(ax=ax, emin=-0.01, emax=emax)



# %%
si = bulk('Si', 'diamond', 5.43)
calc = GPAW(mode=PW(300),
            xc='PBE',
            kpts=(8, 8, 8),
            random=True,  # random guess (needed if many empty bands required)
            occupations=FermiDirac(0.01),
            txt='Si_gs.txt')
si.calc = calc
si.get_potential_energy()
calc.write('Si_gs.gpw')

#%%
# Restart from ground state and fix potential:
calc = GPAW('Si_gs.gpw',
            nbands=16,
            fixdensity=True,
            symmetry='off',
            kpts={'path': 'GXWKL', 'npoints': 60},
            convergence={'bands': 8})

calc.get_potential_energy()
# %%

bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=True, emin=-5 ,emax=15.5)

# %%
from ase.dft.bandgap import bandgap 
gap, p1, p2 = bandgap(calc)
print(gap, p1, p2)
gap, p1, p2 = bandgap(calc, direct=True)
print(gap, p1, p2)
# %%
