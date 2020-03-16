#%%
from ase.build import bulk
from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM
import numpy as np
from matplotlib import pyplot as plt
from gpaw import GPAW, PW

atoms = bulk('Al', 'fcc', a=a_best)

energy = np.zeros(15)

for i in range(15):
    calc = GPAW(mode=PW(300), kpts=(i+1, i+1, i+1))

    atoms.set_calculator(calc)
    energy[i]=atoms.get_potential_energy()

# %%
plt.plot(range(1,16),energy)
plt.xlabel('Grid points per dimension')
plt.ylabel('Energy [eV]')
plt.savefig('Task6.pdf')
# %%

calc = GPAW(mode=PW(50), kpts=(10, 10, 10))

atoms.set_calculator(calc)
energy[i]=atoms.get_potential_energy()

calc.write('al_bulk.gpw')


#%%
from ase.dft.dos import DOS
from gpaw import restart
# from ase.calculators.test import FreeElectrons


atoms, calc = restart('al_bulk.gpw')

# calc = GPAW('groundstate.gpw')
# kpts={'size': (10, 10, 10)}
kpts=(10, 10, 10)
calc.set(kpts = kpts , fixdensity = True)
atoms.set_calculator(calc)
dos = DOS(calc, npts=500, width=0)
energies = dos.get_energies()
weights = dos.get_dos()
plt.figure(figsize=(8,6))
plt.plot(energies, weights,label='Al bulk')

m_e=5.110e5     #eV/c^2
hbar=6.582e-16  #eVs
c_0=3e18        #Å/s
k_bT=25e-3      #eV

freeDOS=atoms.get_volume()/(2*np.pi**2)*(2*m_e/hbar**2/c_0**2)**(3/2)*np.sqrt(energies-energies[0])
fermi_dirac=1/(np.exp((energies-calc.get_fermi_level())/k_bT)+1)
plt.plot(energies,freeDOS*fermi_dirac,'--',label='Free electron gas')
plt.legend()
plt.show()

n=(calc.get_fermi_level()*2*m_e/c_0**2/hbar**2)**(3/2)/(3*np.pi**2)


# %%
