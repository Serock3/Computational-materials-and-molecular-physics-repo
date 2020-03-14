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
atoms, calc = restart('al_bulk.gpw')

# calc = GPAW('groundstate.gpw')
# kpts={'size': (10, 10, 10)}
kpts=(10, 10, 10)
calc.set(kpts = kpts , fixdensity = True)
atoms.set_calculator(calc)
dos = DOS(calc, npts=500, width=0)
energies = dos.get_energies()
weights = dos.get_dos()
plt.plot(energies, weights)
plt.show()



# %%
