#%%
from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM
from ase.optimize import GPMin

import numpy as np
from matplotlib import pyplot as plt

db_rel = connect('Task1/Al-clusters-relaxed.db')

sizes  = np.zeros(10)
E_coh = np.zeros(10)
lattice_param = np.zeros(10)
for i in range(1,11):
    if i==9: 
        continue
    E_coh[i-1] = db_rel[i].data['Cohesive-energy']
    sizes[i-1] = db_rel[i].toatoms().get_number_of_atoms()
    lattice_param[i-1] = db_rel[i].data['Lattice-param']
E_coh[9-1] = db_rel[13].data['Cohesive-energy']
sizes[9-1] = db_rel[13].toatoms().get_number_of_atoms()
lattice_param[9-1] = db_rel[13].data['Lattice-param']
# plt.figure(figsize=(8,6))

plt.scatter(sizes,E_coh)
plt.xlabel('Number of atoms')
plt.ylabel('Cohesive energy per atom [eV]')

plt.tight_layout()
plt.savefig('Task1.pdf')
