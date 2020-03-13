from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM
from ase.optimize import GPMin, MDMin


import numpy as np
from matplotlib import pyplot as plt
#%% Task 1

db = connect('Al-clusters-initial.db')
db_rel = connect('Al-clusters-relaxed.db')
mishin = EAM(potential='HA5_al_potential.alloy')

i=9
print("Started nanoparticle ", i)
atoms = db[i].toatoms()

atoms.set_calculator(mishin)

# dyn = GPMin(atoms, trajectory='relax_ref_'+str(i)+'.traj', logfile='relax_ref_'+str(i)+'.log')
dyn = MDMin(atoms, trajectory='relax_ref_'+str(i)+'.traj', logfile='relax_ref_'+str(i)+'.log')
dyn.run(fmax=0.02, steps=100)

Cohesive_energy=atoms.get_potential_energy()
atoms.get_forces()


size=atoms.get_number_of_atoms()
print('Size:',size)
print('Cohesive energy: ',Cohesive_energy/size)
mean_lattice_param = np.mean([np.min(atoms.get_all_distances()[i,:][atoms.get_all_distances()[i,:]!=0]) for i in range(1,11)])
print('"Lattice param": ', mean_lattice_param*np.sqrt(2))
db_rel.write(atoms,data={'Cohesive-energy': Cohesive_energy/size, 'Lattice-param':mean_lattice_param*np.sqrt(2)})


# %%

db_rel_new = connect('Task1/Al-clusters-relaxed-new.db')
for i in range(1,11):
    if i==9:
        continue
    atoms = db_rel[i].toatoms()
    db_rel_new.write(atoms,data={'Cohesive-energy': db_rel[i].data['Cohesive-energy'], 'Lattice-param':db_rel[i].data['Lattice-param']})
atoms = db_rel[13].toatoms()
db_rel_new.write(atoms,data={'Cohesive-energy': db_rel[13].data['Cohesive-energy'], 'Lattice-param':db_rel[13].data['Lattice-param']})

# %%
