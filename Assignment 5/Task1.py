#%%
from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM
from ase.optimize import GPMin

import numpy as np
from matplotlib import pyplot as plt
#%% Task 1

sizes  = np.zeros(10)
E_coh = np.zeros(10)

db = connect('Al-clusters-initial.db')
db_rel = connect('Al-clusters-relaxed.db')
mishin = EAM(potential='HA5_al_potential.alloy')
for i in range (1,11):
    
    atoms = db[i].toatoms()

    
    #mishin.write_potential('new.eam.alloy')
    
    #plt.figure(figsize=(8,8))
    #plt.tight_layout()
    #mishin.plot()
    atoms.set_calculator(mishin)

    dyn = GPMin(atoms, trajectory='relax_ref.traj', logfile='relax_ref.log')
    dyn.run(fmax=0.02, steps=100)

    Cohesive_energy=atoms.get_potential_energy()
    atoms.get_forces()

    
    size=atoms.get_number_of_atoms()
    print('Size:',size)
    sizes[i-1] = size
    print('Cohesive energy: ',Cohesive_energy/size)
    E_coh[i-1] = Cohesive_energy/size
    mean_lattice_param = np.mean([np.min(atoms.get_all_distances()[i,:][atoms.get_all_distances()[i,:]!=0]) for i in range(1,11)])
    # Lattice_param=np.mean(np.min(atoms.get_all_distances()[atoms.get_all_distances()!=0],axis=0))*np.sqrt(2)
    print('"Lattice param": ', mean_lattice_param*np.sqrt(2))
    db_rel[i].write(atoms,data={'Cohesive-energy': Cohesive_energy/size, 'Lattice-param':mean_lattice_param*np.sqrt(2)})

plt.figure(figsize=(8,6))

plt.scatter(sizes,E_coh)
plt.xlabel('Number of atoms')
plt.ylabel('Cohesive energy per atom [eV]')

plt.tight_layout()
plt.savefig('Task1.pdf')

# %% Task 2

from ase.build import bulk
import numpy as np

n = 201
E = np.zeros(n)
a = np.linspace(3.9,4.2,n)
print('accuracy in a:', a[1]-a[0])
for i in range(n):
    # a = 4.05  # Angstrom lattice spacing
    al = bulk('Al', 'fcc', a=a[i])
    al.set_calculator(mishin)
    
    E[i] = al.get_potential_energy()

plt.plot(a,E)
a_best=a[np.argmin(E)]
E_ground=E[np.argmin(E)]

print('GS Energy: ', E_ground, 'at a=', a_best)
# %% Old task 3

from ase.vibrations import Vibrations
db = connect('Al-clusters-initial.db')
mishin = EAM(potential='HA5_al_potential.alloy')
for i in range (1,11):
    
    atoms = db[i].toatoms()
    al.set_calculator(mishin)
    vibs=Vibrations(atoms)

    freq=vibs.get_frequencies()
    vibs.write_dos('dos.txt')
    dos=np.loadtxt('dos.txt')

    db[i].write(atoms, data = {'frequency' : freq , 'density -of -states' : dos})

# %%
1+1

# %%
