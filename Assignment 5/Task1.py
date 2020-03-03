#%%
from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM

#%% Task 1
for i in range (1,11):
    db = connect('Al-clusters-initial.db')
    atoms = db[i].toatoms()

    mishin = EAM(potential='HA5_al_potential.alloy')
    #mishin.write_potential('new.eam.alloy')
    from matplotlib import pyplot as plt
    #plt.figure(figsize=(8,8))
    #plt.tight_layout()
    #mishin.plot()
    atoms.set_calculator(mishin)
    Cohesive_energy=atoms.get_potential_energy()
    atoms.get_forces()
    size=atoms.get_number_of_atoms()
    print('Size:',size)
    print('Cohesive energy: ',Cohesive_energy/size)
    mean_lattice_param = np.mean([np.min(atoms.get_all_distances()[i,:][atoms.get_all_distances()[i,:]!=0]) for i in range(1,11)])
    # Lattice_param=np.mean(np.min(atoms.get_all_distances()[atoms.get_all_distances()!=0],axis=0))*np.sqrt(2)
    print('"Lattice param": ', mean_lattice_param*np.sqrt(2))

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

