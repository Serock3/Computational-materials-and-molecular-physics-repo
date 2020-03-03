#%%
from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM
import numpy as np
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
# %%
from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM
import numpy as np
from ase.vibrations import Vibrations
import os.path
from os import path

db = connect('../Task3/Al-clusters-initial.db')

new_db=connect('../Task3/Vibrations.db')

mishin = EAM(potential='../Task3/HA5_al_potential.alloy')
for i in range (1,11):
    
    atoms = db[i].toatoms()
    atoms.set_calculator(mishin)
    vibs=Vibrations(atoms,name='vib'+str(i))
    if os.path.exists("vib"+str(i)+".pckl"):
        vibs.split()
    vibs.run()
    vibs.combine()

    freq=vibs.get_frequencies()
    vibs.write_dos('dos.txt',start=0)
    dos=np.loadtxt('dos.txt')

    new_db.write(atoms, data = {'frequency':freq, 'density-of-states':dos})
    print(new_db[i].data['frequency'])
    print(new_db[i].data['density-of-states'])

#%% print freqs and dos
for i in range (1,11):
    print(new_db[i].data['frequency'])
    print(new_db[i].data['density-of-states'])