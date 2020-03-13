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
for i in range(1,11):
    
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

# %% Plot
from matplotlib import pyplot as plt
i=1
plt.plot(new_db[i].data['density-of-states'][:,0],new_db[i].data['density-of-states'][:,1])
plt.xlim(0,400)
plt.figure()
plt.plot(new_db[i].data['frequency'])
# %%
