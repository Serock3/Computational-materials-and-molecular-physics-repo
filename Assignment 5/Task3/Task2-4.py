#%%
from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM
import numpy as np
from matplotlib import pyplot as plt
#%% Task 1

for i in range (1,11):
    db = connect('Al-clusters-initial.db')
    atoms = db[i].toatoms()

    mishin = EAM(potential='HA5_al_potential.alloy')
    #mishin.write_potential('new.eam.alloy')
    
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
plt.xlabel('Lattice parameter [Å]')
plt.ylabel('Energy per atom [eV]')
a_best=a[np.argmin(E)]
E_ground=E[np.argmin(E)]

plt.tight_layout()
plt.savefig('Task2.pdf')
print('GS Energy: ', E_ground, 'at a=', a_best)
# %% Task3
from ase.db import connect
from ase import Atom, atoms
from ase.io import Trajectory, read
from ase.calculators.eam import EAM
import numpy as np
from ase.vibrations import Vibrations
import os.path
from os import path

db = connect('./Task3/Al-clusters-initial.db')

new_db=connect('./Task3/Vibrations.db')

db_rel_new = connect('./Task1/Al-clusters-relaxed-new.db')

mishin = EAM(potential='./Task3/HA5_al_potential.alloy')
for i in range (1,11):
    
    atoms = db_rel_new[i].toatoms()
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

#%% print freqs and dos0
new_db=connect('./Task3/Vibrations.db')
plt.figure(figsize=(8,6))
plt.xlim(0,350)
sizes=np.zeros(10)
for i in range(1,11):
    #print(new_db[i].data['frequency'])
    #print(new_db[i].data['density-of-states'])
    #dos=new_db[i].data['density-of-states']
    sizes[i-1]=np.round(new_db[i].toatoms().get_number_of_atoms())
    #plt.plot(dos[:,0],dos[:,1],label=str(sizes[i-1]))
    
perm=np.argsort(sizes)
sizes=np.round(sizes)

dos=new_db[9+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[9])))

dos=new_db[8+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[8])))

dos=new_db[4+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[4])))

dos=new_db[2+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[2])))

dos=new_db[1+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[1])))

dos=new_db[7+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[7])))

dos=new_db[6+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[6])))

dos=new_db[5+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[5])))

dos=new_db[0+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(sizes[0])))

dos=new_db[3+1].data['density-of-states']
plt.plot(dos[:,0],dos[:,1],label=str(int(int(sizes[3]))))


plt.tick_params(labelsize=12)
plt.legend(fontsize=12)
plt.xlabel('E [cm$^{-1}$]',fontsize=12)
plt.ylabel('DOS [1/cm$^{-1}$]',fontsize=12)
plt.tight_layout()

plt.savefig('DOS.pdf')
#%% Task 3 freq count

freq = [None]*10
freq[9]=new_db[9+1].data['frequency'].shape[0]
freq[8]=new_db[8+1].data['frequency'].shape[0]
freq[4]=new_db[4+1].data['frequency'].shape[0]
freq[2]=new_db[2+1].data['frequency'].shape[0]
freq[1]=new_db[1+1].data['frequency'].shape[0]
freq[7]=new_db[7+1].data['frequency'].shape[0]
freq[6]=new_db[6+1].data['frequency'].shape[0]
freq[5]=new_db[5+1].data['frequency'].shape[0]
freq[0]=new_db[0+1].data['frequency'].shape[0]
freq[3]=new_db[3+1].data['frequency'].shape[0]
plt.scatter(sizes,freq)

plt.tick_params(labelsize=12)
plt.ylabel('Number of modes',fontsize=12)
plt.xlabel('Cluster size',fontsize=12)
plt.tight_layout()
plt.savefig('freqs.pdf')
# %% Task 4

from ase.build import bulk
from ase.calculators.emt import EMT
from ase.phonons import Phonons

# Setup crystal and EMT calculator
atoms = bulk('Al', 'fcc', a=a_best)

# Phonon calculator
N = 7
ph = Phonons(atoms, EMT(), supercell=(N, N, N), delta=0.05)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)
ph.clean()

path = atoms.cell.bandpath('GXWKGLUWLK', npoints=100)
bs = ph.get_band_structure(path)

dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)
#%%
# Plot the band structure and DOS:
import matplotlib.pyplot as plt
fig = plt.figure(1, figsize=(7, 4))
ax = fig.add_axes([.12, .07, .67, .85])

emax = 0.035
PLOT=bs.plot(ax=ax, emin=0.0, emax=emax)



dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(dos.weights[0], dos.energy, y2=0, color='grey',
                   edgecolor='k', lw=1)

dosax.set_ylim(0, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=18)

fig.savefig('Al_phonon.png')

# %%
plt.plot(dos.energy,dos.weights[0])

# %%

# %%
