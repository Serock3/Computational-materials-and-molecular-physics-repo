#%%
from ase.db import connect
from ase.ga.data import DataConnection
from ase.optimize import GPMin
from ase import Atoms, Atom
from gpaw import GPAW, FermiDirac
from ase.io import read, write


db = connect("gadb.db")
atomsBest = db.get('id=41').toatoms()
atomsSecond = db.get('id=218').toatoms()

# Define the calculator
calc = GPAW(nbands=10,
            h=0.25,
            txt='out.txt',
            occupations=FermiDirac(0.05),
            setups={'Na': '1'},
            mode='pw',
            basis='dzp')

atomsBest.set_calculator(calc)
atomsSecond.set_calculator(calc)
print("Best pre -a.get_potential_energy() ",-atomsBest.get_potential_energy())
print("Second pre -a.get_potential_energy() ",-atomsSecond.get_potential_energy())

dyn = GPMin(atomsBest, trajectory='relaxBest_ref.traj', logfile='relaxBest_ref.log')
dyn.run(fmax=0.02, steps=100)
dyn2 = GPMin(atomsSecond, trajectory='relaxSecond_ref.traj', logfile='relaxSecond_ref.log')
dyn2.run(fmax=0.02, steps=100)
print("Best post -a.get_potential_energy() ",-atomsBest.get_potential_energy())
print("Second post -a.get_potential_energy() ",-atomsSecond.get_potential_energy())
print("pw done")