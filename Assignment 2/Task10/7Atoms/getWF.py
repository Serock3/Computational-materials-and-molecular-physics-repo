#%%
from ase.db import connect
from ase.ga.data import DataConnection
from ase.optimize import GPMin
from ase import Atoms, Atom
from gpaw import GPAW, FermiDirac
from ase.io import read, write


db = connect("gadb.db")
atomsBest = db.get('id=72').toatoms()

# a=Atoms(read("guess.xyz", format='xyz'))


# Define the calculator
calc = GPAW(nbands=10,
            h=0.25,
            txt='out.txt',
            occupations=FermiDirac(0.05),
            setups={'Na': '1'},
            mode='lcao',
            basis='dzp')

atomsBest.set_calculator(calc)

print("Best post -a.get_potential_energy() ",-atomsBest.get_potential_energy())
calc.write('WFbest.gpw', mode='all')


# %%
