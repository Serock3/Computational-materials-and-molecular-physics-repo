
#%%

from ase.db import connect
from ase.ga.data import DataConnection
from ase.optimize import GPMin
from ase import Atoms, Atom
from gpaw import GPAW, FermiDirac
from ase.io import read

db = connect('guess.db')
da = DataConnection('guess.db')

a=Atoms(read("guess.xyz", format='xyz'))

# Define the calculator
calc = GPAW(nbands=10,
            h=0.25,
            txt='out.txt',
            occupations=FermiDirac(0.05),
            setups={'Na': '1'},
            mode='lcao',
            basis='dzp')

a.set_calculator(calc)
#print('Relaxing starting candidate {0}'.format(a.info['confid']))
dyn = GPMin(a, trajectory='relax_ref.traj', logfile='relax_ref.log')
dyn.run(fmax=0.02, steps=2)
# %%
