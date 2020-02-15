import time
start_time = time.time()

from random import random
from ase.io import write
from ase.io import read
from ase.db import connect
from ase.optimize import GPMin
from gpaw import GPAW, FermiDirac

from ase.ga.data import DataConnection
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
import argparse

from ase.md.npt import NPT
atoms = read('someStartConfiguration.xyz')

calc = GPAW(
    mode     = 'lcao',
    xc       = 'PBE',
    basis    = 'dzp',
    symmetry= {'point -group ': False}, # Turn  off point -group  symmetry
    charge   = 1, # Charged  system
    txt      = 'gpawOutput.gpaw -out', # Redirects  calculator  output  to thisfile!
    )
atoms.set_calculator(calc)

from ase.units import fs, kB

dyn = NPT( # Some MD  method
    temperature = 350*kB,
    timestep = 1000*fs, # This is not an  appropriate  timestep , I can  tell  youthat!
    ttime = 20*fs, # Don't forget  the fs!
    logfile = 'mdOutput.log', # Outputs  temperature (and  more) to file at eachtimestep
    )
trajectory = Trajectory('someDynamics.traj', 'w', atoms)
dyn.attach(trajectory.write , interval =1) # Write  the  current  positions  etc. to fileeach  timestep

dyn.run (10) # Run 10  steps of MD  simulation