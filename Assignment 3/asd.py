import argparse
import time
from random import random

from ase.db import connect
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.standardmutations import (MirrorMutation, PermutationMutation,
                                      RattleMutation)
from ase.ga.utilities import closest_distances_generator, get_all_atom_types
from ase.io import Trajectory, read, write
from ase.md.npt import NPT
from ase.optimize import GPMin
from ase.units import fs, kB
from gpaw import GPAW, FermiDirac

start_time = time.time()




atoms = read('NaIncert.xyz')

calc = GPAW(
    mode     = 'lcao',
    xc       = 'PBE',
    basis    = 'dzp',
    symmetry= {'point_group': False}, # Turn  off point-group  symmetry
    charge   = 1, # Charged  system
    txt      = 'gpawOutput.gpaw-out', # Redirects  calculator  output  to this file!
    )
atoms.set_calculator(calc)


dyn = NPT( # The constant NPT method
    atoms,
    timestep = 0.5*fs, # This is an  appropriate  timestep , I can  tell  you that!
    temperature = 350*kB,
    externalstress = 0,
    ttime = 20*fs, 
    pfactor = None,
    logfile = 'mdOutput.log', # Outputs  temperature (and  more) to file at eachtimestep
    )
trajectory = Trajectory('someDynamics.traj', 'w', atoms)
dyn.attach(trajectory.write , interval =1) # Write  the  current  positions  etc. to file each  timestep

dyn.run (4000) # Run 4000  steps of MD  simulation
