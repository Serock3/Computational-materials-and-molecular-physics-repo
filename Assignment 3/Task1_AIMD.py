from ase.io import Trajectory, read
from ase.md.npt import NPT
from ase.units import fs, kB
from gpaw import GPAW

atoms = read('NaIncert.xyz')

calc = GPAW(
    mode='lcao',
    xc='PBE',
    basis='dzp',
    symmetry={'point_group': False},  # Turn  off point -group  symmetry
    charge=1,  # Charged  system
    txt='gpawOutput.gpaw-out',  # Redirects  calculator  output  to thisfile!
)
atoms.set_calculator(calc)


dyn = NPT(  # Some MD  method
    atoms,
    timestep=0.5*fs,  # This is not an  appropriate  timestep , I can  tell  you that!
    temperature=350*kB,
    externalstress=0,
    ttime=20*fs,  # Don't forget  the fs!
    pfactor=None,
    # Outputs  temperature (and  more) to file at eachtimestep
    logfile='mdOutput.log',
)
trajectory = Trajectory('someDynamics.traj', 'w', atoms)
# Write  the  current  positions  etc. to file each  timestep
dyn.attach(trajectory.write, interval=1)

dyn.run(4000)  # Run 10  steps of MD  simulation
