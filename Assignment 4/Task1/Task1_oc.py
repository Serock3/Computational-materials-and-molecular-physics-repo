from ase import Atom, atoms
from ase.io import Trajectory, read
from gpaw import GPAW

calc = GPAW('Na_gs_0bands.gpw')
calc.get_eigenvalues()  # Gives  the  eigenvalues  of the KS  equation!
# Gives  occupation  numbers  of  states  included  incalculation!
calc.get_occupation_numbers()
atoms = calc.get_atoms()

calc.set(
    nbands=110,    # Changes number of bands included
    fixdensity=True,  # Fixes density in any new calculator
)

# converge also the empty states (the density is converged already)
calc.set(convergence={'bands': -10},
         fixdensity=True,
         eigensolver='cg')
atoms.get_potential_energy()

print("Eigenvalues: ", calc.get_eigenvalues())
# Save the ground state
calc.write('Na_gs_110bands.gpw', 'all')

