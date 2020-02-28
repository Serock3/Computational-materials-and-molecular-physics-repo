#%%
from ase import Atom, atoms
from ase.io import Trajectory, read
from gpaw import GPAW

# %% Task 1 occupied
atoms = read('our_cluster.xyz')

atoms.center(vacuum=8)  # Apply vacuum

calc = GPAW(mode='fd',
            xc='LDA',
            setups = {'Na':'1'},
            h = 0.3,
            nbands=0,
            txt='ground.gpaw-out',  # This  redirects  the  text  output  to this file!
            )

# Attach calculator to atoms
atoms.set_calculator(calc)

# Calculate the ground state
energy = atoms.get_potential_energy()
print("Energy without empty states: ", energy)
# Save the ground state
calc.write('Na_gs_0bands.gpw', 'all')

#%%

calc = GPAW('myCalc.gpw')
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

# Save the ground state
calc.write('Na_gs_110bands.gpw', 'all')

# %% Save
atoms = read('someAtomsFromSomewhere.xyz')

calc = GPAW(xc='LDA',
            nbands=4,     # Includes 4 bands  in  calculation!
            ...
            txt='outputFile.gpaw -out',  # This  redirects  the  text  output  to this file!
            )
atoms.set_calculator(calc)
atoms.get_potential_energy()  # Performs  the  calculation  on the  system!
calc.write('myCalc.gpw', 'all')  # Writes  everything  it can to file!

# %% Load
calc = GPAW('myCalc.gpw')
calc.get_eigenvalues()  # Gives  the  eigenvalues  of the KS  equation!
# Gives  occupation  numbers  of  states  included  incalculation!
calc.get_occupation_numbers()
atoms = calc.get_atoms()
atoms.get_potential_energy()  # Get  energy  without  recalculating  again!
# You can  also  change  the  settings  from a previous  calculation
# to  continue  working  from  the  same  place:

calc.set(
    nbands=10,    # Changes number of bands included
    fixdensity=True,  # Fixes density in any new calculator
)

atoms.get_potential_energy()  # Recalculate with new parameters

# %%
