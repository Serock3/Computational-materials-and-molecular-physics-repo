#%%
from ase import Atom, atoms
from ase.io import Trajectory, read
from gpaw import GPAW

# %% Task 1 only occupied bands
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


# %%
