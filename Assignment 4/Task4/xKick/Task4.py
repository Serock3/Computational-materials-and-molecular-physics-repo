from gpaw.tddft import *

time_step = 30.0                  # 1 attoseconds = 0.041341 autime
iterations = 1500                # 1500 x 30 as => 45 fs
kick_strength = [1e-5,0.0,0.0]   # Kick to x-direction

# Read ground state
td_calc = TDDFT('Na_gs_0bands.gpw')

# Kick with a delta pulse to z-direction
td_calc.absorption_kick(kick_strength=kick_strength)

# Propagate, save the time-dependent dipole moment to 'na_dm.dat',
# and use 'na_td.gpw' as restart file
td_calc.propagate(time_step, iterations, 'na_dm.dat', 'na_td.gpw')


# Calculate photoabsorption spectrum and write it to 'be_spectrum_z.dat'
photoabsorption_spectrum('na_dm.dat', 'na_spectrum_x.dat',width=0.06)

# %%
