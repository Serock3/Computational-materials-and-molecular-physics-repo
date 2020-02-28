#%%
from matplotlib import pyplot as plt
import numpy as np
from gpaw.tddft import *

e_min = 1
e_max = 6.5
delta_e = 0.03
photoabsorption_spectrum('na_dm_x.dat', 'na_spectrum_x.dat',width=0.06,e_min=e_min,e_max=e_max,delta_e=delta_e)
photoabsorption_spectrum('na_dm_y.dat', 'na_spectrum_y.dat',width=0.06,e_min=e_min,e_max=e_max,delta_e=delta_e)
photoabsorption_spectrum('na_dm_z.dat', 'na_spectrum_z.dat',width=0.06,e_min=e_min,e_max=e_max,delta_e=delta_e)
#%%
fx=np.loadtxt('../Task4/na_spectrum_x.dat')
fy=np.loadtxt('../Task4/na_spectrum_y.dat')
fz=np.loadtxt('../Task4/na_spectrum_z.dat')

f=(fx[:,1]+fy[:,2]+fz[:,3])/3

plt.plot(fx[:,0],f,label='Time-propagation spectrum')
# plt.xlim(0,6)
filename1 = '../Task1/spectrum_w.06eV.dat'
# import data
print("Plotting")
data = np.loadtxt(filename1)
plt.plot(data[:,0],data[:,1],'--', label='LrTDDFT')
plt.xlabel('Transition energy [eV]')
plt.ylabel('Oscillation strength [a.u.]')
plt.legend()

plt.savefig('../Task4/Task4.pdf')

# %%
