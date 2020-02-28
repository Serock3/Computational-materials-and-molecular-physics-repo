from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum

lr = LrTDDFT('LrTDDFTresults.dat')
lr.analyse(range(10))
lr.diagonalize(energy_range =4)
lr.analyse(range(10))
print("Writing the spectrum to file")
photoabsorption_spectrum(lr, 'spectrum_w.04eV.dat', # data file name
                         width=0.06)                # width in eV

lr.write('LrTDDFTresults_4eV.dat')

#%%
import numpy as np
from matplotlib import pyplot as plt
filename1 = '../Task5/spectrum_w.04eV.dat'
# import data
print("Plotting")
data = np.loadtxt(filename1)
plt.plot(data[:,0],data[:,1],label='4 eV cutoff')
plt.xlabel('eV')
plt.ylabel('counts?')
# plt.savefig('test.pdf')
filename1 = '../Task1/spectrum_w.06eV.dat'
# import data
print("Plotting")
data = np.loadtxt(filename1)
plt.plot(data[:,0],data[:,1],'--',label='6 eV cutoff')
plt.xlabel('Transition energy [eV]')
plt.ylabel('Oscillation strength [a.u.]')
plt.legend()

plt.savefig('../Task5/Task5.pdf')


# %%


# %%
