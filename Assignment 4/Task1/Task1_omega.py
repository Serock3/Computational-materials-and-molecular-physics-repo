#%%
from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum
from matplotlib import pyplot as plt
import numpy as np
c = GPAW('Na_gs_110bands.gpw')

dE = 6  # maximal Kohn-Sham transition energy to consider in eV
print("Instantiating lr")
lr = LrTDDFT(c, xc='LDA', energy_range=dE)
print("Diagonalizing lr")
lr.diagonalize()
# write the spectrum to the data file
print("Writing the spectrum to file")
photoabsorption_spectrum(lr, 'spectrum_w.06eV.dat', # data file name
                         width=0.06)                # width in eV

lr.write('LrTDDFTresults.dat')
#%%
filename1 = '../Task1/spectrum_w.06eV.dat'
# import data
print("Plotting")
data = np.loadtxt(filename1)
plt.plot(data[:,0],data[:,1])
plt.xlabel('a.u.')
plt.ylabel('eV')
plt.savefig('test.pdf')


# %%
