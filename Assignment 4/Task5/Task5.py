from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum

lr = LrTDDFT('LrTDDFTresults.dat')
lr.diagonalize(energy_range =4)

print("Writing the spectrum to file")
photoabsorption_spectrum(lr, 'spectrum_w.04eV.dat', # data file name
                         width=0.06)                # width in eV

lr.write('LrTDDFTresults_4eV.dat')