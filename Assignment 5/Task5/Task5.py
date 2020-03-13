#%%
from ase.dft import get_distribution_moment
import matplotlib.pyplot as plt
from ase.dft.dos import DOS
from gpaw import GPAW, PW
db_rel_new = connect('../Task1/Al-clusters-relaxed-new.db')
db_elec = connect('../Task1/Al-clusters-elec.db')
calc = GPAW(xc='LDA',txt='ground.gpaw-out')

# for i in range(4, 5):

i = 1
atoms = db_rel_new[i].toatoms()
# if atoms.get_number_of_atoms() > 50:
#     continue
atoms.set_calculator(calc)
atoms.get_potential_energy()
    

calc.write('myCalc_'+str(int(atoms.get_number_of_atoms()))+'.gpw', 'all')

dos = DOS(calc, width=0.2)
d = dos.get_dos()
e = dos.get_energies()
db_elec.write(atoms,data={'energies':e,'dos':d})
plt.plot(e, d)
plt.xlabel('energy [eV]')
plt.ylabel('DOS')
plt.tight_layout()
plt.savefig('Task5_'+str(int(atoms.get_number_of_atoms()))+'.pdf'))
#plt.show()

volume = get_distribution_moment(e, d)
center, width = get_distribution_moment(e, d, (1, 2))


# %%
