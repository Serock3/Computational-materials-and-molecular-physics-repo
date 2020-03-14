#%%
from ase.dft import get_distribution_moment
import matplotlib.pyplot as plt
from ase.dft.dos import DOS
from gpaw import GPAW, PW
db_rel_new = connect('../Task1/Al-clusters-relaxed-new.db')
db_elec = connect('../Task1/Al-clusters-elec.db')
calc = GPAW(xc='LDA',txt='ground.gpaw-out')

# for i in range(4, 5):

for i in range(1,11):
    atoms = db_rel_new[i].toatoms()
    if atoms.get_number_of_atoms() > 99:
        continue
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
    plt.savefig('Task5_'+str(int(atoms.get_number_of_atoms()))+'.pdf')
    #plt.show()

    volume = get_distribution_moment(e, d)
    center, width = get_distribution_moment(e, d, (1, 2))


# %% Plot
sizes=np.zeros(6)
for i in range(1,7):
    sizes[i-1]=np.round(db_elec[i].toatoms().get_number_of_atoms())
    
perm=np.argsort(sizes)
sizes=np.round(sizes)
plt.figure(figsize=(8,6))
plt.plot(db_elec[2].data['energies'],db_elec[2].data['dos'],'-',label=str(int(sizes[1])))
plt.plot(db_elec[6].data['energies'],db_elec[6].data['dos'],':',linewidth=1.8,label=str(int(sizes[5])))
plt.plot(db_elec[5].data['energies'],db_elec[5].data['dos'],'-.',label=str(int(sizes[4])))
plt.plot(db_elec[4].data['energies'],db_elec[4].data['dos'],'-',label=str(int(sizes[3])))
plt.plot(db_elec[1].data['energies'],db_elec[1].data['dos'],':',linewidth=1.8,label=str(int(sizes[0])))
plt.plot(db_elec[3].data['energies'],db_elec[3].data['dos'],'g-',label=str(int(sizes[2])))
plt.tick_params(labelsize=12)
plt.legend(fontsize=12)
plt.xlabel('energy [eV]',fontsize=12)
plt.ylabel('DOS [1/eV]',fontsize=12)
plt.tight_layout()
plt.savefig('Task5.pdf')
# %%
