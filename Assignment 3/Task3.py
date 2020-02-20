#%%
#Save wavefcn from gadb.db in .gpw format
from ase.io.trajectory import Trajectory
from ase import Atoms, Atom
from matplotlib import pyplot as plt
import numpy as np
traj = Trajectory('./Na-aimd/cluster24.traj')
# traj = Trajectory('someDynamics.traj')
start = 4000
stop = len(traj)
# stop = 1001

distances = np.array([])
for atomToCalc in range(24):
# atomToCalc=5
    
    for i in range(start,stop):
        atoms = traj[i]
        distToCalc=atoms.get_atomic_numbers()==8
        distToCalc[atomToCalc]=False
        distances = np.append(distances,atoms.get_distances(atomToCalc,distToCalc,mic=True))

plt.hist(distances,bins = 100)
plt.show()

hist = np.histogram(distances,bins=120) # hist[0] is occurances, hist[1] is r

dn = hist[0]/(stop-start)/24 # Number (density) of oxygen within [r, r+dr]
dr = hist[1][1]-hist[1][0] #The bins are evenly spaced
r = hist[1][1:]-dr/2

rho = 3.43e-2
gprime = dn/(4*np.pi*r**2*dr*rho)
g= dn/dr

guess = 3.2
width = 4
minapprox = int((guess-r[0])/dr) #Place to start searching for first minimum
minimum = np.argmin(gprime[minapprox-width:minapprox+width])+minapprox-width
plt.show()
plt.figure(figsize=(8, 6))
plt.plot(r,gprime)
plt.xlabel(r"r [Å]", fontsize=18)
plt.ylabel(r"$g_{OO}(r)$", fontsize=18)
plt.tight_layout()
plt.savefig('RDF_water.pdf')
plt.show()

print("solvation shell = ",np.trapz(g[:minimum],dx=dr))
print("asdasd ", sum(hist[0][:minimum])/(stop-start)/24)

# %%
