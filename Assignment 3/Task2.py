import numpy as np
from matplotlib import pyplot as plt
from ase.io.trajectory import Trajectory

#traj = Trajectory('./Na-aimd/NaCluster24.traj')
traj = Trajectory('someDynamics.traj')
start = 1000
stop = len(traj)
distances = np.array([])
for i in range(start,stop):
    atoms = traj[i]
    distances = np.append(distances,atoms.get_distances(72,atoms.get_atomic_numbers()==8,mic=True))

plt.hist(distances,bins = 100)
plt.show()

hist = np.histogram(distances,bins=120) # hist[0] is occurances, hist[1] is r

dn = hist[0]/(stop-start) # Number (density) of oxygen within [r, r+dr]
dr = hist[1][1]-hist[1][0] #The bins are evenly spaced
r = hist[1][1:]-dr/2

rho = 3.43e-2
gprime = dn/(4*np.pi*r**2*dr*rho)
g= dn/dr

guess = 3
width = 10
minapprox = int((guess-r[0])/dr) #Place to start searching for first minimum
minimum = np.argmin(gprime[minapprox-width:minapprox+width])+minapprox-width
plt.show()
plt.plot(r,gprime)
plt.figure(figsize=(8, 6))
plt.plot(r,gprime)
plt.xlabel(r"r [Å]", fontsize=18)
plt.ylabel(r"$g_{Na O}(r)$", fontsize=18)
plt.tight_layout()
plt.savefig('RDF_NaO_ours.pdf')
plt.show()
# print("solvation shell = ",np.trapz(g[:minimum],dx=dr))
print("solvation shell =  ", sum(hist[0][:minimum])/(stop-start))

