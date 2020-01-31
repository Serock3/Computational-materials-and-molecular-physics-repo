# %%
import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import scipy
sns.set('talk')

#Get Hartree potential from electron density
def getVH(ns,N,Rmax):
    h = Rmax/N
    A = np.diag(-2*np.ones(N),0)+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1)
    r=np.linspace(h,Rmax,N)
    tmp = -4*np.pi*h**2*r*ns
    tmp[-1] -= 1
    U = np.linalg.solve(A,tmp)
    return U/r

#Get wave function from hartree potential
def getpsi(V,N,Rmax,Z):
    h = Rmax/N
    r = np.linspace(h,Rmax,N)
    c = 1/h**2-Z/r+V
    A = np.diag(c,0)+np.diag(-np.ones(N-1)/(2*h**2),1)+np.diag(-np.ones(N-1)/(2*h**2),-1)
    
    (E, w) = np.linalg.eig(A)
    eps = np.min(E).real
    f = w[:,np.argmin(E)]
    f = f/np.sqrt(np.trapz(f**2,r))*np.sign(f[0])
    psi = 1/np.sqrt(4*np.pi)*f/r
    return (psi,eps)

# Calculate energy
def getE(eps,VH,Vxc,epsxc,ns,Z):
    return Z*eps-Z*4*np.pi*np.trapz((VH*ns/2+Vxc*ns-epsxc*ns)*r**2,r)


#%% task 2

Rmax = 10 # atomic
Z=1
N=1000                      #Number of grid points
h = Rmax/N                  #Stepsize
r = np.linspace(h,Rmax,N)   
ns=1/np.pi*Z**3*np.exp(-2*Z*r)  #Hydrogen density
VH=getVH(ns,N,Rmax)             #Get Hydrogen Hartree potential
VHanalytic = 1/r-(1+1/r)*np.exp(-2*r)


#Plot
plt.figure(figsize=(8, 6))
plt.plot(r,VH, label = r'Hartree method')
plt.plot(r,VHanalytic, '--' ,label = r'Analytical')
plt.legend()

plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Potential [Ha]", fontsize=18)
plt.savefig('task2.pdf')
#%% task 3

Rmax = 6 # atomic
Z=1
N=1000                  #Number of grid points
h = Rmax/N              #Stepsize
r = np.linspace(h,Rmax,N)

(psi,eps) = getpsi(0,N,Rmax,Z)      #Get hydrogen wave function
print("eps = ", eps)
ns= np.abs(psi)**2
print("E = ",getE(eps,0,0,0,ns,Z))  #Get hydrogen ground state energy

#Plot
plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Radial PDF", fontsize=18)
plt.plot(r,4*np.pi*r**2*ns)
plt.savefig('task3.pdf')
#%% Task 4
Z= 2
N=800               #Number of grid points
Rmax = 6            # atomic
h = Rmax/N          #Stepsize
r = np.linspace(h,Rmax,N)
ns=1/np.pi*Z**3*np.exp(-2*Z*r)      #Guess initial density

MaxIters = 15           #Max number of iterations
E = np.zeros(MaxIters)  

#Solve the self consistency problem
for i in range(MaxIters):
    VH = getVH(ns,N,Rmax)
    (psi,eps) = getpsi(VH,N,Rmax,Z)
    ns = np.abs(psi)**2
    E[i] =getE(eps,VH,0,0,ns,Z)
    if np.abs(E[i]-E[i-1])<1e-5/27.21:
        break
print("eps = ", eps)
print("E = ", E[i])

#Plot
plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Radial PDF", fontsize=18)
plt.plot(r,4*np.pi*r**2*ns)
plt.savefig('task4.pdf')
#%% Task 5
Z= 2
N=1000          #Number of grid points
Rmax = 6        # atomic
h = Rmax/N      #Stepsize
r = np.linspace(h,Rmax,N)

ns=1/np.pi*Z**3*np.exp(-2*Z*r)      #Guess initial density

#Exchange potential
epsx= - 3/4*(3*ns/np.pi)**(1/3)
ndepsx= 1/3*epsx
Vx = epsx + ndepsx

MaxIters = 30           #Max number of iterations
E = np.zeros(MaxIters)

#Solve the self consistency problem
for i in range(MaxIters):
    VH = 2*getVH(ns,N,Rmax)
    (psi,eps) = getpsi(VH+Vx,N,Rmax,Z)
    ns = np.abs(psi)**2
    E[i] =getE(eps,VH,Vx,epsx,ns,Z)
    if np.abs(E[i]-E[i-1])<1e-5/27.21:
        break
print("E = ", E[i])

#Plot
plt.figure()
plt.plot(E[:i])
plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Radial PDF", fontsize=18)
plt.plot(r,4*np.pi*r**2*ns)
plt.savefig('task5.pdf')
# %% Task 6
Z= 2
N=500           #Number of grid points
Rmax = 7        # atomic
h = Rmax/N      #Stepsize
r = np.linspace(h,Rmax,N)

ns=1/np.pi*Z**3*np.exp(-2*Z*r)      #Guess initial density

A = 0.0311
B = -0.048
C = 0.002
D = -0.0116
gamma = -0.1423
beta1 = 1.0529
beta2 = 0.3334

rs=(3/(4*np.pi*ns))**(1/3)

#Exchange and correlation potential
epsx=- 3/4*(3*ns/np.pi)**(1/3)
ndepsx= 1/3*epsx
ndepsc=(rs>=1)*ns*(-(gamma*(beta2+beta1/(2*np.sqrt(rs)))*4*np.pi*rs**4/9)/(1+beta1*np.sqrt(rs)+beta2*rs)**2)
ndepsc=ndepsc+(rs<1)*ns*(A/rs+C+C*np.log(rs)+D)*4*np.pi*rs**4/9

epsc=(rs>=1)*gamma/(1+beta1*np.sqrt(rs)+beta2*rs)
epsc=epsc+(rs<1)*(A*np.log(rs)+B+C*rs*np.log(rs)+D*rs)
epsxc = epsx + epsc

Vx = epsx + ndepsx
Vc = epsc + ndepsc
Vxc = Vx + Vc 


MaxIters = 30           #Max number of iterations
E = np.zeros(MaxIters)

#Solve the self consistency problem
for i in range(MaxIters):
    VH = 2*getVH(ns,N,Rmax)
    (psi,eps) = getpsi(VH+Vxc,N,Rmax,Z)
    ns = np.abs(psi)**2
    E[i] =getE(eps,VH,Vxc,epsxc,ns,Z)
    if np.abs(E[i]-E[i-1])<1e-5/27.21:
        break
print("E = ", E[i])

#Plot
plt.plot(r,4*np.pi*r**2*ns)
plt.figure()
plt.plot(E[:i])


# %%
