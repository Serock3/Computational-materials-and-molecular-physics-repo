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

def solve(ns,MaxIters,N,Rmax,Z):
    #Solve the self consistency problem
    for i in range(MaxIters):
        VH = getVH(ns,N,Rmax)
        (psi,eps) = getpsi(VH,N,Rmax,Z)
        ns = np.abs(psi)**2
        E[i] = getE(eps,VH,0,0,ns,Z)
        if np.abs(E[i]-E[i-1])<1e-5/27.21:
            break
    return (ns, E)

#%% Task 1
# Given alpha values
a = [0.297104, 1.236745, 5.749982, 38.216677]

# Init of matrices from Thijssen 4.3.2
# h matrix
h = np.zeros((4, 4))
for p in range(4):
    for q in range(4):
        h[p, q] = 4*np.pi/(a[p]+a[q])*(3/4*a[q] *
                                       (1-a[q]/(a[p]+a[q]))*np.sqrt(np.pi/(a[p]+a[q]))-1)

# Q matrix (given by eq 4.17 in Thijssen)
Q = np.zeros((4, 4, 4, 4))
for p in range(4):
    for q in range(4):
        for r in range(4):
            for s in range(4):
                Q[p, r, q, s] = 2*np.pi**(5/2)/((a[p]+a[q])
                                                * (a[r]+a[s])*np.sqrt(a[p]+a[q]+a[r]+a[s]))

# S matrix
S = np.zeros((4, 4))
for p in range(4):
    for q in range(4):
        S[p, q] = (np.pi/(a[p]+a[q]))**(3/2)

# Inital values
C = [1, 1, 1, 1]

#Normalizing C according to eq. 4.19 in Thijssen
def normalize(C):
    return C/np.sqrt(np.matmul(C, np.matmul(S, C)))


F = np.zeros((4, 4))

#Evaluating eq. 4.21 in Thijssen
def getEg(C, h, Q):
    def con(Q, C):
        return np.tensordot(Q, C, axes=([0], [0]))
    return 2*np.matmul(C, np.matmul(h, C))+con(con(con(con(Q, C), C), C), C)


MaxIters = 15
E = np.zeros(MaxIters)

for i in range(MaxIters):
    C = normalize(C)
    E[i] = getEg(C, h, Q)
    print(E[i])

    if np.abs(E[i]-E[i-1]) < 1e-5/27.21:
        break
    # Set F(C)
    for p in range(4):
        for q in range(4):
            F[p, q] = h[p, q]+np.matmul(C, np.matmul(Q[p, :, q, :], C))

    # Solve eigenvalue problem
    (eps, w) = scipy.linalg.eig(F, S)

    if any(eps.imag != 0):
        raise Exception('complex eig')

    # Get best C
    C = w[:, np.argmin(eps.real)]

C = normalize(C)
print("Eg = ", getEg(C, h, Q))
print("C= ", C)

# Print the radial PDF
Rmax = 7  # atomic
N = 1000
h = Rmax/N
r = np.linspace(h, Rmax, N)

plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Radial PDF", fontsize=18)
plt.plot(r, 4*np.pi*r**2 *
         np.sum([C[i]*np.exp(-a[i]*r**2) for i in range(4)], axis=0)**2)
plt.savefig('task1.pdf')

task1density = 4*np.pi*r**2 * np.sum([C[i]*np.exp(-a[i]*r**2) for i in range(4)], axis=0)**2
task1r = r
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

task3density = 4*np.pi*r**2*ns
task3r = r
#%% Task 4 Rmax conv
Z= 2
N=800               #Number of grid points
Rmax = 6            # atomic
h = 0.006          #Stepsize
r = np.linspace(h,Rmax,N)
ns=1/np.pi*Z**3*np.exp(-2*Z*r)      #Guess initial density

MaxIters = 15           #Max number of iterations
E = np.zeros(MaxIters)  

conIters = 8
Rmaxlist = np.linspace(3,10,conIters)
Econv = np.zeros(conIters)

#Find best Rmax
for i in range(conIters):
    N = int(np.round(Rmaxlist[i]/h))
    r = np.linspace(h,Rmaxlist[i],N)
    ns=1/np.pi*Z**3*np.exp(-2*Z*r)      #Guess initial density
    E = solve(ns,MaxIters,N,Rmaxlist[i],Z)[1] # Solves the self consistency problem
    Econv[i] = E[E!=0][-1]
    print("E = ", Econv[i], "at rmax", Rmaxlist[i]," and i ",i)
    if np.abs(Econv[i]-Econv[i-1])<1e-5/27.21:
        print("Within convergence criteria at rmax",Rmaxlist[i]," and i ",i)
    
plt.figure(figsize=(6.15, 4.6))
plt.plot(Rmaxlist[:i],Econv[:i])
plt.xlabel(r"$R_{max}$ [$a_0$]", fontsize=18)
plt.ylabel(r"Energy [eV]", fontsize=18)
plt.savefig('task4_rmax.pdf')

ns = solve(ns,MaxIters,N,Rmaxlist[i],Z)[0]
#Plot ns
plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Radial PDF", fontsize=18)
plt.plot(r,4*np.pi*r**2*ns)
# plt.savefig('task4.pdf')

#Plot VH
plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("VH", fontsize=18)
plt.plot(r,getVH(ns,N,Rmaxlist[i]))
# plt.savefig('task4.pdf')
#%% Task 4 h conv
Rmax = 7            # atomic
h = 0.006          #Stepsize
MaxIters = 15           #Max number of iterations
E = np.zeros(MaxIters)  

conIters = 9
hlist = np.logspace(np.log10(0.03),np.log10(0.003),conIters)
Econv = np.zeros(conIters)

#Find best Rmax
for i in range(conIters):
    N = int(np.round(Rmax/hlist[i]))
    r = np.linspace(hlist[i],Rmax,N)
    ns=1/np.pi*Z**3*np.exp(-2*Z*r)      #Guess initial density
    E = solve(ns,MaxIters,N,Rmax,Z)[1] # Solves the self consistency problem
    Econv[i] = E[E!=0][-1]
    print("E = ", Econv[i], "at h", hlist[i]," and i ",i)
    if np.abs(Econv[i]-Econv[i-1])<1e-5/27.21:
        print("Within convergence criteria at h = ",hlist[i]," and i = ",i)
    
plt.figure(figsize=(8, 6))
plt.plot(hlist[:i],Econv[:i])
plt.xlabel(r"h [$a_0$]", fontsize=18)
plt.ylabel("Energy [eV]", fontsize=18)
plt.xlim(0.03,0)
plt.savefig('task4_h.pdf')

ns = solve(ns,MaxIters,N,Rmax,Z)[0]
#Plot ns
plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Radial PDF", fontsize=18)
plt.plot(r,4*np.pi*r**2*ns)
# plt.savefig('task4.pdf')
task4density = 4*np.pi*r**2*ns
task4r = r

#%% plot together
plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Radial PDF", fontsize=18)
plt.plot(task1r,task1density, label ="Helium Simple Hartree Fock")
plt.plot(task3r,task3density, label ="Hydrogen Kohn Sham")
plt.plot(task4r,task4density, label ="Helium Hartree Fock")
plt.legend()
#%% Task 5
Z= 2
N=1000          #Number of grid points
Rmax = 7        # atomic
h = Rmax/N      #Stepsize
r = np.linspace(h,Rmax,N)

ns=1/np.pi*Z**3*np.exp(-2*Z*r)      #Guess initial density

#Exchange potential

def getepsx(n,Z):
    return -3/4*(3*Z*ns/np.pi)**(1/3)

def getVx(n,Z):
    return -1*(3*Z*ns/np.pi)**(1/3)


epsx= getepsx(Z*ns,Z)
#ndepsx= 1/3*epsx
Vx = getVx(Z*ns,Z)

MaxIters = 30           #Max number of iterations
E = np.zeros(MaxIters)
Test = np.zeros(MaxIters)

#Solve the self consistency problem
for i in range(MaxIters):
    VH = Z*getVH(ns,N,Rmax)
    (psi,eps) = getpsi(VH+Vx,N,Rmax,Z)
    ns = np.abs(psi)**2
    
    epsx= getepsx(Z*ns,Z)
    Vx = getVx(Z*ns,Z)
    
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
N=1000           #Number of grid points
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

def getrs(n,Z):
    return (3/(4*np.pi*Z*ns))**(1/3)

def getepsc(n,Z):
    rs=getrs(n,Z)
    epsc=(rs>=1)*gamma/(1+beta1*np.sqrt(rs)+beta2*rs)
    epsc=epsc+(rs<1)*(A*np.log(rs)+B+C*rs*np.log(rs)+D*rs)
    return epsc

def getepsx(n,Z):
    return -3/4*(3*Z*ns/np.pi)**(1/3)

def getVx(n,Z):
    return -1*(3*Z*ns/np.pi)**(1/3)

def getVc(n,Z):
    rs=getrs(n,Z)
    Vc=(rs>=1)*getepsc(n,Z)*(1+7/6*beta1*np.sqrt(rs)+beta2*rs)/(1+beta1*np.sqrt(rs)+beta2*rs)
    Vc+=(rs<1)*(A*np.log(rs)+B-A/3+2/3*C*rs*np.log(rs)+(2*D-C)*rs/3)
    return Vc

#ndepsc=(rs>=1)*Z*ns*((gamma*(beta2+beta1/(2*np.sqrt(rs)))*4*np.pi*rs**4/9)/(1+beta1*np.sqrt(rs)+beta2*rs)**2)
#ndepsc=ndepsc+(rs<1)*-Z*ns*(A/rs+C+C*np.log(rs)+D)*4*np.pi*rs**4/9


epsxc = getepsx(Z*ns,Z)+getepsc(Z*ns,Z)

Vxc = getVx(Z*ns,Z) + getVc(Z*ns,Z) 


MaxIters = 30           #Max number of iterations
E = np.zeros(MaxIters)

#Solve the self consistency problem
for i in range(MaxIters):
    VH = 2*getVH(ns,N,Rmax)
    (psi,eps) = getpsi(VH+Vxc,N,Rmax,Z)
    ns = np.abs(psi)**2
    
    #Update eps and V with new ns
    epsxc = getepsx(Z*ns,Z) + getepsc(Z*ns,Z)
    Vxc = getVx(Z*ns,Z) + getVc(Z*ns,Z)

    E[i] =getE(eps,VH,Vxc,epsxc,ns,Z)
    if np.abs(E[i]-E[i-1])<1e-5/27.21:
        break
print("E = ", E[i])

#Plot
plt.plot(r,4*np.pi*r**2*ns)
plt.figure()
plt.plot(E[:i])


# %%
