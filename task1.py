# %%
import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import scipy
sns.set('talk')

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
Rmax = 5  # atomic
N = 1000
h = Rmax/N
r = np.linspace(h, Rmax, N)

plt.figure(figsize=(8, 6))
plt.xlabel(r"Distance [$a_0$]", fontsize=18)
plt.ylabel("Radial PDF", fontsize=18)
plt.plot(r, 4*np.pi*r**2 *
         np.sum([C[i]*np.exp(-a[i]*r**2) for i in range(4)], axis=0)**2)
plt.savefig('task1.pdf')




# %%
