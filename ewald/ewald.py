import numpy as np
from scipy.special import erfc as erfc
import sys

#***************************************
# USER defined parameters
#
# Select mu in A^-1
#mu = 2.
mu = float(sys.argv[1])

# Body centered cubic CsCl example
a = 4.119
cell = [[a, 0., 0.],
        [0., a, 0.],
        [0., 0., a]]

asym = [[0., 0., 0.],
        [a/2., a/2., a/2.]]

charges = [+1., -1.]
#***************************************

# DO NOT CHANGE THE CODE AFTER THIS UNLESS YOU
# KNOW WHAT YOU'RE DOING

# Fundamental constants
e = 1.60217662e-19
eps0 = 8.8541878128e-12

# Electrostatic prefactor in eV.A (note we removed one e)
pref = e / 8. / np.pi / eps0 * 1.e10

# Process the lattice information
V_cell = np.dot(cell[0], np.cross(cell[1],cell[2]))
cell_rec = np.zeros([3,3])
cell_rec[0] = 2.*np.pi * np.cross(cell[1], cell[2]) \
              / np.dot(cell[0], np.cross(cell[1],cell[2]))
cell_rec[1] = 2.*np.pi * np.cross(cell[2], cell[0]) \
              / np.dot(cell[1], np.cross(cell[2],cell[0]))
cell_rec[2] = 2.*np.pi * np.cross(cell[0], cell[1]) \
              / np.dot(cell[2], np.cross(cell[0],cell[1]))


# Computes the total charge and the energy of the compensating
# background
Q = np.array(charges).sum()
E_HNB = -pref * (np.pi * Q**2) / (V_cell * mu**2)

# Computes the self-energy correction to the reciprocal sum
E_self = pref * (-2.*mu/np.sqrt(np.pi)) * (np.array(charges)**2).sum()

#print(E_HNB, E_self)

# Generate a bunch of reciprocal lattice translations
n_rec = [10,10,10]
Gs = np.zeros([(2*n_rec[0]+1)*(2*n_rec[1]+1)*(2*n_rec[2]+1), 3])
n = 0
# This is very inefficient; I'm sure there's a numpy way to do it faster
for i in range(-n_rec[0], n_rec[0]+1):
    for j in range(-n_rec[1], n_rec[1]+1):
        for k in range(-n_rec[2], n_rec[2]+1):
            Gs[n] = float(i)*cell_rec[0] + \
                    float(j)*cell_rec[1] + \
                    float(k)*cell_rec[2]
            n += 1

# Sort them according to lowest norm
Gs_norm = np.sqrt( Gs[:,0]**2 + Gs[:,1]**2 + Gs[:,2]**2 )
sort_indices = np.argsort(Gs_norm)

# Do the reciprocal space sum
f = open("rec.dat", "w+")
E_rec = 0.
n = 0
for i in sort_indices[1:]:
    G_norm = Gs_norm[i]
    G = Gs[i]
    charge_rec = 0.
    for j in range(0, len(asym)):
        charge_rec += charges[j] * np.exp( -1j* np.dot(G, asym[j]) )
    charge_rec = np.abs(charge_rec)**2
    E_rec += pref * 4.*np.pi/V_cell * np.exp( -G_norm**2/4./mu**2 ) \
            / G_norm**2 * charge_rec
    if n % 1 == 0:
        print("%i %f" % (n, E_rec+E_self), file=f)
    n += 1

f.close()


# Generate a bunch of real lattice translations
n_real = [10,10,10]
Ts = np.zeros([(2*n_real[0]+1)*(2*n_real[1]+1)*(2*n_real[2]+1), 3])
n = 0
# This is very inefficient; I'm sure there's a numpy way to do it faster
cell = np.array(cell)
for i in range(-n_real[0], n_real[0]+1):
    for j in range(-n_real[1], n_real[1]+1):
        for k in range(-n_real[2], n_real[2]+1):
            Ts[n] = float(i)*cell[0] + float(j)*cell[1] + float(k)*cell[2]
            n += 1

# Sort them according to lowest norm
Ts_norm = np.sqrt( Ts[:,0]**2 + Ts[:,1]**2 + Ts[:,2]**2 )
sort_indices = np.argsort(Ts_norm)

# Do the reciprocal space sum
f = open("real.dat", "w+")
E_real = 0.
# Do the T = 0 sum first, skipping the self term
for i in range(0, len(asym)):
    for j in range(0, len(asym)):
        d = np.dot(np.array(asym[j]) - np.array(asym[i]), \
                   np.array(asym[j]) - np.array(asym[i]))**0.5
        if i != j:
            E_real += pref* charges[i]*charges[j] * erfc(mu*d)/d

print("%i %f" % (0, E_real+E_HNB), file=f)

n = 1
# Now all other Ts including the "self term"
for i in sort_indices[1:]:
    T = Ts[i]
#    print(T)
    for j in range(0, len(asym)):
        for k in range(0, len(asym)):
            d = np.dot(np.array(asym[k]) - np.array(asym[j]) + T, \
                       np.array(asym[k]) - np.array(asym[j]) + T)**0.5
            E_real += pref* charges[j]*charges[k] * erfc(mu*d)/d
    if n % 1 == 0:
        print("%i %f" % (n, E_real+E_HNB), file=f)
    n += 1

f.close()
