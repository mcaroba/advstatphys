import numpy as np
import sys

# Hyperparameters
#sigma = 0.001
#sigma2b = 0.5
sigma2b = float(sys.argv[1])
sigma = float(sys.argv[2])
sigma_jitter = 1.e-8
cutoff = 5.

# Number of training data, etc
Ntrain = 150
Nsparse = 15

# Sparse configurations
dOOs = np.arange(1.5, 5.0001, 3.5/float(Nsparse-1))
dOHs = np.arange(1.5, 5.0001, 3.5/float(Nsparse-1))
dHHs = np.arange(1.5, 5.0001, 3.5/float(Nsparse-1))

# Vacuum permittivity in SI units
eps0 = 8.8541878e-12
e0 = 1.60217662e-19

# Electrostatic pontetial in eV
def elec(Z1, Z2, d12):
    return 1.e10 / 4. / np.pi / eps0 * e0 * Z1 * Z2 / d12

# 2-body kernel:
def k2b(d1, d2, sigma2b, cutoff):
    if d1 > cutoff or d2 > cutoff:
        return 0.
    else:
        return np.exp(-0.5*(d1-d2)**2/sigma2b**2)

# Partial charges
ZH = 0.32
ZO = -2.*ZH

# Reference energy
E0 = -14.620

# Read in training data from file
f = open("random_dimer.dat", "r")
raw_energies = []
descriptors = []
for line in f:
    raw_energies.append(float(line.split()[0]))
    descriptors.append([float(a) for a in line.split()[1:]])

f.close()

# Transform from raw energy to just the repulsion part
energies = []
for i in range(0, len(raw_energies)):
    raw_energy = raw_energies[i]
    dOO, dOH1, dOH2, dOH3, dOH4, dHH1, dHH2, dHH3, dHH4 = descriptors[i]
    energy = raw_energy - 2.*E0
    energy -= elec(ZO,ZO,dOO)
    energy -= elec(ZO,ZH,dOH1) + elec(ZO,ZH,dOH2)
    energy -= elec(ZO,ZH,dOH3) + elec(ZO,ZH,dOH4)
    energy -= elec(ZH,ZH,dHH1) + elec(ZH,ZH,dHH2)
    energy -= elec(ZH,ZH,dHH3) + elec(ZH,ZH,dHH4)
    energies.append(energy)


rmse = 0.
f = open("cross_validation.dat", "w+")
for ex_config in range(0,Ntrain):
#   Build the covariance matrix (this could be done only once but whatever)
    CTS = np.zeros([Ntrain, 3*Nsparse])
    for t in range(0, Ntrain):
        dOO, dOH1, dOH2, dOH3, dOH4, dHH1, dHH2, dHH3, dHH4 = descriptors[t]
#       OO GAP
        for s in range(0, Nsparse):
            ds = dOOs[s]
            CTS[t,s] = k2b(ds, dOO, sigma2b, cutoff)
#       OH GAP
        for s in range(0, Nsparse):
            s0 = Nsparse
            ds = dOHs[s]
            CTS[t,s+s0] = k2b(ds, dOH1, sigma2b, cutoff)
            CTS[t,s+s0] += k2b(ds, dOH2, sigma2b, cutoff)
            CTS[t,s+s0] += k2b(ds, dOH3, sigma2b, cutoff)
            CTS[t,s+s0] += k2b(ds, dOH4, sigma2b, cutoff)
#       HH GAP
        for s in range(0, Nsparse):
            s0 = 2*Nsparse
            ds = dHHs[s]
            CTS[t,s+s0] = k2b(ds, dHH1, sigma2b, cutoff)
            CTS[t,s+s0] += k2b(ds, dHH2, sigma2b, cutoff)
            CTS[t,s+s0] += k2b(ds, dHH3, sigma2b, cutoff)
            CTS[t,s+s0] += k2b(ds, dHH4, sigma2b, cutoff)

    CTS_train = np.delete(CTS, (ex_config), axis=0)

#   This is the regularizing matrix
    LTT_inv = np.zeros([Ntrain-1,Ntrain-1])
    for t in range(0, Ntrain-1):
        LTT_inv[t,t] = 1./sigma

#   This is CSS
    CSS = np.zeros([3*Nsparse, 3*Nsparse])
    for s1 in range(0, Nsparse):
        CSS[s1,s1] = 1. + sigma_jitter
        ds1 = dOOs[s1]
        for s2 in range(s1+1, Nsparse):
            ds2 = dOOs[s2]
            CSS[s1,s2] = k2b(ds1, ds2, sigma2b, cutoff)
            CSS[s2,s1] = CSS[s1,s2]

    for s1 in range(0, Nsparse):
        CSS[s1+Nsparse,s1+Nsparse] = 1. + sigma_jitter
        ds1 = dOHs[s1]
        for s2 in range(s1+1, Nsparse):
            ds2 = dOHs[s2]
            CSS[s1+Nsparse,s2+Nsparse] = k2b(ds1, ds2, sigma2b, cutoff)
            CSS[s2+Nsparse,s1+Nsparse] = CSS[s1+Nsparse,s2+Nsparse]

    for s1 in range(0, Nsparse):
        CSS[s1+2*Nsparse,s1+2*Nsparse] = 1. + sigma_jitter
        ds1 = dHHs[s1]
        for s2 in range(s1+1, Nsparse):
            ds2 = dHHs[s2]
            CSS[s1+2*Nsparse,s2+2*Nsparse] = k2b(ds1, ds2, sigma2b, cutoff)
            CSS[s2+2*Nsparse,s1+2*Nsparse] = CSS[s1+2*Nsparse,s2+2*Nsparse]

    this_energies = np.delete(energies, ex_config)
#   These are the alphas
    C_temp = np.dot(LTT_inv, CTS_train)
    C_temp = np.dot(np.transpose(CTS_train), C_temp)
    C_temp = CSS + C_temp
    C_temp = np.linalg.inv(C_temp)
    C_temp = np.dot(C_temp, np.transpose(CTS_train))
    C_temp = np.dot(C_temp, LTT_inv)
    alphas = np.dot(C_temp, this_energies)

#   Now we make the predictions
    pred = np.dot(CTS[ex_config,:], alphas)
    rmse += (energies[ex_config]-pred)**2
    print("%f %f" % (energies[ex_config], pred), file=f)

f.close()

rmse = np.sqrt( rmse/float(Ntrain) )

# Print hyperparameters and error
print("%f %f %f" % (sigma2b, sigma, rmse))
