from ase import Atoms
from gpaw import GPAW
import numpy as np

# Equilibrium for PBE:
d = 0.974
a = 104.4

calc = GPAW(xc="PBE", txt = "out")

positions = [[d*np.cos(a/2./180.*np.pi), d*np.sin(a/2./180.*np.pi), 0.],
             [d*np.cos(a/2./180.*np.pi), - d*np.sin(a/2./180.*np.pi), 0.],
             [0., 0., 0.]
            ]
atoms = Atoms("H2O", positions = positions)
atoms.center(vacuum = 4.)
atoms.set_calculator(calc)

e = atoms.get_potential_energy()

print("%f %f %f" % (d, a, e))


