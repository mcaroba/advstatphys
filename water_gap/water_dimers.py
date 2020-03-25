from ase import Atoms
from gpaw import GPAW
import numpy as np
from ase.visualize import view

# True if you want the system's ball and stick visual model to pop up
visualize = True

print_to_file = False

# Equilibrium for PBE:
d = 0.974
a = 104.4

calc = GPAW(xc="PBE", txt = "out")

# Write data to this file. We open in append mode in case we have previously
# computed more configurations
if print_to_file:
    f = open("random_dimer.dat", "a+")

# Select number of configurations
n_config = 150

# We generate random separations (d_OO) and orientations
d12_list = np.random.sample(n_config)
rot_list = np.random.sample([n_config,3])


for config in range(0, n_config):
#   Positions for molecule 1
    pos1 = [[-d*np.cos(a/2./180.*np.pi), d*np.sin(a/2./180.*np.pi), 0.],
            [-d*np.cos(a/2./180.*np.pi), - d*np.sin(a/2./180.*np.pi), 0.],
            [0., 0., 0.],
           ]
#   Positions for molecule 2 before translation and rotation
    pos2 = [[d*np.cos(a/2./180.*np.pi), d*np.sin(a/2./180.*np.pi), 0.],
            [d*np.cos(a/2./180.*np.pi), - d*np.sin(a/2./180.*np.pi), 0.],
            [0., 0., 0.]
           ]
    atoms1 = Atoms("H2O", positions = pos1)
    atoms2 = Atoms("H2O", positions = pos2)

#   Oxygens are separated by at least 1.5 Angstrom and up to 5 Angstrom
    d12 = 1.5 + 3.5*d12_list[config]

#   Rotation angles for the second molecule
    rotx = rot_list[config][0]
    roty = rot_list[config][1]
    rotz = rot_list[config][2]

#   Rotate molecule 2
    atoms2.rotate( rotx * 360., "x")
    atoms2.rotate( roty * 360., "y")
    atoms2.rotate( rotz * 360., "z")

#   Translate molecule 2
    atoms2.translate([d12, 0., 0.])

#   System is made up of both molecules
    atoms = atoms1 + atoms2
    atoms.center(vacuum = 4.)

    if visualize:
        view(atoms)

#   Get potential energy
    atoms.set_calculator(calc)
    e = atoms.get_potential_energy()

#   Get descriptors
    dOO = atoms.get_distance(2,5)
    dOH1 = atoms.get_distance(2,3)
    dOH2 = atoms.get_distance(2,4)
    dOH3 = atoms.get_distance(0,5)
    dOH4 = atoms.get_distance(1,5)
    dHH1 = atoms.get_distance(0,3)
    dHH2 = atoms.get_distance(0,4)
    dHH3 = atoms.get_distance(1,3)
    dHH4 = atoms.get_distance(1,4)

    if print_to_file:
        print("%f %f %f %f %f %f %f %f %f %f" % (e, dOO, dOH1, dOH2, dOH3, dOH4,
                                              dHH1, dHH2, dHH3, dHH4), file=f )
    else:
        print("%f %f %f %f %f %f %f %f %f %f" % (e, dOO, dOH1, dOH2, dOH3, dOH4,
                                              dHH1, dHH2, dHH3, dHH4) )
