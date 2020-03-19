import numpy as np
from gpaw import GPAW
from ase import Atoms
from ase.parallel import parprint

# We loop through different structural parameter values
for d in np.arange(0.92,1.03,0.02):
    for theta in np.arange(90.,120.,5.):

#       Convert angle to radians
        theta *= np.pi/180.

#        Positions of the H2O atoms (in that order)
        pos = [[d*np.cos(theta/2.), d*np.sin(theta/2.), 0.],
               [d*np.cos(theta/2.), -d*np.sin(theta/2.), 0.],
               [0., 0., 0.]]

#       Let us instantiate an Atoms object and add some vacuum around the molecule
        atoms = Atoms("H2O", positions=pos, pbc=False)
        atoms.center(vacuum=4.)

        calc = GPAW(xc="PBE", txt="gpaw.out")
        atoms.set_calculator(calc)

#       Let us run a GPAW calculation to get the cohesive energy of this configuration
        e_water = atoms.get_potential_energy()

        parprint("%f %f %f" % (d, theta*180./np.pi, e_water))
