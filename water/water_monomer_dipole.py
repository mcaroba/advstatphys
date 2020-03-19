import numpy as np
from gpaw import GPAW
from ase import Atoms
from ase.visualize import view

# We set this to True if we want lots of info on the screen,
# otherwise False
print_to_screen=False

# Initial structural parameters
d = 0.974
theta = 104.4 * np.pi/180.

# Positions of the H2O atoms (in that order)
pos = [[d*np.cos(theta/2.), d*np.sin(theta/2.), 0.],
       [d*np.cos(theta/2.), -d*np.sin(theta/2.), 0.],
       [0., 0., 0.]]

# Let us instantiate an Atoms object and add some vacuum around the molecule
atoms = Atoms("H2O", positions=pos, pbc=False)
atoms.center(vacuum=4.)

# Let's view the molecule
view(atoms)

# We create a GPAW calculator object and assign it to atoms
if print_to_screen:
    calc = GPAW(xc="PBE")
else:
    calc = GPAW(xc="PBE", txt="gpaw.out")

atoms.set_calculator(calc)

# Let us run a GPAW calculation to get the cohesive energy of this configuration
e_water = atoms.get_potential_energy()
dip_water = atoms.get_dipole_moment()

print("Dipole of a H2O monomer: (%f, %f, %f) electrons/Angstrom" % \
      (dip_water[0], dip_water[1], dip_water[2]))
