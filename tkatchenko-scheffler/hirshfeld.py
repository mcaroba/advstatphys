from ase.io import read
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from gpaw import GPAW

# Read in atomic configuration
system = "Benzene-methane_complex"
# If taken from the Github repo, the XYZ file should be availabel locally
# but make it fail safe
try:
    atoms = read("s22/" + system + ".xyz")
except:
    from ase.data.s22 import data as s22
    from ase import Atoms
    atoms = Atoms(s22[system]["symbols"], s22[system]["positions"])
    atoms.center(vacuum = 4.)

# Do a GPAW calculation for this molecule
calc = GPAW(xc="PBE")
atoms.set_calculator(calc)
e = atoms.get_potential_energy()

# Do the Hirshfeld partitioning of the charge density
hp = HirshfeldPartitioning(calc).get_effective_volume_ratios()

# Print the ratios for the atoms in the system
for atom in range(0, len(atoms)):
    print(atoms.symbols[atom], hp[atom])
