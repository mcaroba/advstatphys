from ase.io import read
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from gpaw.analyse.vdwradii import vdWradii
from gpaw import GPAW
from ase.calculators.vdwcorrection import vdWTkatchenko09prl
import numpy as np

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

# Now let's shift the position of the methane molecule relative
# to the benzene
for z_shift in np.arange(-1.,2.,0.1):
    atoms_shifted = atoms.copy()
    for i in range(12,17):
        atoms_shifted.positions[i] += [0, 0, z_shift]
    atoms_shifted.center(vacuum = 4.)
    calc = GPAW(xc="PBE", txt="gpaw.out")
    atoms_shifted.set_calculator(calc)
    e = atoms_shifted.get_potential_energy()
    calc_ts = vdWTkatchenko09prl(HirshfeldPartitioning(calc),
                 vdWradii(atoms_shifted.get_chemical_symbols(), 'PBE'))
    atoms_shifted.set_calculator(calc_ts)
    e_ts = atoms_shifted.get_potential_energy()
    print("%f %f %f" % (z_shift, e, e_ts))
