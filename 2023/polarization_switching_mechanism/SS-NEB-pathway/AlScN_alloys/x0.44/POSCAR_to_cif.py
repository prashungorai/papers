from glob import iglob
from pymatgen.core import Lattice, Structure


for name in iglob('POSCAR-*'):
    strc = Structure.from_file(name)
    strc.to(filename=name+".cif")
