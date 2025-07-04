from pymatgen.core.structure import Structure
from pymatgen.analysis.ferroelectricity.polarization import Polarization, calc_ionic, zval_dict_from_potcar, get_total_ionic_dipole
from pymatgen.io.vasp.inputs import Poscar, Potcar
from pymatgen.io.vasp.outputs import Outcar
import numpy as np
from glob import glob


num_structs = 41 #39 + 2 end points
poscars = []
outcars = []
for i in sorted(glob('*/')):
    poscars.append(Structure.from_file(i+"POSCAR"))
    outcars.append(Outcar(i+"OUTCAR"))

potcar=Potcar.from_file("POTCAR")

for i in range(num_structs):
    outcars[i].zval_dict = zval_dict_from_potcar(potcar)


pol_from_out_struct_method = Polarization.from_outcars_and_structures(outcars,poscars,
    calc_ionic_from_zval = True)
print(pol_from_out_struct_method.get_same_branch_polarization_data(convert_to_muC_per_cm2=True,all_in_polar=True))
np.savetxt('polarization_diff_quanta.dat', pol_from_out_struct_method.get_same_branch_polarization_data(convert_to_muC_per_cm2=True,all_in_polar=True)[:,2])
np.savetxt('quanta.dat', pol_from_out_struct_method.get_lattice_quanta(convert_to_muC_per_cm2=True,all_in_polar=True)[:,2])
