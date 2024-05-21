from pymatgen.core.structure import Structure
from pymatgen.analysis.ferroelectricity.polarization import Polarization, calc_ionic, zval_dict_from_potcar, get_total_ionic_dipole
from pymatgen.io.vasp.inputs import Poscar, Potcar
from pymatgen.io.vasp.outputs import Outcar
import numpy as np

### Input preparation 
# for 36 intermediate structures
num_structs = 38   # 36 + two end points
poscars = [Structure.from_file("POSCAR-%02g"%(i)) for i in range(num_structs)] # POSCAR named as POSCAR-00, POSCAR-01, ...., POSCAR-37
outcars = [Outcar("OUTCAR-%02g"%(i)) for i in range(num_structs)]  # OUTCAR names as OUTCAR-00, ..., OUTCAR-37
potcar=Potcar.from_file("POTCAR")
for i in range(num_structs):
    outcars[i].zval_dict = zval_dict_from_potcar(potcar)

### Using VASP outputs to calute (need POSCARs, OUTCARs, and the associated POTCAR)
# calc_ionic_from_zval = True: Use point charge to calculate ionic contribution
pol_from_out_struct_method = Polarization.from_outcars_and_structures(outcars,poscars,
    calc_ionic_from_zval = True)

# all_in_polar=True: use the polar structure to calculate the lattice quanta 
np.savetxt('polarization.dat', pol_from_out_struct_method.get_same_branch_polarization_data(convert_to_muC_per_cm2=True,all_in_polar=True)[:,2])
np.savetxt('quanta.dat', pol_from_out_struct_method.get_lattice_quanta(convert_to_muC_per_cm2=True,all_in_polar=True)[:,2])
