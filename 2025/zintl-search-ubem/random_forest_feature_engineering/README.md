# Random Forest Model to Interpret Zintl Phase Stability

This contains training datasets and script for calculating the chemical and structural input features. 

## Contents
 
- `rf_prediction_dft_validation_set.json.gz`
Compressed JSON file containing 3,087 DFT fully-relaxed structures of Zintl phases that are predicted to be stable by UBEM and M3GNet approaches, and ICSD phases.
- `rf_prediction_dft_validation.csv`
DFT fully-relaxed total energies and formation energies, using chemical potentials listed in `ubem_prediction_validation/fere_chemical_potentials.csv`. This file also includes RF-predicted formation energies for the same dataset.
- `feature.py`
Functions for computing all the features used in the RF model. Each feature is described in detail in the Supporting Information of the paper. The feature functions utilize the [pylada-light](https://github.com/pylada/pylada-light) package. 

An example demonstration:

```python
from pylada.crystal import read           # Pylada reader
import feature                            # this repo’s feature module

# Load a POSCAR/CONTCAR file
structure = read.poscar("CONTCAR_a1b1c1_33_As1Cd1Na1_icsd_009571_Na_As_Hg")
atom_type_counts = feature.get_atom_types_counts(structure)
print(atom_type_counts)

