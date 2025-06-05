# UBEM-Predicted Thermodynamic Stability of Hypothetical Zintl Structures

This folder contains datasets generated with the UBEM graph neural network.

## Contents
 
- `ubem_prediction.csv`
Predicted volume-relaxed energies and decomposition energies for 72,696 hypothetical structures without Yb or Eu, using chemical potentials that are listed in `fere_chemical_potentials.csv`

- `ubem_prediction_dft_validation.csv`
DFT fully-relaxed total energies and decomposition energies for 1,809 structures predicted stable by UBEM. The DFT calculations follow the [NREL Materials Database](https://materials.nrel.gov) workflow. Duplicate ICSD entries have been removed to ensure uniqueness.

- `ubem_prediction_dft_validation_set.json.gz`
Compressed JSON file with the 1,809 DFT fully-relaxed structures referenced above. 

