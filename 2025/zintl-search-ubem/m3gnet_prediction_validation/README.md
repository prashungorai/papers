# M3GNet-Predicted Thermodynamic Stability of Hypothetical Zintl Structures

This folder contains datasets generated with the [M3GNet](https://github.com/materialsvirtuallab/m3gnet) graph neutral network.

##Contents
 
- `m3gnet_relaxed_structures.tar.gz`
Compressed file containing 90,805 crystal structures relaxed using the M3GNet model.
- `m3gnet_prediction.csv`
Predicted total energies and decomposition energies of the 90,805 M3GNet-relaxed structures. Decomposition energies are calculated based on competing phase from the Materials Project.
- `m3gnet_prediction_dft_validation.csv`
DFT fully-relaxed total energies and decomposition energies of 1,996 structures predicted stable by M3GNet. The DFT calculations follow the DFT setup of Materials Project. Structures containing Eu and Yb are included in this dataset. Duplicates corresponding to ICSD entries have been removed to ensure uniqueness.

