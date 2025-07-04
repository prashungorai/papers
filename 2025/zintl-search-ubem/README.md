# Charting the Chemical Space of Zintl Phases with Graph Neural Networks and Bonding Insights

A GNN for predicting thermodynamic stability of hypothetical Zintl structures using the upper bound energy minimization approach (UBEM) and Random Forest (RF) model for chemical interpretation

This package provides the following functionalities:

- Predict volume-relaxed energy of hypothetical Zintl structures using pre-trained GNN model (trained on Zintl structures)
- Train GNN model with custom datasets to predict volume-relaxed energy
- Python script to calculate chemical and structural features of Zintl phases
  
## Installation

Most dependencies for this project are installable via mamba/conda, with the exception of [nfp](https://github.com/NREL/nfp), which can be installed via pip.

```yaml
channels:
  - conda-forge
  - defaults

dependencies:
  - python=3.7
  - jupyterlab
  - seaborn
  - pandas
  - scikit-learn
  - jupyter
  - notebook
  - pymatgen
  - tqdm
  - tensorflow-gpu
  - pip
    - pip:
    - nfp >= 0.3.12
    - tensorflow-addons
```

## Contents

Each structure is associated with an `id`, which begins with `POSCAR` if it is a hypothetical structure, or with `CONTCAR` if it is a relaxed structure (via DFT or machine learning). The energy associated with each `id` — whether volume-relaxed energy, total energy, decomposition energy, or formation energy — is expressed in eV/atom.

- Hypothetical structures: `decorations_prototypes.json.gz`
  
- Preprocess structures and train a new scale-invariant GNN: `ubem_training_pretrained_model`

- Training structures and energies: 
  - UBEM: `ubem_training_pretrained_model`
  - RF: `random_forest_feature_engineering`

- Prediction and validation:
  - UBEM: `ubem_prediction_validation`
  - M3GNet: `m3gnet_prediction_validation`
  - RF: `random_forest_feature_engineering`

- Calculate features for RF model: `random_forest_feature_engineering`



The files in this project are curated and maintained by

* [Rinkumoni Chaliha](mailto:chalir[at]rpi[dot]edu)
* [Manish Kumar Kothakonda](mailto:mkothako[at]stanford[dot]edu)
* [Cheng-Wei Lee](mailto:clee2[at]mines[dot]edu)
* [Jeffrey N. Law](mailto:jeffrey[dot]law[at]nrel[dot]gov)
* [Prashun Gorai](mailto:goraip[at]rpi[dot]edu)

### Cite
```
@article{chaliha2025,
  title={Charting the Chemical Space of Zintl Phases with Graph Neural Networks and Bonding Insights},
  author={Chaliha, Rinkumoni and Kothakonda, Manish K and Lee, Cheng-Wei and Law, Jeffrey N and Yang, Qian and Bobev, Svilen and Gorai, Prashun},
  journal={ChemRxiv},
  year={2025},
  doi = {dx.doi.org/10.xxx/xxx}
}
```

### License

This package is released under the MIT License.
