# UBEM: Preprocessing and Training and Pretrained Model

This folder contains scripts and datasets for preprocessing and training UBEM graph neural network (GNN). 

## Preprocess structures
Preprocessing is needed to train the GNN with a different dataset. The input to the `preprocess.py` script is a .csv file of volume-relaxed energies with at least the two columns `id` and `energyperatom` (e.g, `inputs/vol_ref_energies.csv`), and the pymatgen structures in a `.json.gz` (e.g, `inputs/dft_vol_relaxed_structures.json.gz`) file where the key is the structure `id` and the value is the pymatgen structure. Note the `id` in the two files must match.

As an example for how to run `preprocess.py`, here is the command used to preprocess the DFT volume-relaxed dataset:

```
python preprocess.py \
    --dataset inputs/dft_vol_relaxed_structures.json.gz inputs/vol_rel_energies.csv zintl
```

This will create a pandas dataframe as a pickle file with an `inputs` column that contains the preprocessed structures.

## Train GNN
To train the GNN, use the following command:

```
python train_model.py --inputs inputs/preprocessed/scaled_inputs_zintl.p
```

A folder `outputs/` will be created containing the following files:
   - log.csv: contains information on epochs, training and validation losses
   - `best_model.hdf5`: contains model weights
   - `predicted_energies.csv.gz`: contains predicted energy of the test crystal structures. Predicted energy for the test structures are provided in `ubem_training_testing.csv`.

## Pretrained model
- See [250527_zintl_run.ipynb](https://github.com/Rinkumoni119/Zintl_search_UBEM/blob/main/notebooks/250527_zintl_run.ipynb) for how to load the pretrained model (`best_model.hdf5`) and make energy predictions

- The Python utilities and dataset of competing phases needed to calculate decomposition energies are available in the [upper-bound-energy-gnn repository](https://github.com/jlaw9/upper-bound-energy-gnn/)
