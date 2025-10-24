# Antibody Sequence Prediction Model

This project is built on ESM (Evolutionary Scale Modeling) and PyTorch for training and inference of antibody repair models.

## Environment Requirements

Before running this project, please ensure the following dependencies are installed:

- Python 3.8
- PyTorch (>=1.9.0)
- ESM (`fair-esm`)

## Data Preparation
The data used in this project is sourced from the OAS database, pre-processed and randomly split into training and test sets:

- Training data: `data_OAS/train.pkl`

- Test data: `data_OAS/test.pkl`

## Quick Start
### 1. Model Training
Run the training script to train the model:

```bash
python train.py
```
After training completes, the model weights file `model.ckpt` will be saved in the corresponding subdirectory under `output/`.

### 2. Model Inference
**Step 1: Prepare Prediction Data**
Create a pickle file containing the antibody sequences to be predicted. The file should be formatted as a list of dictionaries, with each dictionary containing the following keys:
```python
[
    {
        'heavy_chain_id': 'your_heavy_chain_id',
        'light_chain_id': 'your_light_chain_id', 
        'heavy_chain_seq': 'heavy_chain_amino_acid_sequence',
        'light_chain_seq': 'light_chain_amino_acid_sequence'
    },
    # More sequences...
]
```
Important Notes:

- Ensure all heavy and light chains you wish to predict are included in this list

- Initial pairing does not affect prediction

- Sequences should use standard amino acid single-letter codes

**Step 2: Run Inference**
Use the following command for prediction:

```bash
python inference.py --test_data_path path/to/your_test_data.pkl
```
Replace `path/to/your_test_data.pkl` with the actual path to your data file.
