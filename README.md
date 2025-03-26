# KaML-ESM Inference Pipeline

An end-to-end pKa prediction pipeline using state-of-the-art ESM embeddings and multi-layer perceptron (MLP) inference. This tool accepts multiple input types and offers configurable model choices for acidic and basic channels.

## Overview

The KaML-ESM Inference Pipeline performs the following tasks:
- **Embedding Extraction:**  
  Extracts protein embeddings using either ESM2 or ESMC models. 
- **Inference:**  
  Uses ensemble MLP models (AcidicMLP and BasicMLP) for pKa prediction with parallel processing.
- **Structure Processing:**  
  Optionally folds and optimizes protein structures.
- **Model Selection:**  
  Use the optional arguments `--acidic` and `--basic` to select the model for each channel. By default, acidic uses esm2 and basic uses esmC. If both channels use the same model, embeddings are extracted only once.

## Environment Setup

The repository includes an environment setup script located in the `env` directory. This script supports both Conda and Python virtual environments.

To set up a **Conda** environment, run:

    cd env
    ./setup_envs.sh conda

Or, to set up a **venv** environment, run:

    cd env
    ./setup_envs.sh venv

These scripts do the following:
- Create the environment (e.g., "KaML-ESM") with Python 3.10.12.
- Install all dependencies listed in the respective `*_env.txt` files.
- Modify the activation script to add the repository's `/bin` directory to your PATH so that you can call KaML-ESM from anywhere as long as the environment is active.
- Set up environment variables and alias masking for the different models.

Make sure you have a file named `forge_token.txt` (containing your ESM Forge token) in the repository root, or adjust the path as needed.

## ESM Token Setup

Both esmC and esm3 models require a valid token from the ESM platform to operate. You must obtain a token through ESM and assign it to an environment variable. For example, if your token is saved in a file at `~/wrk/tmp/plm_playground/env/forge_token.txt`, run:

    export PLMPG_FORGE_TOKEN=$(cat ~/wrk/tmp/plm_playground/env/forge_token.txt)

Ensure that this environment variable is set in your shell session before running the pipeline. Keep your token secure and do not share it publicly.

## Usage

After activating the appropriate environment (via Conda or venv), run the pipeline script as follows:

    p53=MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
    KaML-ESM --seq ${p53}

Other input options include:
- `--pdb` for specifying a PDB file,
- `--pdbid` for fetching a PDB structure by ID, or
- `--uniprot` for fetching a sequence from UniProt.

To change the default model selection for each channel, use the `--acidic` and `--basic` options:

    KaML-ESM --seq "YOUR_SEQUENCE" --acidic esm2 --basic esmC

By default, acidic uses esm2 (with the AcidicMLP architecture) and basic uses esmC (with the BasicMLP architecture). If both channels use the same model, embeddings are extracted only once.

The pipeline outputs:
- A CSV file (`predictions.csv`) with predicted pKa values.
- An updated PDB file (if structure output is enabled).
- A detailed log file (`pipeline.log`) in the specified output directory.

## Code Structure

- **bin/KaML-ESM**  
  Main executable for embedding extraction, inference, and structure updates.
  
- **env/setup_envs.sh**  
  Shell script to set up Conda or venv environments and install dependencies.

- **src/plmpg/esm2**  
Vendored ESM2 (i.e. fair-esm) with a wrapper to enable side-by-side esm2 and esm3/C.

- **wts**  
 saved model weights for esm2 and esmC acidic and basic channels

- **README.md**  
  This documentation file.

## Citation

If you use the KaML-ESM Inference Pipeline in your research, please cite:

Protein Electrostatics are Fine-Tuned by Evolution, Shen, M. Dayhoff II G.W, Shen, J. 2025

## License

This project is licensed under the MIT License.

