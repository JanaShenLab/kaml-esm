# KaML-ESM Inference Pipeline

An end-to-end pKa prediction pipeline using state-of-the-art ESM embeddings
and ensemble MLP inference. Supports single‑sequence, PDB, UniProt and
multi‑FASTA inputs, with configurable channels.

## Overview

KaML‑ESM runs:
- Embedding extraction (ESM2 or ESMC)
- Ensemble inference (AcidicMLP & BasicMLP)
- Structure folding (ESM3-medium)
- Confomer-based inference (KaML-CBTree)
- Results export (CSV, beta-factor-labeled PDB, logs)

## Features

- **Inputs**:  
    • `--seq` (raw sequence)  
    • `--pdb` (local PDB file)  
    • `--pdbid` (fetch PDB by ID)  
    • `--uniprot` (fetch UniProt sequence)  
    • `--fasta` (multi‑FASTA batch)

- **Channels**:  
    • Acidic: `--acidic` (`esm2` or `esmC`)  
    • Basic:  `--basic`  (`esm2` or `esmC`)

- **Structure**:  
    • Default: Forge folding* via `ESM_FORGE_TOKEN`  
    • `--nofold` (do not fold a structure, will disable CBTree unless pdbid/pdb provided)  
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\* in our tests folding takes approximately 20s on average

- **Safety**:  
    • `--skip_safety` to bypass ESM safety filter (permission required)  

- **CBTREE**:  
    • Disabled with `--nocbtree` and `--nofold` unless a pdb or pdbID is input 

## Requirements

- Python 3.10+
- Internet (for Forge API, UniProt, PDB fetch)
- Python packages described in: env/KaML-ESM_env.txt

## Environment Setup

Use the script in `env/` for venv:

    cd env
    ./setup_envs.sh venv

This creates a Python 3.10.12 environment, installs dependencies, and
adds `bin/` to your `PATH`.

## ESM Forge Token

KaML‑ESM requires `ESM_FORGE_TOKEN` for Forge and ESMC:

    export ESM_FORGE_TOKEN=$(cat path/to/forge_token.txt)

Keep this token private.

## Usage Examples

Single sequence:

    kaml-esm --seq "MEEPQSDPSV..." --outdir results/seq1

Fetch by UniProt:

    kaml-esm --uniprot P04637 --outdir results/p53

Fetch PDB:

    kaml-esm --pdbid 1CRN --outdir results/1crn

Multi‑FASTA:

    kaml-esm --fasta proteins.fasta --nproc 4 --outdir results/all

Skip safety filter (requires permission):

    kaml-esm --seq "MEEPQSDPSV..." --skip_safety

Skip structure folding:

    kaml-esm --seq "MEEPQSDPSV..." --nofold

Disable CBTREE:

    kaml-esm --seq "MEEPQSDPSV..." --nocbtree

## Outputs

Default `--outdir` is `output/`, containing:

- `predictions.csv`  — per‑residue pKa, shift, error, optional CBTREE  
- `predicted_structure.pdb`  — updated B‑factors in PDB  
- `pipeline.log`  — debug/info log  
- subfolders for multi‑FASTA runs  

## Code Layout

    bin/kaml-esm            # main CLI script
    bin/kaml-cbtree         # KaML-CBtree helper script
    bin/rida                # rida binary (req. by CBtree)
    bin/mkdssp              # dssp binary (req. by CBtree)
    env/setup_envs.sh       # env setup for python virtual enviroment (venv)
    src/plmpg/esm2          # vendored ESM2 code
    wts/                    # pretrained weights
    README.md               # this file

## Citation

If you use KaML‑ESM, please cite:

> Protein Electrostatic Properties are Fine‑Tuned Through Evolution  
> Shen M., Dayhoff II G.W., Shen J. 2025 (In‑Review)

## License

MIT License

