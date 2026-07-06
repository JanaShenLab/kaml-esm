# KaML: Acid-Base Constant Prediction From Protein Sequence (ESM) or Structure (CBT)

Companion repository for:
-  Shen, Mingzhe, Dayhoff II, Guy W., et al. "Protein Electrostatic Properties are Fine-Tuned Through Evolution." bioRxiv (2025): 2025-09. (https://www.biorxiv.org/content/10.1101/2025.04.17.649309v1)

## How to cite KaML

If you use the KaML models, command-line interface, or precomputed human-wide pKa predictions
in your work, please cite:

    Shen, M., Dayhoff II, G.W., et al., 2025. Protein Electrostatic 
    Properties are Fine-Tuned Through Evolution. bioRxiv (2025)

    Shen, Mingzhe, et al. "KaMLs for Predicting Protein pKa Values and 
    Ionization States: Are Trees All You Need?." JCTC 21.3 (2025): 1446-1458.

BibTeX:

    @article{shen2025kamlesm,
      title={Protein Electrostatic Properties are Fine-Tuned Through Evolution},
      author={Shen, Mingzhe and Dayhoff II, Guy W and Shen, Jana},
      journal={bioRxiv},
      pages={2025--04},
      year={2025},
      publisher={Cold Spring Harbor Laboratory}
    }

    @article{shen2025kamls,
      title={KaMLs for Predicting Protein p K a Values and Ionization States: Are Trees All You Need?},
      author={Shen, Mingzhe and Kortzak, Daniel and Ambrozak, Simon and Bhatnagar, Shubham and Buchanan, Ian and Liu, Ruibin and Shen, Jana},
      journal={Journal of Chemical Theory and Computation},
      volume={21},
      number={3},
      pages={1446--1458},
      year={2025},
      publisher={ACS Publications}
    }

## License ## 
**License** Content © 2025 Mingzhe Shen, Guy W. Dayhoff II and Jana Shen, licensed under
[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

Unless otherwise noted, this repository’s content is © 2025 Guy W. Dayhoff II (on behalf of all authors) and is licensed under Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0).

Plain English: You may copy, adapt, and share the repository’s content for non-commercial purposes as long as you provide proper attribution. For any commercial use, please contact the authors for permission.

# Web-based Inference and Human proteome-wide pKas

We provide a browser-based KaML interface for both on-demand
predictions and interactive exploration of a precomputed human
proteome atlas:

    Home:      https://kaml.computchem.org
    Atlas:     https://kaml.computchem.org/human

The web tools run the same KaML models described in the manuscript and
use ESM-2 (`esm2_t33_650M_UR50D`) and/or ESM-C (`esmc-6b-2024-12`) via 
the EvolutionaryScale Forge API for protein embeddings where on-demand 
inference is required.

------------------------------------------------------------
Web-based inference: https://kaml.computchem.org
------------------------------------------------------------

The web-base interface provides a simple front end to the KaML inference
pipeline for user-supplied sequences.

Requirements for KaML-ESM2:
- A modern web browser.

Requirements for KaML-ESMC:
- An active ESM Forge API token issued by EvolutionaryScale.

The KaML web interface does not issue ESM Forge tokens. Each user must
obtain and manage their own token directly from EvolutionaryScale.

Obtaining an ESM Forge API token:

1. Visit the EvolutionaryScale Forge site:

       https://forge.evolutionaryscale.ai

2. Sign up or sign in with your account.
3. In the leftmost menu, under the API header select 'API Keys'
4. In the textbox with the 'API Key Name' placehold text, name your key, e.g. esm
5. Create a new Forge API token and copy the token string.

Treat this token as a secret (similar to a password); do not publish
or commit it.

Using the token in the KaML web interface:

1. Open:

       https://kaml.computchem.org

2. Select ESMC for either channel, i.e. Acidic or Basic.
2. Enter your protein sequence, UniProtID, or PDBID-CHAINID into the input box.
3. Paste your ESM Forge API token into the token field.
4. Submit to run predictions.

The Forge token is used only to obtain ESM-C embeddings for your
sequences via the Forge API; KaML then applies the pre-trained KaML-ESMC
heads to produce residue-level predictions.

Please follow EvolutionaryScale’s terms of use and your institution’s
policies when requesting, storing, and using Forge API tokens.

------------------------------------------------------------
Precomputed human proteome-wide atlas: kaml.computchem.org/human
------------------------------------------------------------

The `/human` view provides interactive access to a precomputed
KaML “human atlas” of predictions across the human proteome.

Key points:

- Predictions in this atlas were precomputed using the KaML-ESM2 models and
  ESM-2 embeddings described in the manuscript.
- Exploration of the atlas (searching, browsing, viewing scores) does
  not require an ESM Forge token, because no new embeddings or model
  inference are run client-side.
- The atlas is intended as a convenient starting point for exploring
  residue-level predictions without installing local software.

Typical usage:

1. Open:

       https://kaml.computchem.org/human

2. Search or browse for a protein of interest by UniProtKB accession
   or gene name.
3. Inspect the per-residue KaML-ESM2 predictions. 

The precomputed human proteome-wide atlas integrates three key types of information, which correspond
to sections in the web interface:

# Command-line Inference Interface 
Command-line interface for running KaML-ESM2, KaML-ESMC, and KaML-CBT2 
residue-level predictions on protein sequences using either ESM-2 or ESM-C
embeddings and pre-trained weights.

---

## Repository layout

- kamlCLI.py
  Main command-line interface.

- env/wizard.sh
  Simple installation wizard that creates a Python virtual environment,
  installs dependencies, and downloads pretrained weights from Zenodo.

- env/wts/
  Default location for model weight directories. The CLI expects
  task-specific subdirectories here (i.e. esm2, esmC, CBtree2).

---

## Requirements

- NVIDIA GPU with >= 24GB VRAM (e.g. RTX 4090)
- 40 GB disk space
- 128 GB system memory
- Python 3.10 or newer
- POSIX-like environment (Linux / macOS)
- Packages: $(project_root)/env/KaML-ESM_env.txt

You can install these manually:

    pip install -r $(path_to)/KaML-ESM_env.txt

note: if you do a manual installation, you must add $(project_root)/bin to your PATH variable.

alternatively (recommended): use the provided wizard.

---

## Installation

Clone the repository:

    git clone https://github.com/JanaShenLab/KaML-ESM.git
    cd KaML-ESM

Run the installation wizard:

    bash env/wizard.sh

The wizard:

- creates a virtual environment named "KaML"
- activates it
- upgrades pip
- installs core Python dependencies
- is the place to add commands to download weights from Zenodo into
  env/wts/

To re-activate the environment later:

    source activate

---

## Weights

By default, kamlCLI.py looks for model weights under:

    env/wts/

with task-specific subdirectories (for example):

- env/wts/esm2/acidic
- env/wts/esm2/basic
- env/wts/esmC/acidic
- env/wts/esmC/basic
- env/wts/CBtree2

You can override the root directory for weights with:

    export KAML_WTS_DIR=$(path_to)/wts

and the CLI will use that directory instead of env/wts/.

---

## Forge token

NOTE: This section ONLY pertains to KaML-ESMC models.

The CLI requires an ESM Forge token to compute ESM-C embeddings.

- To obtain an ESM Forge token see Web-based inference.

Token handling:

1. First run:
   - Provide a token via --forge-token, either as the raw string or a
     path to a file that contains the token.
   - The token is cached to a user-level file (by default:
     ~/.esm_forge_token).

2. Subsequent runs:
   - If --forge-token is omitted, the cached token is used.
   - If no cached token exists, the CLI will prompt for one
     interactively (input is hidden) and then cache it.

You can override the cache location with:

    export KAML_FORGE_TOKEN_FILE=$(path_to)/token_file

---

## Basic usage

### Activate the environment

```bash
source env/KaML/bin/activate
```

### CLI help (all options)

```bash
python kamlCLI.py --help
```

### Forge token requirement (important)

If you do **not** pass `--nofold`, KaML may need to fold a structure. Folding uses **ESM Forge**, so you must have a Forge token available (via `--forge-token` or the cached token mechanism) **unless** you provide a structure:

- Provide a structure directly: `--pdb ...` or `--pdbid ...`
- Or provide a precomputed structure directory: `--structs ...` containing `<unique_id>.pdb`

> Practical rule: If you are running on `--seq`, `--uniprot`, or `--fasta` and you **don’t** pass `--nofold`, you should assume a Forge token is required unless you are sure a matching structure will be found via `--structs`.

---

### Inputs (choose exactly one)

- `--seq <AASEQ>` : amino-acid sequence string
- `--uniprot <UNIPROT_ID>` : fetch sequence from UniProt
- `--pdb <PATH>` : use a local PDB file
- `--pdbid <PDB_ID>` : fetch a PDB from RCSB
- `--fasta <PATH>` : process a multi-FASTA file (one output subdir per record)

---

### Common options

- `--outdir <DIR>` : output directory (default: `output/`)
- `--structs <DIR>` : directory of precomputed structures; KaML looks for `<unique_id>.pdb` here before folding
- `--nofold` : skip folding (sequence-only inference when no structure is provided)
- `--nproc <INT>` : parallel workers for multi-FASTA (default: 1)
- `--nthreads <INT>` : threads per sequence for ensemble inference (default: 1)
- `--acidic {esm2,esmC}` : acidic channel model preference (default: `esm2`)
- `--basic  {esm2,esmC}` : basic channel model (default: `esm2`)
- `--nocbtree` : skip CBTree2 predictions
- `--debug` : enable debug logging
- `--skip_safety` : skip Forge safety filter (requires permission)
- `--forge-token <TOKEN_OR_FILE>` : Forge token string or a path to a file containing it (used/cached when Forge is needed)

---

### Examples

#### 1) Single sequence, **no folding** (no Forge required)

```bash
python kamlCLI.py \\
--seq "ACDEFGHIKLMNPQRSTVWY" \\
--nofold \\
--outdir output
```

#### 2) Single sequence, **fold via Forge** (Forge token required)

```bash
python kamlCLI.py \\
--seq "ACDEFGHIKLMNPQRSTVWY" \\
--forge-token /path/to/forge_token.txt \\
--outdir output
```

#### 3) Provide a structure explicitly (no folding; no Forge required)

```bash
python kamlCLI.py \\
--pdb /path/to/structure.pdb \\
--outdir output
```

#### 4) Fetch a structure by PDB ID (no folding; no Forge required)

```bash
python kamlCLI.py \\
--pdbid 1I0E \\
--outdir output
```

#### 5) Multi-FASTA, **no folding** (no Forge required)

```bash
python kamlCLI.py \\
--fasta proteins.fasta \\
--nofold \\
--nproc 4 \\
--outdir output
```

#### 6) Multi-FASTA, **fold via Forge** (Forge token required unless `--structs` covers every record)

```bash
python kamlCLI.py \\
--fasta proteins.fasta \\
--nproc 4 \\
--forge-token /path/to/forge_token.txt \\
--outdir output
```

#### 7) Multi-FASTA using precomputed structures (no folding if every `<unique_id>.pdb` exists)

```bash
python kamlCLI.py \\
--fasta proteins.fasta \\
--structs /path/to/structs_dir \\
--nproc 4 \\
--outdir output
```

#### 8) Prefer ESM-C embeddings (Forge token required)

```bash
python kamlCLI.py \\
--seq "ACDEFGHIKLMNPQRSTVWY" \\
--nofold \\
--basic esmC \\
--acidic esmC \\
--forge-token /path/to/forge_token.txt \\
--outdir output
```

#### 9) Disable CBTree2

```bash
python kamlCLI.py \\
--seq "ACDEFGHIKLMNPQRSTVWY" \\
--nofold \\
--nocbtree \\
--outdir output
```

## Output

KaML writes outputs to an output directory (default: `./output/`). It does **not** print a tab-separated prediction table to standard output.

You control the output location with `--outdir`:

```bash
python kamlCLI.py \\
--seq "MQLKPMEINPEMLNKVLSRLGVAGQWRFVDVLGLEEESLGSVPAPACALLLLFPLTAQHENFRKKQIEELKGQEVSPKVYFMKQTIGNSCGTIGLIHAVANNQDKLGFEDGSVLKQFLSETEKMSPEDRAKCFEKNEAIQAAHDAVAQEGQCRVDDKVNFHFILFNNVDGHLYELDGRMPFPVNHGASSEDTLLKDAAKVCREFTEREQGEVRFSAVALCKAA" \\
--nofold \\
--outdir output
```

### Output files

For a single sequence run, KaML writes:

- `<outdir>/predictions.csv`
- `<outdir>/predicted_structure.pdb` (only if a structure was provided or generated)

For multi-FASTA runs (`--fasta`), KaML creates one subdirectory per sequence under `<outdir>/`, using `<unique_id>` from the FASTA header (first whitespace-delimited token), and writes the same files inside each subdirectory.

### Notes

- If you do **not** pass `--nofold` and you do not provide a structure (`--pdb`, `--pdbid`, or `--structs`), KaML will attempt to fold via ESM Forge (Forge token required).
- At the end of the run, KaML prints a completion message indicating the output directory, e.g.:

```text
KaML run complete. Outputs written to ./output/
```

The first line of the file is a header:

    Residue_ID    Pred_pKa  Pred_Shift     Error  Conf_pKa
    ...

Each subsequent line corresponds to a residue position. Note: Conf_pKa only present when CBTree2 has been invoked.

---

------------------------------------------------------------
# Archived artifacts (immutable, citable)

Pre-trained Ensemble Weights:
- **KaML-ESM2 ensemble weights v1.0.0 (Zenodo)** — DOI: `https://doi.org/10.5281/zenodo.17943825`
- **KaML-ESMC ensemble weights v1.0.0 (Zenodo)** — DOI: `https://doi.org/10.5281/zenodo.17943447`
- **KaML-CBT2 weights v1.0.0 (Zenodo)** - DOI: `https://doi.org/10.5281/zenodo.17943947`
