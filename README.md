# mRNA Target Selection Pipeline

This repository provides a flexible, modular pipeline for extracting, processing, and evaluating protein antigens and epitopes for vaccine target selection, with a focus on pathogens such as *Staphylococcus aureus* and influenza. The pipeline integrates data from IEDB, UniProt, and literature, and supports downstream immunoinformatics analyses.

---

## Features
- **Automated data retrieval** from IEDB and UniProt
- **Antigen compilation** from IEDB, literature, and patents
- **Protein sequence fetching** and random negative control generation
- **Epitope prediction** (MHC-I, MHC-II, B-cell)
- **Antigen/epitope clustering and conservation analysis** (MMseqs2)
- **Antigen feature prediction** (SignalP, TargetP)
- **Epitope evaluation** (Allergenicity, Population Coverage)
- **Statistical analysis** (KS test, feature comparison)
- **Highly modular and parallelizable**

---

## Project Structure
```
mRNA_Target_Selection/
├── bin/
│   ├── IEDB_fetch.py                # Fetch antigens/epitopes from IEDB
│   ├── compile_antigens.py          # Compile antigens from IEDB/literature
│   ├── fetch_sequences_Uniprot.py   # Fetch protein sequences from UniProt
│   ├── generate_random_sequences.py # Fetch random non-antigen proteins
│   ├── align_antigens_mmseqs.py     # Align antigens to genomes (MMseqs2)
│   ├── antigen_analysis.py          # Run SignalP/TargetP on antigens
│   ├── calculate_features_kstest.py # Feature extraction and KS test
│   ├── evaluate_epitopes.py         # Evaluate epitopes (allergenicity, pop coverage)
│   ├── sanity_check_antigen_seq.py  # Validate epitope mapping to antigens
│   └── tools/
│       ├── common.py                # Shared utilities/constants
│       ├── extract_epitopes.py      # Epitope extraction from predictions
│       ├── run_algpred.py           # Allergenicity prediction (AlgPred2)
│       ├── run_bcell.py             # B-cell epitope prediction
│       ├── run_cluster.py           # MMseqs2 clustering
│       ├── run_mhci.py              # MHC-I epitope prediction
│       ├── run_mhcii.py             # MHC-II epitope prediction
│       ├── run_popcoverage.py       # Population coverage analysis
│       ├── run_signalp.py           # SignalP wrapper
│       └── run_targetp.py           # TargetP wrapper
├── data/
│   └── <pathogen_subfolders>/
│       ├── <organism_tag>_IEDB_antigens.csv
│       ├── <organism_tag>_IEDB_epitope.csv
│       ├── <organism_tag>_compiled_antigens.csv
│       ├── <organism_tag>_compiled_proteins.csv
│       ├── random_compiled_proteins.csv
│       └── ... (outputs, features, results)
├── notes/                           # Meeting notes, literature, questions
├── requirements.txt                 # Python dependencies
├── algpred2_dependencies.yml        # Conda env for AlgPred2
└── README.md                        # This file
```

---

## Pipeline Overview

1. **Data Collection**
    - `IEDB_fetch.py`: Download antigen and epitope data from IEDB.
    - `compile_antigens.py`: Merge IEDB and literature antigens.
    - `fetch_sequences_Uniprot.py`: Fetch protein sequences for antigens.
    - `generate_random_sequences.py`: Fetch random non-antigen proteins for controls.

2. **Sequence Analysis**
    - `align_antigens_mmseqs.py`: Align antigens to strain genomes (MMseqs2).
    - `antigen_analysis.py`: Run SignalP/TargetP on antigens.
    - `sanity_check_antigen_seq.py`: Validate epitope mapping to antigens.

3. **Epitope Prediction**
    - `tools/run_mhci.py`, `tools/run_mhcii.py`: Predict MHC-I/II epitopes.
    - `tools/run_bcell.py`: Predict B-cell epitopes.
    - `tools/run_cluster.py`: Cluster proteins/epitopes (MMseqs2).

4. **Epitope Evaluation**
    - `evaluate_epitopes.py`: Run allergenicity (AlgPred2) and population coverage.
    - `tools/extract_epitopes.py`: Extract and filter predicted epitopes.
    - `calculate_features_kstest.py`: Extract features and compare positive/negative sets (KS test).

5. **Utilities**
    - `tools/common.py`: Directory, file, and allele panel utilities.

---

## Dependencies

- Python 3.7+
- pandas, requests, scipy, numpy, scikit-learn, joblib
- [NCBI Entrez Direct (EDirect)](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [Seqkit](https://bioinf.shenwei.me/seqkit/usage/#translate)
- [MMseqs2](https://github.com/soedinglab/MMseqs2) (via conda)
- [AlgPred2](https://github.com/masashitsuboi/AlgPred2) (via conda/pip, see `algpred2_dependencies.yml`)
- [SignalP](https://services.healthtech.dtu.dk/services/SignalP-5.0/)
- [TargetP](https://services.healthtech.dtu.dk/services/TargetP-2.0/)
- [IEDB MHCI Epitope Prediction](https://nextgen-tools.iedb.org/pipeline?tool=tc1)
- [IEDB MHCII Epitope Prediction](https://nextgen-tools.iedb.org/pipeline?tool=tc2)
- [IEDB Bcell Epitope Prediction](http://tools.iedb.org/bcell/)
- [IEDB Population Coverage](http://tools.iedb.org/population/)
- [DeepLocPro1.0](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/)
- [Convert_CIF_to_PDB](https://github.com/SDMscript/Convert_CIF_to_PDB)

Install Python dependencies:
```sh
pip install -r requirements.txt
```

Install EDirect:
```sh
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
# or
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
export PATH=${HOME}/edirect:${PATH}
```

Install Seqkit:
```sh
conda install -c bioconda seqkit
```

Install MMseqs2: See [MMseqs2 GitHub](https://github.com/soedinglab/MMseqs2)
```sh
# install via conda
conda install -c conda-forge -c bioconda mmseqs2
```

Install DeeplocPro: See [Deeplocpro GitHub](https://github.com/Jaimomar99/deeplocpro)
```sh
git clone https://github.com/Jaimomar99/deeplocpro.git
cd DeepLocPro
pip install .

```

Install CIF to PDB converter: See [Github](https://github.com/SDMscript/Convert_CIF_to_PDB)
```sh
git clone https://github.com/SDMscript/Convert_CIF_to_PDB.git
```
---

## Usage

Each script is self-documented and can be run as a standalone tool. See the docstring at the top of each script for arguments and usage examples. Typical workflow:

1. Fetch and compile antigens/epitopes:
    ```sh
    python bin/IEDB_fetch.py s_aureus "Staphylococcus aureus"
    python bin/compile_antigens.py s_aureus "Staphylococcus aureus"
    python bin/fetch_sequences_Uniprot.py s_aureus "Staphylococcus aureus"
    python bin/generate_random_sequences.py s_aureus "Staphylococcus aureus"
    ```
2. Predict and analyze features/epitopes:
    ```sh
    python bin/antigen_analysis.py s_aureus proteins --tool-root bin/tools
    python bin/evaluate_epitopes.py s_aureus epitope_outputs --tool-root bin/tools
    python bin/calculate_features_kstest.py s_aureus
    ```

---

## Notes
- See `notes/` for meeting minutes, literature, and project questions.
- All scripts are modular and can be adapted for other pathogens or datasets.
- For details on each tool, see the docstring in the corresponding script.

---

## Citation
If you use this pipeline, please cite the IEDB and relevant tool authors as appropriate.