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
│   ├── IEDB_epitope.py              # Additional IEDB epitope processing
│   ├── compile_antigens.py          # Compile antigens from IEDB/literature
│   ├── fetch_sequences_Uniprot.py   # Fetch protein sequences from UniProt
│   ├── generate_random_sequences.py # Fetch random non-antigen proteins
│   ├── align_antigens_mmseqs.py     # Align antigens to genomes (MMseqs2)
│   ├── antigen_analysis.py          # Run SignalP/TargetP on antigens
│   ├── calculate_features_kstest.py # Feature extraction and KS test
│   ├── evaluate_epitopes.py         # Evaluate epitopes (allergenicity, pop coverage)
│   ├── sanity_check_antigen_seq.py  # Validate epitope mapping to antigens
│   ├── fetch_NCBI_strain_genome.sh  # Fetch strain genomes from NCBI
│   ├── fetch_PDB_structure.py       # Fetch protein structures from PDB
│   ├── parse_patent_to_fasta.py     # Parse patent sequences to FASTA format
│   ├── visualise_cluster.py         # Visualize clustering results
│   ├── R/
│   │   ├── R.Rproj                  # R project file
│   │   └── visualize_features.R     # R script for feature visualization
│   └── tools/
│       ├── __init__.py              # Python package initialization
│       ├── common.py                # Shared utilities/constants
│       ├── extract_epitopes.py      # Epitope extraction from predictions
│       ├── run_algpred.py           # Allergenicity prediction (AlgPred2)
│       ├── run_bcell.py             # B-cell epitope prediction
│       ├── run_cluster.py           # MMseqs2 clustering
│       ├── run_deeplocpro.py        # Subcellular localization prediction
│       ├── run_ellipro.py           # ElliPro B-cell epitope prediction
│       ├── run_ifnepitope2.py       # IFN-gamma epitope prediction
│       ├── run_mhci.py              # MHC-I epitope prediction
│       ├── run_mhcii.py             # MHC-II epitope prediction
│       ├── run_mixmhc2pred.py       # MixMHC2pred epitope prediction
│       ├── run_popcoverage.py       # Population coverage analysis
│       ├── run_signalp.py           # SignalP wrapper
│       └── run_targetp.py           # TargetP wrapper
├── data/
│   ├── literature-sourced/          # Literature-derived antigen data
│   ├── patents/                     # Patent-derived sequences and data
│   └── <pathogen_subfolders>/       # Pathogen-specific data folders
│       ├── <organism_tag>_IEDB_antigens.csv
│       ├── <organism_tag>_IEDB_epitope.csv
│       ├── <organism_tag>_compiled_antigens.csv
│       ├── <organism_tag>_compiled_proteins.csv
│       ├── random_compiled_proteins.csv
│       ├── epitope_outputs/         # Epitope prediction results
│       ├── evaluation_outputs/      # Epitope evaluation results
│       ├── mmseqs_results/          # MMseqs2 clustering results
│       ├── pdb_results/             # PDB structure data
│       ├── strain_genomes/          # Strain genome sequences
│       └── ... (features, analysis results)
├── notes/                           # Meeting notes, literature, questions
│   ├── external_resources/          # External documentation and resources
│   ├── literature_search/           # Literature search notes
│   └── meetings/                    # Meeting minutes and discussions
├── results/                         # Analysis results and visualizations
│   └── <pathogen_subfolders>/       # Pathogen-specific results
├── requirements.txt                 # Python dependencies
├── algpred2_dependencies.yml        # Conda env for AlgPred2 and IFNepitope2
└── README.md                        # This file
```

---

## Pipeline Overview

1. **Data Collection**
    - `IEDB_fetch.py`: Download antigen and epitope data from IEDB.
    - `IEDB_epitope.py`: Additional IEDB epitope processing and analysis.
    - `compile_antigens.py`: Merge IEDB and literature antigens.
    - `fetch_sequences_Uniprot.py`: Fetch protein sequences for antigens.
    - `generate_random_sequences.py`: Fetch random non-antigen proteins for controls.
    - `fetch_NCBI_strain_genome.sh`: Fetch strain genomes from NCBI databases.
    - `fetch_PDB_structure.py`: Retrieve protein structures from PDB.
    - `parse_patent_to_fasta.py`: Parse patent sequences to FASTA format.

2. **Sequence Analysis**
    - `align_antigens_mmseqs.py`: Align antigens to strain genomes (MMseqs2).
    - `antigen_analysis.py`: Run SignalP/TargetP on antigens.
    - `sanity_check_antigen_seq.py`: Validate epitope mapping to antigens.
    - `visualise_cluster.py`: Visualize clustering and alignment results.

3. **Epitope Prediction**
    - `tools/run_mhci.py`, `tools/run_mhcii.py`: Predict MHC-I/II epitopes.
    - `tools/run_bcell.py`: Predict B-cell epitopes.
    - `tools/run_ellipro.py`: ElliPro B-cell epitope prediction.
    - `tools/run_mixmhc2pred.py`: MixMHC2pred epitope prediction.
    - `tools/run_ifnepitope2.py`: IFN-gamma epitope prediction.
    - `tools/run_cluster.py`: Cluster proteins/epitopes (MMseqs2).

4. **Epitope Evaluation**
    - `evaluate_epitopes.py`: Run allergenicity (AlgPred2) and population coverage.
    - `tools/run_algpred.py`: Allergenicity prediction using AlgPred2.
    - `tools/run_popcoverage.py`: Population coverage analysis.
    - `tools/extract_epitopes.py`: Extract and filter predicted epitopes.
    - `calculate_features_kstest.py`: Extract features and compare positive/negative sets (KS test).

5. **Additional Analysis Tools**
    - `tools/run_deeplocpro.py`: Subcellular localization prediction.
    - `tools/run_signalp.py`: Signal peptide prediction.
    - `tools/run_targetp.py`: Protein targeting prediction.
    - `R/visualize_features.R`: R script for statistical visualization of features.

6. **Utilities**
    - `tools/common.py`: Directory, file, and allele panel utilities.

---

## Dependencies

- Python 3.7+
- pandas, requests, scipy, numpy, scikit-learn, joblib, matplotlib, networkx, tqdm, pybiolib, torch, biopython, gemmi
- [NCBI Entrez Direct (EDirect)](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [Seqkit](https://bioinf.shenwei.me/seqkit/usage/#translate)
- [MMseqs2](https://github.com/soedinglab/MMseqs2) (via conda)
- [AlgPred2](https://github.com/masashitsuboi/AlgPred2) (via conda/pip, see `algpred2_dependencies.yml`)
- [IFNepitope2](https://webs.iiitd.edu.in/raghava/ifnepitope2/) (via pip)
- [SignalP](https://services.healthtech.dtu.dk/services/SignalP-5.0/)
- [TargetP](https://services.healthtech.dtu.dk/services/TargetP-2.0/)
- [DeepLocPro1.0](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/)
- [IEDB MHCI Epitope Prediction](https://nextgen-tools.iedb.org/pipeline?tool=tc1)
- [IEDB MHCII Epitope Prediction](https://nextgen-tools.iedb.org/pipeline?tool=tc2)
- [IEDB Bcell Epitope Prediction](http://tools.iedb.org/bcell/)
- [IEDB Population Coverage](http://tools.iedb.org/population/)
- [IEDB Ellipro](http://tools.iedb.org/ellipro/)
- [MixMHC2pred-2.0](https://github.com/GfellerLab/MixMHC2pred)
- R (for statistical visualization)

Install Python dependencies:
```sh
pip install -r requirements.txt
```

Install AlgPred2 and IFNepitope2 environment:
```sh
conda env create -f algpred2_dependencies.yml
conda activate algpred2_env
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

Install MixMHC2pred-2.0 in [MixMHC2pred GitHub](https://github.com/GfellerLab/MixMHC2pred/releases)

Install R and required packages:
```r
# In R console
install.packages(c("ggplot2", "dplyr", "tidyr", "readr"))
```

## Usage

Each script is self-documented and can be run as a standalone tool. See the docstring at the top of each script for arguments and usage examples. Typical workflow:

1. Fetch and compile antigens/epitopes:
    ```sh
    python bin/IEDB_fetch.py s_aureus "Staphylococcus aureus"
    python bin/IEDB_epitope.py s_aureus "Staphylococcus aureus"
    python bin/compile_antigens.py s_aureus "Staphylococcus aureus"
    python bin/fetch_sequences_Uniprot.py s_aureus "Staphylococcus aureus"
    python bin/generate_random_sequences.py s_aureus "Staphylococcus aureus"
    ```

2. Fetch additional data (optional):
    ```sh
    bash bin/fetch_NCBI_strain_genome.sh s_aureus strain_list.csv
    python bin/fetch_PDB_structure.py s_aureus
    python bin/parse_patent_to_fasta.py patents/sequences/
    ```

3. Predict and analyze features/epitopes:
    ```sh
    python bin/antigen_analysis.py s_aureus proteins --tool-root bin/tools
    python bin/evaluate_epitopes.py s_aureus epitope_outputs --tool-root bin/tools
    python bin/calculate_features_kstest.py s_aureus
    ```

4. Visualization and clustering:
    ```sh
    python bin/visualise_cluster.py s_aureus
    Rscript bin/R/visualize_features.R  # Run from project root
    ```

---

## Notes
- See `notes/` for meeting minutes, literature, and project questions.
- All scripts are modular and can be adapted for other pathogens or datasets.
- For details on each tool, see the docstring in the corresponding script.

---

## Citation
If you use this pipeline, please cite the IEDB and relevant tool authors as appropriate.