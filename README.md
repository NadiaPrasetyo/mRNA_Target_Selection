# mRNA Target Selection Pipeline

This repository provides a comprehensive, modular pipeline for extracting, processing, and evaluating protein antigens and epitopes for vaccine target selection, with a focus on pathogens such as *Staphylococcus aureus* and influenza. The pipeline integrates data from IEDB, UniProt, literature, and patents, and supports downstream immunoinformatics analyses with statistical evaluation.

---

## Features

- **Automated data retrieval** from IEDB, UniProt, and literature sources
- **Comprehensive antigen compilation** from multiple data sources including patents
- **Multi-strain genome analysis** with automated sequence fetching and translation
- **Advanced epitope prediction** (MHC-I, MHC-II, B-cell) using state-of-the-art tools
- **Protein sequence clustering and conservation analysis** using MMseqs2 and Rate4Site
- **Subcellular localization prediction** (SignalP, TargetP, DeepLocPro, DeepTMHMM)
- **Immunological feature evaluation** (allergenicity, antigenicity, population coverage)
- **Phylogenetic analysis** with MAFFT and evolutionary rate estimation
- **Statistical comparison** between antigen and random protein sets using KS-tests and AUROC
- **Visualization and reporting** with ROC curves and feature distribution plots
- **Highly modular and parallelizable** architecture

---

## Project Structure

```
mRNA_Target_Selection/
├── bin/
│   ├── IEDB_fetch.py                # Fetch antigens/epitopes from IEDB
│   ├── IEDB_epitope.py              # IEDB epitope prediction pipeline
│   ├── compile_antigens.py          # Compile antigens from multiple sources
│   ├── fetch_sequences_Uniprot.py   # Fetch protein sequences from UniProt
│   ├── generate_random_sequences.py # Generate random non-antigen controls
│   ├── align_antigens_mmseqs.py     # Align antigens to strain genomes
│   ├── antigen_analysis.py          # Comprehensive antigen feature analysis
│   ├── calculate_features_kstest.py # Statistical feature comparison
│   ├── evaluate_epitopes.py         # Epitope evaluation pipeline
│   ├── sanity_check_antigen_seq.py  # Validate epitope-antigen mapping
│   ├── fetch_NCBI_strain_genome.sh  # Fetch and translate strain genomes
│   ├── fetch_PDB_structure.py       # Fetch protein structures from PDB
│   ├── parse_patent_to_fasta.py     # Parse patent sequences to FASTA
│   ├── visualise_cluster.py         # Visualize clustering results
│   ├── R/
│   │   ├── R.Rproj                  # R project file
│   │   └── visualize_features.R     # Statistical visualization in R
│   └── tools/
│       ├── common.py                # Shared utilities and constants
│       ├── extract_epitopes.py      # Extract epitopes from predictions
│       ├── run_algpred.py           # Allergenicity prediction (AlgPred2)
│       ├── run_bcell.py             # B-cell epitope prediction
│       ├── run_cluster.py           # MMseqs2 clustering
│       ├── run_deeplocpro.py        # Subcellular localization (DeepLocPro)
│       ├── run_deeptmhmm.py         # Transmembrane topology (DeepTMHMM)
│       ├── run_ellipro.py           # ElliPro B-cell epitope prediction
│       ├── run_ifnepitope2.py       # IFN-gamma epitope prediction
│       ├── run_mafft.py             # Multiple sequence alignment
│       ├── run_mhci.py              # MHC-I epitope prediction
│       ├── run_mhcii.py             # MHC-II epitope prediction
│       ├── run_mixmhc2pred.py       # MixMHC2pred epitope prediction
│       ├── run_popcoverage.py       # Population coverage analysis
│       ├── run_rate4site.py         # Evolutionary rate estimation
│       ├── run_signalp.py           # Signal peptide prediction
│       └── run_targetp.py           # Protein targeting prediction
├── data/
│   ├── literature-sourced/          # Literature-derived antigen data
│   ├── patents/                     # Patent-derived sequences
│   └── <pathogen_subfolders>/       # Pathogen-specific data
│       ├── <organism_tag>_IEDB_antigens.csv
│       ├── <organism_tag>_IEDB_epitope.csv
│       ├── <organism_tag>_compiled_antigens.csv
│       ├── <organism_tag>_compiled_proteins.csv
│       ├── random_compiled_proteins.csv
│       ├── strain_genomes/          # Strain genome sequences
│       ├── epitope_outputs/         # Prediction results
│       │   ├── mhci/, mhcii/, bcell/, ellipro/
│       │   ├── algpred/, ifnepitope2/, mixmhc2pred/
│       │   ├── signalp/, targetp/, deeplocpro/, deeptmhmm/
│       │   ├── cluster/, mafft_rate4site/
│       │   └── evaluation_outputs/
│       ├── random_analysis/         # Random protein analysis
│       ├── mmseqs_results/          # Sequence alignment results
│       └── pdb_results/             # PDB structure data
├── results/
│   └── <pathogen_subfolders>/       # Analysis results and visualizations
│       ├── ks_test_results.csv      # Statistical comparison results
│       ├── roc_plots/               # ROC curve visualizations
│       └── feature_distributions/   # Feature distribution plots
├── notes/                           # Documentation and meeting notes
├── requirements.txt                 # Python dependencies
├── ext_tools_dependencies.yml       # Conda environment definition
└── README.md                        # This file
```

---

## Pipeline Overview

### 1. Data Collection
- **IEDB Integration**: Automated fetching of experimental epitope and antigen data
- **Literature Mining**: Compilation of antigens from research publications
- **Patent Analysis**: Extraction of sequences from patent databases
- **Strain Genome Retrieval**: Multi-strain genome fetching and translation
- **Random Controls**: Generation of matched random protein sets

### 2. Sequence Analysis
- **Homology Mapping**: Antigen alignment to strain genomes using MMseqs2
- **Clustering Analysis**: Protein sequence clustering and similarity scoring
- **Conservation Analysis**: Evolutionary rate estimation with Rate4Site
- **Phylogenetic Analysis**: Multiple sequence alignment with MAFFT

### 3. Feature Prediction
- **Subcellular Localization**: SignalP, TargetP, DeepLocPro predictions
- **Membrane Topology**: DeepTMHMM transmembrane region prediction
- **Allergenicity Assessment**: AlgPred2 allergenicity prediction
- **Antigenicity Prediction**: IFN-gamma epitope prediction with IFNepitope2

### 4. Epitope Prediction
- **MHC-I/II Binding**: NetMHCpan and MixMHC2pred predictions
- **B-cell Epitopes**: Bepipred and ElliPro conformational predictions
- **Population Coverage**: HLA allele frequency analysis
- **Epitope Filtering**: Threshold-based epitope selection

### 5. Statistical Analysis
- **Feature Comparison**: KS-test analysis between antigen and random sets
- **Performance Metrics**: AUROC calculation for all features
- **Visualization**: ROC curves and distribution plots
- **Feature Ranking**: Identification of discriminative features

---

## Dependencies

### Core Requirements
- Python 3.7+
- pandas, requests, scipy, numpy, scikit-learn, matplotlib, networkx, tqdm
- biopython, gemmi, pybiolib, torch

### External Tools
- [MMseqs2](https://github.com/soedinglab/MMseqs2) (via conda)
- [NCBI EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [Seqkit](https://bioinf.shenwei.me/seqkit/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (via conda)
- [Rate4Site](https://www.tau.ac.il/~itaymay/cp/rate4site.html) (via conda)

### Specialized Tools
- [SignalP](https://services.healthtech.dtu.dk/services/SignalP-5.0/)
- [TargetP](https://services.healthtech.dtu.dk/services/TargetP-2.0/)
- [DeepLocPro](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/)
- [DeepTMHMM](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/)
- [AlgPred2](https://github.com/masashitsuboi/AlgPred2) (via conda)
- [IFNepitope2](https://webs.iiitd.edu.in/raghava/ifnepitope2/) (via pip)
- [MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred)

### IEDB Tools
- [IEDB MHCI/II Prediction](https://nextgen-tools.iedb.org/)
- [IEDB B-cell Prediction](http://tools.iedb.org/bcell/)
- [IEDB Population Coverage](http://tools.iedb.org/population/)
- [IEDB ElliPro](http://tools.iedb.org/ellipro/)

### Installation

Install Python dependencies:
```bash
pip install -r requirements.txt
```

Create external tools environment:
```bash
conda env create -f ext_tools_dependencies.yml
conda activate external_tools_env
```

Install core bioinformatics tools:
```bash
# EDirect
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=${HOME}/edirect:${PATH}

# Seqkit and MMseqs2
conda install -c bioconda seqkit mmseqs2

# DeepLocPro
git clone https://github.com/Jaimomar99/deeplocpro.git
cd deeplocpro && pip install .
```

Install R and visualization packages:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr"))
```

---

## Usage

### Complete Workflow Example

```bash
# 1. Fetch and compile data
python bin/IEDB_fetch.py s_aureus "Staphylococcus aureus"
python bin/IEDB_epitope.py s_aureus "Staphylococcus aureus"
python bin/compile_antigens.py s_aureus "Staphylococcus aureus"
python bin/fetch_sequences_Uniprot.py s_aureus "Staphylococcus aureus"
python bin/generate_random_sequences.py s_aureus "Staphylococcus aureus"

# 2. Fetch strain genomes and structures
bash bin/fetch_NCBI_strain_genome.sh s_aureus S.aureus_strains.csv
python bin/fetch_PDB_structure.py s_aureus proteins --threads 8

# 3. Run comprehensive analysis
python bin/antigen_analysis.py s_aureus proteins --tool-root /path/to/tools --threads 8
python bin/antigen_analysis.py s_aureus random_proteins --tool-root /path/to/tools --threads 8

# 4. Epitope prediction and evaluation
python bin/IEDB_epitope.py s_aureus proteins --tool-root /path/to/iedb_tools --threads 8
python bin/evaluate_epitopes.py s_aureus epitope_outputs --tool-root /path/to/tools

# 5. Statistical analysis and visualization
python bin/calculate_features_kstest.py s_aureus --threads 4 --verbose --write-raw
python bin/visualise_cluster.py s_aureus epitope_outputs --split-clusters
Rscript bin/R/visualize_features.R
```

### Key Analysis Results

The pipeline generates comprehensive results including:

- **Feature discrimination**: KS-test p-values and AUROC scores for 50+ features
- **ROC curves**: Performance visualization for each feature
- **Conservation scores**: Evolutionary analysis across multiple strains
- **Population coverage**: HLA allele frequency analysis
- **Clustering results**: Sequence similarity and conservation mapping

---

## Case Study: *Staphylococcus aureus*

The pipeline has been extensively tested on *S. aureus* data, analyzing:
- **245 antigen proteins** from IEDB and literature
- **359 random control proteins** from UniProt
- **6 representative strain genomes** (CC1, CC5, CC8, CC22, CC30, CC93)
- **Statistical comparison** across multiple immunological features

Key findings demonstrate the pipeline's ability to identify discriminative features for vaccine target selection.

---

## Configuration

The pipeline supports extensive configuration through:
- **Allele panels**: Default and extended HLA allele sets
- **Thresholds**: Customizable IC50, percentile, and similarity cutoffs
- **Tool selection**: Modular tool execution with configurable parameters
- **Output formats**: JSON, CSV, FASTA, and visualization outputs

---

## Notes

- See [`notes/`](notes/) for meeting minutes, literature references, and project documentation
- All scripts are modular and can be adapted for other pathogens or datasets
- The pipeline supports both academic and research use cases
- Extensive logging and error handling for robust operation

---

## Citation

If you use this pipeline, please cite the relevant tool authors:
- IEDB: Vita et al. (2019)
- MMseqs2: Steinegger & Söding (2017)
- MAFFT: Katoh & Standley (2013)
- And other tools as appropriate for your analysis

---

## Contact

For questions or contributions, please refer to the project documentation in [`notes/`](notes/) or contact the development team.