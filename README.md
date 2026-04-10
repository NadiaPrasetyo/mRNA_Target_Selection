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
- **Statistical comparison** between antigen and random/human protein sets using KS-tests and AUROC
- **Visualization and reporting** with ROC curves and feature distribution plots
- **Highly modular and parallelizable** architecture
- **Machine learning analysis** with Random Forest and PCA for antigen prediction
- **Complete pipeline orchestration** with master run script

## Categories

### Subcellular Localisation
- Predicts protein localization using tools like SignalP, TargetP, and DeepLocPro.
- Identifies membrane-bound and secreted proteins for vaccine target prioritization.

### Allergenicity
- Evaluates potential allergenicity of antigens using AlgPred2.
- Filters out proteins with high allergenic potential to ensure safety.

### Immunogenicity
- Assesses antigenicity and population coverage using tools like IFNepitope2 and IEDB Population Coverage.
- Prioritizes proteins with broad immunogenic potential across diverse HLA alleles.

### Conservation Analysis
- Analyzes sequence conservation across multiple strains using MMseqs2 and Rate4Site.
- Identifies conserved regions critical for vaccine design.

### Epitope Prediction
- Predicts MHC-I, MHC-II, and B-cell epitopes using tools like NetMHCpan, MixMHC2pred, and ElliPro.
- Filters epitopes based on binding affinity, population coverage, and immunogenicity.

### Structure/Composition Analysis
- Analyzes protein structure and composition using PDB data and Pfam domain annotations.
- Identifies structural features relevant for antigenicity and vaccine design.

---

## Project Structure

```
mRNA_Target_Selection/
├── bin/
│   ├── run_whole_pipeline.py       # Master pipeline orchestrator
│   ├── IEDB_fetch.py              # Fetch antigens/epitopes from IEDB
│   ├── IEDB_epitope.py            # IEDB epitope prediction pipeline
│   ├── compile_antigens.py        # Compile antigens from multiple sources
│   ├── fetch_sequences_Uniprot.py # Fetch protein sequences from UniProt
│   ├── generate_random_sequences.py # Generate random/human non-antigen controls
│   ├── fetch_NCBI_strain_genomes.py # Fetch and translate strain genomes (Python version)
│   ├── align_antigens_mmseqs.py   # Align antigens to strain genomes
│   ├── antigen_analysis.py        # Comprehensive antigen feature analysis
│   ├── calculate_features_kstest.py # Statistical feature comparison
│   ├── evaluate_epitopes.py       # Epitope evaluation pipeline
│   ├── sanity_check_antigen_seq.py # Validate epitope-antigen mapping
│   ├── fetch_PDB_structure.py     # Fetch protein structures from PDB
│   ├── fetch_pfam_hmmer.py        # Fetch Pfam HMM profiles
│   ├── parse_patent_to_fasta.py   # Parse patent sequences to FASTA
│   ├── visualise_cluster.py       # Visualize clustering results
│   ├── rf_antigen_pipeline.py     # Random Forest antigen prediction
│   ├── pca_analysis.py           # PCA and statistical analysis
│   ├── R/
│   │   ├── R.Rproj               # R project file
│   │   └── visualize_features.R  # Statistical visualization in R
│   └── tools/
│       ├── common.py             # Shared utilities and constants
│       ├── run_algpred.py        # Allergenicity prediction (AlgPred2)
│       ├── run_bcell.py          # B-cell epitope prediction (BepiPred-3.0)
│       ├── run_cluster.py        # MMseqs2 clustering
│       ├── run_deeplocpro.py     # Subcellular localization (DeepLocPro)
│       ├── run_deeptmhmm.py      # Transmembrane topology (DeepTMHMM)
│       ├── run_discotope.py      # DiscoTope B-cell epitope prediction
│       ├── run_dnds.py           # dN/dS analysis with MACSE/HyPhy
│       ├── run_dssp.py           # Secondary structure (PyDSSP)
│       ├── run_ellipro.py        # ElliPro B-cell epitope prediction
│       ├── run_ifnepitope2.py    # IFN-gamma epitope prediction
│       ├── run_mafft.py          # Multiple sequence alignment
│       ├── run_mhci.py           # MHC-I epitope prediction (NetMHCpan-4.2)
│       ├── run_mhcii.py          # MHC-II epitope prediction (NetMHCIIpan-4.3)
│       ├── run_mixmhc2pred.py    # MixMHC2pred epitope prediction
│       ├── run_protlearn.py      # Protein feature extraction
│       ├── run_rate4site.py      # Evolutionary rate estimation
│       ├── run_signalp.py        # Signal peptide prediction
│       └── run_targetp.py        # Protein targeting prediction
├── data/
│   ├── literature-sourced/        # Literature-derived antigen data
│   ├── patents/                   # Patent-derived sequences
│   └── <pathogen_subfolders>/     # Pathogen-specific data
│       ├── <organism_tag>_IEDB_antigens.csv
│       ├── <organism_tag>_IEDB_epitope.csv
│       ├── <organism_tag>_compiled_antigens.csv
│       ├── <organism_tag>_compiled_proteins.csv
│       ├── random_compiled_proteins.csv
│       ├── human_compiled_proteins.csv
│       ├── strain_genomes/        # Strain genome sequences
│       ├── epitope_outputs/       # Prediction results
│       │   ├── mhci/, mhcii/, bcell/, ellipro/
│       │   ├── algpred/, ifnepitope2/, mixmhc2pred/
│       │   ├── signalp/, targetp/, deeplocpro/, deeptmhmm/
│       │   ├── cluster/, mafft_rate4site/
│       │   ├── dssp/, protlearn/, discotope/
│       │   └── evaluation_outputs/
│       ├── random_analysis/       # Random protein analysis
│       ├── human_analysis/        # Human protein analysis
│       ├── mmseqs_protein/        # Protein sequence alignment results
│       ├── mmseqs_nucleotide/     # Nucleotide sequence alignment results
│       ├── pdb_sequences/         # PDB structure data
│       └── pfam_hmms/            # Pfam HMM profiles
├── results/
│   └── <pathogen_subfolders>/     # Analysis results and visualizations
│       ├── ks_test_results_random.csv # Statistical comparison results (random)
│       ├── ks_test_results_human.csv  # Statistical comparison results (human)
│       ├── roc_plots/             # ROC curve visualizations
│       ├── feature_distributions/ # Feature distribution plots
│       └── raw_data/             # Raw feature data for ML analysis
├── notes/                         # Documentation and meeting notes
├── requirements.txt               # Python dependencies
├── ext_tools_dependencies.yml     # Conda environment definition
├── discotope_tools_dependencies.yml # DiscoTope-specific environment
└── README.md                     # This file
```

---

## Pipeline Overview

### 1. Data Collection
- **IEDB Integration**: Automated fetching of experimental epitope and antigen data
- **Literature Mining**: Compilation of antigens from research publications
- **Patent Analysis**: Extraction of sequences from patent databases
- **Strain Genome Retrieval**: Multi-strain genome fetching and translation
- **Random/Human Controls**: Generation of matched control protein sets

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
- **B-cell Epitopes**: BepiPred-3.0 and ElliPro conformational predictions
- **Discontinuous Epitopes**: DiscoTope conformational B-cell epitope prediction
- **Population Coverage**: HLA allele frequency analysis
- **Epitope Filtering**: Threshold-based epitope selection

### 5. Statistical Analysis
- **Feature Comparison**: KS-test analysis between antigen and random/human sets
- **Performance Metrics**: AUROC calculation for all features
- **Visualization**: ROC curves and distribution plots
- **Feature Ranking**: Identification of discriminative features
- **Machine Learning**: Random Forest and PCA analysis for antigen prediction

---

## Dependencies

### Core Requirements
- Python 3.8+
- pandas, requests, scipy, numpy, scikit-learn, matplotlib, networkx, tqdm
- biopython, gemmi, pybiolib, torch, pydssp
- adjustText (for PCA plots)

### External Tools
- [MMseqs2](https://github.com/soedinglab/MMseqs2) (via conda)
- [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (via conda)
- [Rate4Site](https://www.tau.ac.il/~itaymay/cp/rate4site.html) (via conda)
- [HMMER](http://hmmer.org/) (via conda)
- [FastTree](http://www.microbesonline.org/fasttree/) (via conda)
- [MACSE](https://bioweb.supagro.inra.fr/macse/) (via conda)
- [HyPhy](http://hyphy.org/) (via conda)
- [Pfam database](https://www.ebi.ac.uk/interpro/download/pfam/)
- [EMBOSS backtranseq](https://emboss.bioinformatics.nl/cgi-bin/emboss/backtranseq)

### Specialized Tools
- [SignalP-5.0](https://services.healthtech.dtu.dk/services/SignalP-5.0/)
- [TargetP-2.0](https://services.healthtech.dtu.dk/services/TargetP-2.0/)
- [DeepLocPro](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/)
- [DeepTMHMM](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/) (via BioLib)
- [AlgPred2](https://github.com/masashitsuboi/AlgPred2) (via conda)
- [IFNepitope2](https://webs.iiitd.edu.in/raghava/ifnepitope2/) (via pip)
- [MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred)
- [BepiPred-3.0](https://services.healthtech.dtu.dk/services/BepiPred-3.0/)
- [DiscoTope-3.0](https://github.com/mnielLab/discotope3_web)
- [NetMHCpan-4.2](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/)
- [NetMHCIIpan-4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/)
- [IEDB B-cell Prediction](https://services.healthtech.dtu.dk/services/BepiPred-3.0/)
- [IEDB ElliPro](http://tools.iedb.org/ellipro/)

### Python Libraries for Structure Analysis
- [PyDSSP](https://github.com/ShintaroMinami/PyDSSP) (secondary structure)
- [ProtLearn](https://github.com/tadorfer/protlearn) (protein features)
- [Gemmi](https://github.com/project-gemmi/gemmi) (structure parsing)

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

Create DiscoTope environment (if using DiscoTope):
```bash
conda env create -f discotope_tools_dependencies.yml
conda activate discotope_tools_env
```
### Install Core Bioinformatics Tools

#### Required Tools and Resources
Ensure the following tools and files are available in the `tool-root` directory:

- `BepiPred3_src`
- `deeplocpro`
- `discotope3_web`
- `MixMHC2pred-2.0`
- `netMHCIIpan-4.3`
- `netMHCpan-4.2`
- `signalp-5.0b`
- `targetp-2.0`
- `ElliPro.jar`
- `Pfam-A.hmm`

#### Download Instructions
- **DTU Tools**: Download from their respective websites.
- **MixMHC2pred**: Download the `MixMHC2pred-2.0.zip` file from [GitHub Releases](https://github.com/GfellerLab/MixMHC2pred/releases).
- **ElliPro**: Download from [IEDB Tools](http://tools.iedb.org/ellipro/download/).
- **DeepLocPro**: Clone repository at [github](https://github.com/Jaimomar99/deeplocpro.git)
- Other tools are available as conda/pip packages

#### Install Additional Tools via Conda
```bash
conda install -c bioconda mmseqs2 mafft fasttree macse hyphy
```

#### Install DeepLocPro
```bash
git clone https://github.com/Jaimomar99/deeplocpro.git
cd deeplocpro && pip install .
```

#### Install Pfam Database
```bash
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

#### Install R and Visualization Packages
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "ggpattern", "glue"))
```

## Usage

### Complete Workflow Using Master Pipeline

```bash
# Complete pipeline for SARS-CoV-2
python bin/run_whole_pipeline.py -pd sars_cov_2 -n "SARS-CoV-2" \
    -tr /opt/iedb_tools --threads 8 --verbose \
    --pfam-hmm /path/to/Pfam-A.hmm

# With human negative controls instead of random
python bin/run_whole_pipeline.py -pd sars_cov_2 -n "SARS-CoV-2" \
    -tr /opt/iedb_tools --threads 8 --human-negative

# Run specific steps only
python bin/run_whole_pipeline.py -pd sars_cov_2 -n "SARS-CoV-2" \
    -tr /opt/iedb_tools --steps iedb compile random align analysis

# Skip certain steps
python bin/run_whole_pipeline.py -pd sars_cov_2 -n "SARS-CoV-2" \
    -tr /opt/iedb_tools --skip genomes pdb --threads 4

# Dry run to preview commands
python bin/run_whole_pipeline.py -pd sars_cov_2 -n "SARS-CoV-2" \
    -tr /opt/iedb_tools --dry-run
```
### Manual Step-by-Step Workflow

```bash
# 1. Fetch IEDB data
python bin/IEDB_fetch.py s_aureus "Staphylococcus aureus"

# 2. Compile antigens
python bin/compile_antigens.py s_aureus "Staphylococcus aureus"

# 3. Fetch UniProt sequences
python bin/fetch_sequences_Uniprot.py s_aureus "Staphylococcus aureus"

# 4. Generate random sequences
python bin/generate_random_sequences.py s_aureus "Staphylococcus aureus"

# 5. Fetch strain genomes
python bin/fetch_NCBI_strain_genomes.py --random "Staphylococcus aureus" \
    --random-num 6 s_aureus

# 6. Align sequences (Protein and Nucleotide)
python bin/align_antigens_mmseqs.py s_aureus "Staphylococcus aureus" \
    --threads 8 --mode protein
python bin/align_antigens_mmseqs.py s_aureus random \
    --threads 8 --mode protein --output-dir random_mmseqs_protein
python bin/align_antigens_mmseqs.py s_aureus "Staphylococcus aureus" \
    --threads 8 --mode nucleotide --output-dir mmseqs_nucleotide
python bin/align_antigens_mmseqs.py s_aureus random \
    --threads 8 --mode nucleotide --output-dir random_mmseqs_nucleotide

# 7. Fetch PDB structures
python bin/fetch_PDB_structure.py s_aureus mmseqs_protein --threads 8
python bin/fetch_PDB_structure.py s_aureus random_mmseqs_protein \
    --threads 8 --output-dir random_pdb_sequences

# 8. Fetch Pfam HMM profiles
python bin/fetch_pfam_hmmer.py s_aureus --pathogen_name "staphylococcus_aureus" \
    --pfam_hmm /path/to/Pfam-A.hmm

# 9. Run comprehensive antigen analysis
python bin/antigen_analysis.py s_aureus mmseqs_protein \
    --tool-root /path/to/tools --threads 8
python bin/antigen_analysis.py s_aureus random_mmseqs_protein \
    --tool-root /path/to/tools --threads 8 --output-dir random_analysis
python bin/antigen_analysis.py s_aureus mmseqs_nucleotide \
    --tool-root /path/to/tools --threads 8 --tools DNDS
python bin/antigen_analysis.py s_aureus random_mmseqs_nucleotide \
    --tool-root /path/to/tools --threads 8 --tools DNDS --output-dir random_analysis

# 10. Epitope prediction and evaluation
python bin/IEDB_epitope.py s_aureus mmseqs_protein \
    --tool-root /path/to/iedb_tools --threads 8
python bin/IEDB_epitope.py s_aureus random_mmseqs_protein \
    --tool-root /path/to/iedb_tools --threads 8 --output-dir random_analysis

# 11. Structure analysis
python bin/fetch_PDB_structure.py s_aureus mmseqs_protein --threads 8
python bin/fetch_PDB_structure.py s_aureus random_mmseqs_protein \
    --threads 8 --output-dir random_pdb_sequences
python bin/IEDB_epitope.py s_aureus pdb_sequences \
    --tool-root /path/to/iedb_tools --threads 8 --tools DSSP ProtLearn ElliPro
python bin/IEDB_epitope.py s_aureus random_pdb_sequences \
    --tool-root /path/to/iedb_tools --threads 8 --tools DSSP ProtLearn ElliPro \
    --output-dir random_analysis

# 12. Statistical analysis and visualization
python bin/calculate_features_kstest.py s_aureus --threads 4 --verbose --write-raw
python bin/visualise_cluster.py s_aureus epitope_outputs --split-clusters
Rscript bin/R/visualize_features.R
```

### Machine Learning Analysis

```bash
# Random Forest antigen prediction
python bin/rf_antigen_pipeline.py --base-dir results \
    --input-dir S.aureus E.coli P.aeruginosa --test-bacterium S.aureus

# PCA analysis
python bin/pca_analysis.py --base-dir results \
    --input-dir S.aureus E.coli P.aeruginosa --verbose
```

### Filtering Analysis Results
```bash
# Filter predicted antigens for high allergenicity and human homologs
python bin/filter_predicted_antigens.py --input-raw results/S_aureus_analysis_features_with_probs.csv --input-pred results/S_aureus_analysis_predictions.csv -o results/filtered_S_aureus_analysis_predictions.csv --input-fasta compiled_mmseqs_matches.fasta
```


### Key Analysis Results

The pipeline generates comprehensive results including:

- **Feature discrimination**: KS-test p-values and AUROC scores for 50+ features
- **ROC curves**: Performance visualization for each feature
- **Conservation scores**: Evolutionary analysis across multiple strains
- **Epitope predictions**: MHC-I/II binding predictions and B-cell epitopes
- **Population coverage**: HLA allele frequency analysis
- **Machine learning models**: Random Forest feature importance and PCA loadings
- **Statistical visualizations**: Distribution plots and summary statistics

### Pipeline Steps Available

The master pipeline supports these modular steps:

- `iedb`: Fetch IEDB antigen and epitope data
- `compile`: Compile antigens from multiple sources
- `uniprot`: Fetch protein sequences from UniProt
- `random`: Generate random/human control sequences
- `genomes`: Fetch strain genomes from NCBI
- `align`: Align antigens to strain genomes (protein/nucleotide)
- `pdb`: Fetch PDB structures and run structure-based predictions
- `pfam`: Fetch Pfam HMM profiles
- `analysis`: Run comprehensive feature analysis tools
- `epitopes`: Run IEDB epitope prediction tools
- `features`: Statistical feature comparison and visualization

---

## Output Files

### Data Files
- `<pathogen>_IEDB_antigens.csv`: IEDB antigen records
- `<pathogen>_compiled_proteins.csv`: Compiled protein sequences
- `random_compiled_proteins.csv`: Random control sequences
- `human_compiled_proteins.csv`: Human control sequences (if --human-negative)

### Analysis Results
- `ks_test_results_random.csv`: Statistical comparison with random controls
- `ks_test_results_human.csv`: Statistical comparison with human controls
- Feature-specific prediction files in tool subdirectories
- ROC plots and distribution visualizations
- Machine learning model outputs and feature importance

---

## Support

For questions or contributions, please refer to the project documentation in [`notes/`](notes/) or contact the development team.

---
