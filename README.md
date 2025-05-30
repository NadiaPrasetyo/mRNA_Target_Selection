# General Protein Sequence Pipeline

This repository provides a flexible and scalable pipeline designed to extract, process, and retrieve protein sequences from IEDB and literature sources for various pathogens.

## Dependencies:
1. Python
   - Python 3.7 or higher
   - pandas==2.2.3
   - Requests==2.32.3

2. NCBI Entrez Direct
To install Entrez Direct (EDirect), open a terminal and run one of the following commands:

```sh
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
```
or
```sh
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```

This will create an `edirect` folder in your home directory and may suggest adding EDirect to your `PATH`. You can do this by running:

```sh
echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bash_profile
```

After installation, set the `PATH` for your current session:

```sh
export PATH=${HOME}/edirect:${PATH}
```

3. Seqkit
Install using Conda or Mamba
```
# conda or mamba
conda install -c bioconda seqkit
```

Usage: https://bioinf.shenwei.me/seqkit/usage/#translate

4. MMseqs2 (github: https://github.com/soedinglab/MMseqs2)
Installing options:
```sh
# install by brew
brew install mmseqs2

# install via conda
conda install -c conda-forge -c bioconda mmseqs2

# install docker
docker pull ghcr.io/soedinglab/mmseqs2

# MMseqs2-GPU mostly-static AVX2 build requiring glibc >= 2.29 and nvidia driver >=525.60.13 (see below)
wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz; tar xvfz mmseqs-linux-gpu.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

# static build with AVX2 (fastest)
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

# static build with SSE4.1
wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz; tar xvfz mmseqs-linux-sse41.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

# static build with SSE2 (slowest, for very old systems)
wget https://mmseqs.com/latest/mmseqs-linux-sse2.tar.gz; tar xvfz mmseqs-linux-sse2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

# macOS (universal, works on Apple Silicon and Intel Macs)
wget https://mmseqs.com/latest/mmseqs-osx-universal.tar.gz; tar xvfz mmseqs-osx-universal.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
```
    
## Project Structure
```
mRNA_Target_Selection/
├── bin/
│   ├── IEDB_fetch.py
│   ├── compile_antigens.py
│   ├── fetch_sequences_Uniprot.py
│   └── generate_random_sequences.py
├── data/
│   └── <pathogen_subfolders>/
│       ├── <organism_tag>_IEDB_antigens.csv
│       ├── <organism_tag>_IEDB_epitope.csv
│       ├── <organism_tag>_compiled_antigens.csv
│       ├── <organism_tag>_compiled_proteins.csv
│       └── <organism_tag>_random_proteins.csv
├── literature/
│   └── <literature_excel_files>.xlsx
├── README.md
└── requirements.txt
```

## 📁 Key Pipeline Scripts (`bin/`)

### 1. `IEDB_fetch.py`
**Description:**  
Fetches antigen and epitope data for a given organism from the IEDB API. Saves results as CSV files in the appropriate data subfolder.

**Usage:**  
```bash
python bin/IEDB_fetch.py <output_folder> <source_organism>
```
**Example:**  
```bash
python bin/IEDB_fetch.py S.aureus "Staphylococcus aureus"
```
**Output:**  
- `data/<output_folder>/<organism_tag>_IEDB_antigens.csv`
- `data/<output_folder>/<organism_tag>_IEDB_epitope.csv`

---

### 2. `compile_antigens.py`
**Description:**  
Compiles antigen data from IEDB CSV files and literature Excel files for a given organism. Produces a unified, standardized antigen list.

**Usage:**  
```bash
python bin/compile_antigens.py <short_name> <long_name>
```
**Example:**  
```bash
python bin/compile_antigens.py S.aureus "Staphylococcus aureus"
```
**Output:**  
- `data/<short_name>/<organism_tag>_compiled_antigens.csv`

---

### 3. `fetch_sequences_Uniprot.py`
**Description:**  
Fetches protein sequence and metadata from UniProt for each antigen in the compiled antigen list. Produces a CSV with detailed protein information.

**Usage:**  
```bash
python bin/fetch_sequences_Uniprot.py <pathogen> <organism>
```
**Example:**  
```bash
python bin/fetch_sequences_Uniprot.py S.aureus "Staphylococcus aureus"
```
**Output:**  
- `data/<pathogen>/<organism_tag>_compiled_proteins.csv`

---

### 4. `generate_random_sequences.py`
**Description:**  
Samples random reviewed UniProt protein entries for the organism, excluding known antigens and matching sequence length bounds. Used to generate negative/control protein sets.

**Usage:**  
```bash
python bin/generate_random_sequences.py <pathogen_subfolder> <organism_name>
```
**Example:**  
```bash
python bin/generate_random_sequences.py S.aureus "Staphylococcus aureus"
```
**Output:**  
- `data/<pathogen_subfolder>/<organism_tag>_random_proteins.csv`

---

## 🧭 Pipeline Overview

1. **Fetch IEDB antigen/epitope data:**  
   `IEDB_fetch.py`
2. **Compile antigens from IEDB and literature:**  
   `compile_antigens.py`
3. **Fetch UniProt protein sequences for antigens:**  
   `fetch_sequences_Uniprot.py`
4. **Generate random non-antigen protein set:**  
   `generate_random_sequences.py`

---


## Future Updates
- **Expand Pathogens**: Add support for new pathogen directories like `Influenza`, `Hepatitis-B`, etc.
- **Automation Tools**: Develop scripts to identify available Excel files and automate the pipeline execution.

## References
- `{ReferenceDocumentation}` 
- Kans J. Entrez Direct: E-utilities on the Unix Command Line. 2013 Apr 23 [Updated 2025 Mar 25]. In: Entrez Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK179288/
- Vita R, Blazeska N, Marrama D; IEDB Curation Team Members; Duesing S, Bennett J, Greenbaum J, De Almeida Mendes M, Mahita J, Wheeler DK, Cantrell JR, Overton JA, Natale DA, Sette A, Peters B. The Immune Epitope Database (IEDB): 2024 update. Nucleic Acids Res. 2025 Jan 6;53(D1):D436-D443. doi: 10.1093/nar/gkae1092. PMID: 39558162; PMCID: PMC11701597. Available from: [www.iedb.org](https://www.iedb.org/)
