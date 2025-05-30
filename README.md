# General Protein Sequence Pipeline

This repository provides a flexible and scalable pipeline designed to extract, process, and retrieve protein sequences from IEDB and literature sources for various pathogens.

## Dependencies:
1. Python
   - Python 3.7 or higher
   - pandas==2.2.3
   - Requests==2.32.3


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
- 
- Vita R, Blazeska N, Marrama D; IEDB Curation Team Members; Duesing S, Bennett J, Greenbaum J, De Almeida Mendes M, Mahita J, Wheeler DK, Cantrell JR, Overton JA, Natale DA, Sette A, Peters B. The Immune Epitope Database (IEDB): 2024 update. Nucleic Acids Res. 2025 Jan 6;53(D1):D436-D443. doi: 10.1093/nar/gkae1092. PMID: 39558162; PMCID: PMC11701597. Available from: [www.iedb.org](https://www.iedb.org/)
