# General Protein Sequence Pipeline

This repository provides a flexible and scalable pipeline designed to extract, process, and retrieve protein sequences from IEDB and literature sources for various pathogens.

## Dependencies:
1. Python
   - Python 3.7 or higher
   - pandas
   - openpyxl
   - requests

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
  
## 🧭 Overview
The `get_iedb_antigens.sh` pipeline includes the following key steps:
1. **Extract UniProt IDs from Excel tables (IEDB or literature)**
2. **Fetch corresponding protein sequences in FASTA format from NCBI**
3. **Organize and clean antigen sequence data** @TODO

## 📁 Repository Structure
The repository is structured to handle multiple pathogens, with each pathogen having its own subdirectory:
- `bin/`: Contains shell scripts and Python scripts for processing.
  - `fetch_entrez_antigen_fasta.sh`: Script to fetch sequences of antigens with only names no identifiers/ ID numbers.
  - `IEDB_antigen_to_fasta.sh`: Pipeline script to batch fetch sequences of antigens with known Uniprot Identifiers
  - `parse_xlsx_txt.py`: Extracts UniProt IDs from Excel files.
  - `remove_duplicates_xlsx.py`: Removes duplicated entries in Excel files, creates queries to fetch sequences from NCBI protein database
- `data/`: Data storage section with subdirectories for each pathogen and files containing results.
  - `{Pathogen}/`: Directory specific to the target {Pathogen}, e.g., `S.aureus/`, `Influenza`.
    - Example files:
      - Input Excel files (`*.xlsx`).
      - FASTA output files (`antigens.fasta`).
      - UniProt ID list (`uniprot_ids.txt`).
- `notes/`: Reference documents and literature notes.
- `results/`: Placeholder for result files.

## Usage
There are a couple of commands that can be used:
```bash
bin/fetch_more_antigens.sh -i S{ExcelFile1}.xlsx {ExcelFile2}.xlsx... -o {OutputFile}.fasta
```
To fetch antigen sequences from IEDB or literature antigen list that contains UniProt identifiers:
```bash
bin/IEDB_antigen_to_fasta.sh [curl|wget] {input}.xlsx -o {output}.fasta
```

## Future Updates
- **Expand Pathogens**: Add support for new pathogen directories like `Influenza`, `Hepatitis-B`, etc.
- **Automation Tools**: Develop scripts to identify available Excel files and automate the pipeline execution.

## References
- `{ReferenceDocumentation}` 
- Kans J. Entrez Direct: E-utilities on the Unix Command Line. 2013 Apr 23 [Updated 2025 Mar 25]. In: Entrez Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK179288/
- Vita R, Blazeska N, Marrama D; IEDB Curation Team Members; Duesing S, Bennett J, Greenbaum J, De Almeida Mendes M, Mahita J, Wheeler DK, Cantrell JR, Overton JA, Natale DA, Sette A, Peters B. The Immune Epitope Database (IEDB): 2024 update. Nucleic Acids Res. 2025 Jan 6;53(D1):D436-D443. doi: 10.1093/nar/gkae1092. PMID: 39558162; PMCID: PMC11701597. Available from: [www.iedb.org](https://www.iedb.org/)
