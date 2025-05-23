# General Protein Sequence Pipeline

This repository provides a flexible and scalable pipeline designed to extract, process, and retrieve protein sequences from IEDB and literature sources for various pathogens.

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
bin/fetch_more_antigens.sh {ExcelFile1}.xlsx {ExcelFile2}.xlsx... -o {OutputDirectory} @TODO
```
To fetch antigen sequences from IEDB or literature antigen list that contains UniProt identifiers:
```bash
bin/IEDB_antigen_to_fasta.sh [curl|wget] {input}.xlsx -o {output}.fasta
```

## Future Updates
- **Expand Pathogens**: Add support for new pathogen directories like `Influenza`, `Hepatitis-B`, etc.
- **Automation Tools**: Develop scripts to identify available Excel files and automate the pipeline execution.

## References
- `{ReferenceDocumentation}` @Todo add later NCBI, IEDB, 