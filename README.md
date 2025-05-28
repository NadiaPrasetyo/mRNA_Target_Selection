# General Protein Sequence Pipeline

This repository provides a flexible and scalable pipeline designed to extract, process, and retrieve protein sequences from IEDB and literature sources for various pathogens.

## Dependencies:
1. Python
   - Python 3.7 or higher
   - pandas
   - openpyxl
   - requests
   - biopython

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


  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:09:28 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 1\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 1"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:09:31 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 1\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_6836624c0020c933cf0cf135</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus HO 5096 0412"[All Fields] AND "staphylococcal superantigen-like protein 1"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 1"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 1
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 3
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 5
 WARNING:  FAILURE ( Wed May 28 01:09:37 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 6\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 6"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:09:39 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 6\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 6"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:09:42 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 6\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_68366257622d6aeced02a6c6</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus HO 5096 0412"[All Fields] AND "staphylococcal superantigen-like protein 6"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 6"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 6
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 7
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 10
 WARNING:  FAILURE ( Wed May 28 01:09:48 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 11\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 11"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:09:51 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 11\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 11"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:09:53 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 11\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_68366262f3b71a8bf102015b</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus HO 5096 0412"[All Fields] AND "staphylococcal superantigen-like protein 11"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 11"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 11
 WARNING:  FAILURE ( Wed May 28 01:09:57 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 13\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 13"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:09:59 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 13\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 13"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:10:01 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphylococcal superantigen-like protein 13\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_6836626a5cef84b1a90f0479</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus HO 5096 0412"[All Fields] AND "staphylococcal superantigen-like protein 13"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 13"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 13
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal peroxidase inhibitor
 WARNING:  FAILURE ( Wed May 28 01:10:25 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphopain cysteine protease\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphopain cysteine protease"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:10:28 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphopain cysteine protease\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphopain cysteine protease"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:10:30 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus HO 5096 0412\"[All Fields] AND \"staphopain cysteine protease\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_68366287d1214b1dc3025a59</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus HO 5096 0412"[All Fields] AND "staphopain cysteine protease"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphopain cysteine protease"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphopain cysteine protease
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: toxic shock syndrome toxin-1
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: von willebrand factor-binding protein
  [INFO] Found 0 protein entries
  [INFO] Matched 259 proteins so far
[INFO] Processing strain MRSA252 (EMBL: BX571856) [5]
  [INFO] Found 888 protein entries
  [INFO] Matched 339 proteins so far
[INFO] Processing strain JKD6159 (EMBL: CP002114) [6]
[WARN] No UniProt entries found for EMBL ID CP002114, trying NCBI...
[WARN] Empty GenPept result for strain: JKD6159, antigen: adenosine synthase a
[WARN] Empty GenPept result for strain: JKD6159, antigen: clumping factor a
[WARN] Empty GenPept result for strain: JKD6159, antigen: collagen adhesin
[WARN] Empty GenPept result for strain: JKD6159, antigen: coagulase
 WARNING:  FAILURE ( Wed May 28 01:11:19 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"extracellular adherence proteins\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"extracellular adherence proteins"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:11:21 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"extracellular adherence proteins\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"extracellular adherence proteins"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:11:23 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"extracellular adherence proteins\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_683662bc8fb3876d37092b1c</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "extracellular adherence proteins"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"extracellular adherence proteins"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: extracellular adherence proteins
[WARN] Empty GenPept result for strain: JKD6159, antigen: enterotoxins b
[WARN] Empty GenPept result for strain: JKD6159, antigen: flipr
[WARN] Empty GenPept result for strain: JKD6159, antigen: fibronectin-binding protein a
[WARN] Empty GenPept result for strain: JKD6159, antigen: fibronectin-binding protein b
 WARNING:  FAILURE ( Wed May 28 01:11:37 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"haemolysin ab and cb\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"haemolysin ab and cb"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:11:40 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"haemolysin ab and cb\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"haemolysin ab and cb"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:11:42 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"haemolysin ab and cb\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_683662cf65991160b8026345</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "haemolysin ab and cb"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"haemolysin ab and cb"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: haemolysin ab and cb
[WARN] Empty GenPept result for strain: JKD6159, antigen: leukocidin ab
[WARN] Empty GenPept result for strain: JKD6159, antigen: leukocidin ed
 WARNING:  FAILURE ( Wed May 28 01:11:48 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"leukocidin sf-pv\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"leukocidin sf-pv"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:11:51 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"leukocidin sf-pv\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"leukocidin sf-pv"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:11:53 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"leukocidin sf-pv\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_683662da258d99af5301705c</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "leukocidin sf-pv"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"leukocidin sf-pv"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: leukocidin sf-pv
 WARNING:  FAILURE ( Wed May 28 01:12:19 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal binder of immunoglobulin\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal binder of immunoglobulin"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:12:22 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal binder of immunoglobulin\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal binder of immunoglobulin"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:12:24 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal binder of immunoglobulin\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_683662f9d02f31b42c0438db</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "staphylococcal binder of immunoglobulin"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal binder of immunoglobulin"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal binder of immunoglobulin
[WARN] Empty GenPept result for strain: JKD6159, antigen: serine-aspartate repeat protein d
[WARN] Empty GenPept result for strain: JKD6159, antigen: serine-aspartate repeat protein e
 WARNING:  FAILURE ( Wed May 28 01:12:30 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"seiw\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"seiw"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:12:33 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"seiw\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"seiw"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:12:35 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"seiw\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_683663044f596be59e0c5d84</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "seiw"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"seiw"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: seiw
 WARNING:  FAILURE ( Wed May 28 01:12:43 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 1\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 1"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:12:45 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 1\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 1"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:12:48 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 1\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_68366311711bf41c0f0ec424</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "staphylococcal superantigen-like protein 1"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 1"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 1
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 3
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 5
 WARNING:  FAILURE ( Wed May 28 01:12:54 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 6\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 6"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:12:56 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 6\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 6"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:12:59 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 6\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_6836631ca43c202a350e5f75</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "staphylococcal superantigen-like protein 6"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 6"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 6
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 7
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 10
 WARNING:  FAILURE ( Wed May 28 01:13:05 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 11\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 11"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:13:07 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 11\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 11"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:13:10 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 11\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_68366327014270edd707ca1b</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "staphylococcal superantigen-like protein 11"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 11"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 11
 WARNING:  FAILURE ( Wed May 28 01:13:13 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 13\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 13"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:13:16 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 13\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 13"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:13:18 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphylococcal superantigen-like protein 13\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_6836632f3c1aec03e7053f48</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "staphylococcal superantigen-like protein 13"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphylococcal superantigen-like protein 13"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 13
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal peroxidase inhibitor
 WARNING:  FAILURE ( Wed May 28 01:13:37 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphopain cysteine protease\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphopain cysteine protease"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
SECOND ATTEMPT
 WARNING:  FAILURE ( Wed May 28 01:13:39 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphopain cysteine protease\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<WarningList>
  <QuotedPhraseNotFound>"staphopain cysteine protease"[All Fields]</QuotedPhraseNotFound>
  <OutputMessage>No items found.</OutputMessage>
<WarningList>
LAST ATTEMPT
 ERROR:  FAILURE ( Wed May 28 01:13:42 PM NZST 2025 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ esearch.fcgi -retmax 0 -usehistory y -db protein -term "\"Staphylococcus aureus subsp. aureus JKD6159\"[All Fields] AND \"staphopain cysteine protease\"[All Fields]" -tool edirect -edirect 24.0 -edirect_os Linux -email prana47p@aoraki07
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult>
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <QueryKey>1</QueryKey>
  <WebEnv>MCID_683663470d0a4658a30dbc1c</WebEnv>
  <QueryTranslation>"Staphylococcus aureus subsp. aureus JKD6159"[All Fields] AND "staphopain cysteine protease"[All Fields]</QueryTranslation>
  <WarningList>
    <QuotedPhraseNotFound>"staphopain cysteine protease"[All Fields]</QuotedPhraseNotFound>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>
QUERY FAILURE
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphopain cysteine protease
[WARN] Empty GenPept result for strain: JKD6159, antigen: toxic shock syndrome toxin-1
  [INFO] Found 0 protein entries
  [INFO] Matched 375 proteins so far
[DONE] Wrote 375 proteins to: data/S.aureus/S.aureus_compiled_proteins.csv
(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ ^C
(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ ^C
(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ git pull
remote: Enumerating objects: 7, done.
remote: Counting objects: 100% (7/7), done.
remote: Compressing objects: 100% (1/1), done.
remote: Total 4 (delta 3), reused 4 (delta 3), pack-reused 0 (from 0)
Unpacking objects: 100% (4/4), 536 bytes | 15.00 KiB/s, done.
From https://github.com/NadiaPrasetyo/mRNA_Target_Selection
   5820154..3bb2a18  main       -> origin/main
Updating 5820154..3bb2a18
Fast-forward
 bin/fetch_sequences_Uniprot.py | 10 ++++++++--
 1 file changed, 8 insertions(+), 2 deletions(-)
(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ python bin/fetch_sequences_Uniprot.py S.aureus
[INFO] Loaded 92 antigen keyword patterns (UniProt)
[INFO] Loaded 39 cleaned antigen names (NCBI)
[INFO] Processing strain MSSA476 (EMBL: BX571857) [1]
  [INFO] Found 890 protein entries
  [INFO] Matched 80 proteins so far
[INFO] Processing strain N315 (EMBL: BA000018) [2]
  [INFO] Found 928 protein entries
  [INFO] Matched 162 proteins so far
[INFO] Processing strain NCTC 8325 (EMBL: CP000253) [3]
  [INFO] Found 810 protein entries
  [INFO] Matched 238 proteins so far
[INFO] Processing strain HO 5096 0412 (EMBL: HE681097) [4]
[WARN] No UniProt entries found for EMBL ID HE681097, trying NCBI...
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: adenosine synthase a
^CTraceback (most recent call last):
  File "/projects/health_sciences/bms/biochemistry/gardner_group/nadia/mRNA_Target_Selection/bin/fetch_sequences_Uniprot.py", line 307, in <module>
    main(sys.argv[1])
    ~~~~^^^^^^^^^^^^^
  File "/projects/health_sciences/bms/biochemistry/gardner_group/nadia/mRNA_Target_Selection/bin/fetch_sequences_Uniprot.py", line 274, in main
    gp_text = fetch_protein_data_ncbi(strain, antigen)
  File "/projects/health_sciences/bms/biochemistry/gardner_group/nadia/mRNA_Target_Selection/bin/fetch_sequences_Uniprot.py", line 51, in fetch_protein_data_ncbi
    output = subprocess.check_output(cmd, shell=True, text=True)
  File "/home/prana47p/miniconda3/lib/python3.13/subprocess.py", line 474, in check_output
    return run(*popenargs, stdout=PIPE, timeout=timeout, check=True,
           ~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
               **kwargs).stdout
               ^^^^^^^^^
  File "/home/prana47p/miniconda3/lib/python3.13/subprocess.py", line 558, in run
    stdout, stderr = process.communicate(input, timeout=timeout)
                     ~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/prana47p/miniconda3/lib/python3.13/subprocess.py", line 1208, in communicate
    stdout = self.stdout.read()
KeyboardInterrupt

(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ python bin/fetch_sequences_Uniprot.py S.aureus
[INFO] Loaded 92 antigen keyword patterns (UniProt)
[INFO] Loaded 39 cleaned antigen names (NCBI)
[INFO] Processing strain MSSA476 (EMBL: BX571857) [1]
  [INFO] Found 890 protein entries
  [INFO] Matched 80 proteins so far
[INFO] Processing strain N315 (EMBL: BA000018) [2]
  [INFO] Found 928 protein entries
  [INFO] Matched 162 proteins so far
[INFO] Processing strain NCTC 8325 (EMBL: CP000253) [3]
  [INFO] Found 810 protein entries
  [INFO] Matched 238 proteins so far
[INFO] Processing strain HO 5096 0412 (EMBL: HE681097) [4]
[WARN] No UniProt entries found for EMBL ID HE681097, trying NCBI...
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: adenosine synthase a
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: chemotaxi inhibitory protein
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: clumping factor a
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: coagulase
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: extracellular adherence protein
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: ecb
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: flipr
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: fibronectin-binding protein a
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: haemolysin ab and cb
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: leukocidin ab
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: leukocidin ed
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: leukocidin sf-pv
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: protein a
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal binder of immunoglobulin
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: serine-aspartate repeat protein e
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: seiw
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal enterotoxin-like toxin x
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 1
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 3
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 5
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 6
^CTraceback (most recent call last):
  File "/projects/health_sciences/bms/biochemistry/gardner_group/nadia/mRNA_Target_Selection/bin/fetch_sequences_Uniprot.py", line 307, in <module>
    main(sys.argv[1])
    ~~~~^^^^^^^^^^^^^
  File "/projects/health_sciences/bms/biochemistry/gardner_group/nadia/mRNA_Target_Selection/bin/fetch_sequences_Uniprot.py", line 274, in main
    gp_text = fetch_protein_data_ncbi(strain, antigen)
  File "/projects/health_sciences/bms/biochemistry/gardner_group/nadia/mRNA_Target_Selection/bin/fetch_sequences_Uniprot.py", line 51, in fetch_protein_data_ncbi
    output = subprocess.check_output(cmd, shell=True, text=True)
  File "/home/prana47p/miniconda3/lib/python3.13/subprocess.py", line 474, in check_output
    return run(*popenargs, stdout=PIPE, timeout=timeout, check=True,
           ~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
               **kwargs).stdout
               ^^^^^^^^^
  File "/home/prana47p/miniconda3/lib/python3.13/subprocess.py", line 558, in run
    stdout, stderr = process.communicate(input, timeout=timeout)
                     ~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/prana47p/miniconda3/lib/python3.13/subprocess.py", line 1208, in communicate
    stdout = self.stdout.read()
KeyboardInterrupt

(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ git pull
remote: Enumerating objects: 7, done.
remote: Counting objects: 100% (7/7), done.
remote: Compressing objects: 100% (1/1), done.
remote: Total 4 (delta 3), reused 4 (delta 3), pack-reused 0 (from 0)
Unpacking objects: 100% (4/4), 343 bytes | 9.00 KiB/s, done.
From https://github.com/NadiaPrasetyo/mRNA_Target_Selection
   3bb2a18..94cf4ff  main       -> origin/main
Updating 3bb2a18..94cf4ff
Fast-forward
 bin/fetch_sequences_Uniprot.py | 3 ---
 1 file changed, 3 deletions(-)
(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ python bin/fetch_sequences_Uniprot.py S.aureus
[INFO] Loaded 92 antigen keyword patterns (UniProt)
[INFO] Loaded 39 cleaned antigen names (NCBI)
[INFO] Processing strain MSSA476 (EMBL: BX571857) [1]
  [INFO] Found 890 protein entries
  [INFO] Matched 80 proteins so far
[INFO] Processing strain N315 (EMBL: BA000018) [2]
  [INFO] Found 928 protein entries
  [INFO] Matched 162 proteins so far
[INFO] Processing strain NCTC 8325 (EMBL: CP000253) [3]
  [INFO] Found 810 protein entries
  [INFO] Matched 238 proteins so far
[INFO] Processing strain HO 5096 0412 (EMBL: HE681097) [4]
[WARN] No UniProt entries found for EMBL ID HE681097, trying NCBI...
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: adenosine synthase a
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: clumping factor a
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: coagulase
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: extracellular adherence proteins
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: enterotoxins b
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: ecb
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: flipr
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: fibronectin-binding protein a
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: haemolysin ab and cb
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: leukocidin ab
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: leukocidin ed
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: leukocidin sf-pv
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: protein a
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal binder of immunoglobulin
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: serine-aspartate repeat protein e
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: seiw
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal enterotoxin-like toxin x
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 1
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 3
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 5
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 6
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 7
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 10
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 11
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal superantigen-like protein 13
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphylococcal peroxidase inhibitor
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: staphopain cysteine protease
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: toxic shock syndrome toxin-1
[WARN] Empty GenPept result for strain: HO 5096 0412, antigen: von willebrand factor-binding protein
  [INFO] Found 0 protein entries
  [INFO] Matched 253 proteins so far
[INFO] Processing strain MRSA252 (EMBL: BX571856) [5]
  [INFO] Found 888 protein entries
  [INFO] Matched 333 proteins so far
[INFO] Processing strain JKD6159 (EMBL: CP002114) [6]
[WARN] No UniProt entries found for EMBL ID CP002114, trying NCBI...
[WARN] Empty GenPept result for strain: JKD6159, antigen: adenosine synthase a
[WARN] Empty GenPept result for strain: JKD6159, antigen: clumping factor a
[WARN] Empty GenPept result for strain: JKD6159, antigen: collagen adhesin
[WARN] Empty GenPept result for strain: JKD6159, antigen: coagulase
[WARN] Empty GenPept result for strain: JKD6159, antigen: extracellular adherence proteins
[WARN] Empty GenPept result for strain: JKD6159, antigen: enterotoxins b
[WARN] Empty GenPept result for strain: JKD6159, antigen: flipr
[WARN] Empty GenPept result for strain: JKD6159, antigen: fibronectin-binding protein a
[WARN] Empty GenPept result for strain: JKD6159, antigen: haemolysin ab and cb
[WARN] Empty GenPept result for strain: JKD6159, antigen: leukocidin ab
[WARN] Empty GenPept result for strain: JKD6159, antigen: leukocidin ed
[WARN] Empty GenPept result for strain: JKD6159, antigen: leukocidin sf-pv
[WARN] Empty GenPept result for strain: JKD6159, antigen: protein a
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal binder of immunoglobulin
[WARN] Empty GenPept result for strain: JKD6159, antigen: serine-aspartate repeat protein e
[WARN] Empty GenPept result for strain: JKD6159, antigen: seiw
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 1
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 3
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 5
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 6
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 7
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 10
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 11
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal superantigen-like protein 13
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphylococcal peroxidase inhibitor
[WARN] Empty GenPept result for strain: JKD6159, antigen: staphopain cysteine protease
[WARN] Empty GenPept result for strain: JKD6159, antigen: toxic shock syndrome toxin-1
  [INFO] Found 0 protein entries
  [INFO] Matched 363 proteins so far
[DONE] Wrote 363 proteins to: data/S.aureus/S.aureus_compiled_proteins.csv
(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ git status
On branch main
Your branch is up to date with 'origin/main'.

Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git restore <file>..." to discard changes in working directory)
	modified:   data/S.aureus/S.aureus_compiled_proteins.csv

no changes added to commit (use "git add" and/or "git commit -a")
(base) [prana47p@rtis-hpc-r07 mRNA_Target_Selection]$ 
