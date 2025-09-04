### Done:
- Extended DSSP, dN/dS, and ProtLearn to the other 2 bacteria:
    - *S. aureus*
    - *S. pneumoniae*
- Compiled a list of antigens for *Chlamydia trachomatis*.
- Created the script to run the whole pipeline in steps.
- Created a script to compute and plot a PCA (still running due to large memory load).

### To Do:
1. **Trim down ProtLearn features to informative ones:**
     - ARGP820101: Hydrophobicity index (Argos et al., 1982)
     - JOND750101: Hydrophobicity (Jones, 1975)
     - BHAR880101: Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)
     - CHOC750101: Average volume of buried residue (Chothia, 1975)
     - DAYM780101: Amino acid composition (Dayhoff et al., 1978a)
     - DAYM780201: Relative mutability (Dayhoff et al., 1978b)
     - GRAR740101: Composition (Grantham, 1974)
     - GRAR740102: Polarity (Grantham, 1974)
     - GRAR740103: Volume (Grantham, 1974)
     - JOND750102: pK (-COOH) (Jones, 1975)
     - KYTJ820101: Hydropathy index (Kyte-Doolittle, 1982)

2. **Pull human protein negative sets** to see if it would change the performance of epitope prediction tools.

3. **Take statistical summaries** of each feature for each antigen:
     - Ensure memory load remains low.

4. **Adding glycoprotein antigen epitope prediction tool**: [SEPPA3.0](https://academic.oup.com/nar/article/47/W1/W388/5494762)
