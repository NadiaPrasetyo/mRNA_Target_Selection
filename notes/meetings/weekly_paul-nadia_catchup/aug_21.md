## Done

- Implemented fetching nucleotide sequence and Pfam ID fetching from `fetch_sequence_Uniprot`.
- Implemented fetching the HMM files of Pfams associated with all the antigens and proteins.
- Ran *S. pyogenes* analysis.
- Implemented DNDS script using [HyPhy](https://hyphy.org/hbl-reference/).

## To Do

1. Fetch all the domains of interest from Pfam and create an HMM concatenation of all; run `hmmscan` on the protein sequences.
2. Analyze data from AlphaFold PDB files using [ProtLearn](https://protlearn.readthedocs.io/en/latest/dimensionality_reduction.html).
3. Run HyPhy DNDS statistical analysis.
4. Perform statistical tests:
    - Run t-test to understand directional differences.
    - Run KS test and AUROC.
5. Implement and run [DSSP](https://pdb-redo.eu/dssp/about).
6. Extend analysis to other bacteria (Gram-negative).
7. Create a summary of NetMHCpan methods to understand why it is failing.

---