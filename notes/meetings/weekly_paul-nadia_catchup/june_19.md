## To Report

- Completed epitope evaluation (script currently running)
- Retrieved feature data for randoms
- Obtained KS test results comparing features between random and positive sets for *S. aureus*

## Next Steps

- Add new literature to the reading list and review relevant papers
- Try out MARIA for both MHCI and MHCII analysis (as an alternative to NetMHCpan)
- Replace current clustering approach with MMseqs2 for epitope analysis
    - Use MMseqs2 to obtain similarity scores (max similarity between query and each database sequence)
    - Calculate percentage identity
- Use KS tests to compare performance of different analysis pipelines (e.g., MMseqs2 vs. previous clustering)
- Evaluate and choose between SignalP and TargetP (TMHMM is redundant and will be scrapped)
- Work on completing datasets, including adding epitope patents for accuracy evaluation
