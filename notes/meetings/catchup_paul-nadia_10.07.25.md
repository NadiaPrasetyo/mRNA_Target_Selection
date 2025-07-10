**Completed Tasks:**
- Interpreted cluster results
- Analyzed AlgPred results
- Reviewed PopCov results
- Integrated all above features into the KS test
- Generated visualization data to compare feature differences between random and positive sets
- Developed a script to extract datasets and summarize them into: sequence, label, feature, subfeature, value

**To Do:**
- Incorporate antigenicity and IFN-gamma prediction into the pipeline (commonly included in vaccine target selection pipelines)
- Investigate patent data for positive evaluation (in progress)
- Update conservation scores to utilize both sum of percent identities and log₁₀(max e-value) - dN/dS to check specific clusters
    - Note: Bit score depends on evolutionary rate and sequence length; longer, slower-evolving sequences yield higher bit scores
- Adjust cluster sensitivity parameters
- Explore new tools for extracellular prediction, e.g., DeepLocPro for subcellular localization of prokaryotic proteins
- Evaluate additional relevant tools