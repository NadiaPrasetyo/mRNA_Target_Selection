## To Report

- **Epitope analysis scripts completed:**  
    - Results obtained for MHC-I, MHC-II, and B-cell epitopes  
    - Includes binding scores and percentiles  
    - Used 27 MHC-I alleles, 7 main MHC-II alleles, and various prediction methods

- **Antigen analysis scripts completed:**  
    - Results from SignalP, TargetP, and TMHMM  
    - Includes locations, association probability scores, types, and plots

---

## Next Steps: get more features

- **Evaluate epitope quality:**  
    - Allergenicity  
    - Antigenicity  
    - Population coverage
    - Conservancy

- **Simulate protein-protein interaction:**  
    - Molecular docking with MHC-I

---

## Data Visualization

- Fill in feature scores maps for each sequence to understand which is best as a filter  
- Identify best features for predicting antigens and epitopes from pathogen genome

---

## Controls

- **Positive controls:**  
    - Known interactors with MHC-I, MHC-II, B-cell

- **Negative controls:**  
    - Random sequences

---

## Analysis

- Fill in the feature matrix for each sequence
- Compare negative and positive controls:
    - Perform KS test for each feature to assess informativeness

- Develop a pipeline to analyze conservation across strains:
    - Highly conserved: likely highly functional
    - Not conserved: most positively conserved, likely to interact with immune system