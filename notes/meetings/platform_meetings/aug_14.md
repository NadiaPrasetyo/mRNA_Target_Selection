## Fold AI Progress

### Indication Agnostic Analysis
- Using decision-making trees but hit a wall, waiting on an indication.

### Methods
1. **Manual Approach** (Holly's approach):
  - Standard prioritization based on conservation and regions of interest.
2. **AI Approach** (Ryosuke's approach):
  - Developing a machine learning pipeline to predict gene targets.
  - Focused on the proteome of interest: most probable antigens.
  - Working on an extended rendition to identify epitopes.
  - **Benchmarking**:
    - Positive and negative datasets curated from ImmunoDB proteins.
    - Negative dataset includes non-antigen/non-epitope annotated sequences.
  - **Improvements Needed**:
    - Refine T-cell and B-cell epitope separation.
    - Identify protein domains with both T-cell and B-cell epitopes.
    - Potential to create novel proteins with correct T-cell and B-cell epitopes.

### Desired AI Model Output
- List of epitopes.
- Recombinant synthetic/novel proteins including the epitopes.
- HLA allele binding information.
- Population coverage analysis.

---

## Lisa's Progress

### Influenza HA Stalk Modifications
- Stripped the influenza HA to just the stalk.
- Modified the stalk using AlphaFold to improve stability:
  - Achieved a 10x fold response improvement of the modified stalk vs. wildtype stalk.
- Attached a single epitope to the modified stalk:
  - Achieved a 10x fold response improvement of the modified stalk + epitope vs. just the modified stalk.

### Current Work
- Exploring other approaches:
  - Adding multiple epitopes or full domains (limited by protein size constraints).
- Creating head-only constructs (unknown head).