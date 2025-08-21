

### NetMHCpan (v4.1) Overview

NetMHCpan‑4.1 predicts peptide binding to MHC‑I molecules by training a feed‑forward neural network on large, mixed datasets of binding‑affinity (BA) measurements and mass‑spectrometry eluted‑ligand (EL) observations (both single‑allele SA and multi‑allele MA). The NNAlign‑MA framework first learns a pan‑specific model on SA data (burn‑in), then iteratively annotates MA peptides to the most likely allele and retrains, allowing simultaneous motif deconvolution and prediction.

#### Key Questions

1. **How do they split their training data into positive and negative sets?**
    - Positive EL peptides are experimentally observed ligands.
    - Negatives are random natural peptides of matching length, added at a 5‑fold excess for EL (SA) and at a 5‑fold excess of the most represented positive length for EL (MA).
    - For BA, peptides with IC₅₀ ≤ 500 nM are labeled positive; all others are negative.

2. **What are considered positive and negative sets?**
    - **Positive:** Binding‑affinity hits (≤ 500 nM) or EL ligands.
    - **Negative:** Synthetic decoys drawn from the UniProt proteome (random peptides).

3. **How is the neural network defining the input and creating the outputs?**
    - **Input:** One‑hot encoded peptide sequence (length‑binned for class I) and MHC pseudo‑sequence (pan‑specific encoding).
    - **Output:** Two neurons:
      - Predicts BA (continuous affinity).
      - Predicts EL‑likelihood (probability of being a presented ligand).

4. **What features/fields are used to determine the labels and outputs?**
    - **Features:** Peptide amino‑acid identity, peptide length bin, and allele‑specific pseudo‑sequence residues.
    - **Labels:** Quantitative BA values or binary EL presence/absence.

5. **What do the outputs mean, and how are they relevant biologically?**
    - **BA Output:** Converted to predicted affinity (nM) and %Rank.
    - **EL Output:** Likelihood score and %Rank EL, indicating how strongly a peptide is expected to be naturally presented.
    - These scores guide epitope selection for vaccine or immunotherapy design.

6. **How are the outputs and the machine learning model evaluated?**
    - **Evaluation Metrics:**
      - 5‑fold cross‑validation using AUC, AUC₀.₁ (up to 10 % FPR), Pearson correlation (PCC), and Positive Predictive Value (PPV).
      - Independent epitope benchmarks report median FRANK values (e.g., 0.00220 for CD8⁺ epitopes) and PPV (≈0.82 for HLA‑DR).

---

### Positive and Negative Training Sets in NetMHCpan‑4.1

| **Aspect**               | **What Was Done**                                                                                     | **Why It Makes Sense Biologically**                                                                                                                                                                                                 |
|--------------------------|-------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Positive BA Instances** | Peptides with experimentally measured binding affinity ≤ 500 nM were labeled positive.               | A ≤ 500 nM IC₅₀ is widely accepted as “strong binder” for MHC‑I; such peptides are truly capable of forming stable peptide‑MHC complexes in vivo.                                                                                 |
| **Positive EL Instances** | All peptides identified by mass‑spectrometry as eluted ligands (SA and MA) were labeled positive.    | EL peptides have survived the entire antigen‑processing pathway and were physically captured from the cell surface, representing genuine presented epitopes.                                                                      |
| **Negative BA Instances** | 100 random peptides were drawn from the UniProt proteome and assigned a target value of 0.01.        | Random natural peptides are unlikely to bind strongly; adding a small non‑zero value prevents division‑by‑zero errors while preserving the “non‑binder” nature of the examples.                                                  |
| **Negative EL Instances** | Random peptides of the same length distribution were sampled from UniProt and added at 5× the count of the most abundant positive length. | Oversampling negatives creates a realistic background that mimics the pool of intracellular peptides that never reach the cell surface, improving the network’s ability to discriminate true ligands from noise.                  |
| **Source Proteomes**      | Random sequences were extracted from the UniProt database.                                           | UniProt contains proteins from all organisms, providing a diverse, unbiased representation of the peptide space that a cell could theoretically generate. Using this broad background avoids systematic bias toward any species. |
| **Length Filtering**      | Training peptides were limited to 8–14 aa for class I and 13–21 aa for class II.                    | These ranges correspond to the natural length preferences of MHC‑I and MHC‑II binding grooves, ensuring that both positives and negatives are biologically plausible candidates for presentation.                                  |

#### Why the Sets Are Appropriate

- Positives come from direct binding assays (BA) or actual presentation (EL), ensuring the model learns true binding motifs.
- Negatives are synthetic decoys drawn from the same organismal proteome pool, preserving amino‑acid composition and length distribution while representing peptides that never bind or are never presented.
- The 5‑fold enrichment of negatives for EL mirrors the real cellular environment where most peptides are not presented, helping the network learn realistic decision boundaries.

This combination of experimentally verified binders and carefully constructed decoys provides a balanced yet biologically realistic training landscape, allowing the NNAlign‑MA architecture to capture genuine peptide‑MHC interaction signals while minimizing overfitting to artifacts.



NetMHCpan‑4.1 does not predict whether a whole protein is an “antigen” or a “non‑antigen”. It predicts, for every 9‑mer (or 8‑14 aa) peptide derived from a sequence, the likelihood that the peptide will bind a given HLA‑I molecule (EL‑likelihood) and, optionally, a quantitative affinity (BA) value. The training data used by the server consist of individual peptides that were either (i) measured as binders in vitro (BA ≤ 500 nM) or (ii) observed as naturally eluted ligands by mass‑spectrometry (EL = 1). Negatives are random natural peptides drawn from the UniProt proteome and added at a 5‑fold excess for EL data (or 100 random peptides per allele for BA data)  
. Because the model learns only the sequence motifs that favour binding, any peptide that happens to contain those motifs will receive a high EL score, regardless of the biological context of the source protein.
When you feed full‑length bacterial proteins (positive set) and length‑matched random proteins (negative set) to NetMHCpan, the server first slides a window across each protein and scores all overlapping 9‑mers. Random proteins contain a large number of peptides that, by chance, match the learned HLA motifs, so their aggregate EL‑score distribution is almost identical to that of bona‑fide antigens. Moreover, the model was trained almost exclusively on human‑derived peptides (human BA and EL data 
‑45); bacterial peptide repertoires are under‑represented, so the network cannot distinguish bacterial “true” epitopes from background peptides. Finally, antigenicity depends on many steps that are outside the scope of binding prediction (protein expression level, proteasomal cleavage, TAP transport, peptide stability, T‑cell receptor recognition). NetMHCpan therefore cannot be expected to separate “antigen” versus “non‑antigen” protein sets; it only separates binding versus non‑binding peptides.
To observe discrimination you would need a peptide‑level benchmark (e.g., known bacterial T‑cell epitopes vs. random 9‑mers) and, if desired, combine NetMHCpan scores with additional processing predictors (e.g., proteasomal cleavage, TAP transport) before aggregating to the protein level. This explains the near‑random AUROC (≈0.5) you obtained.

**TLDR**

**NetMHCpan (v4.1) Overview:**
- Predicts peptide binding to MHC‑I molecules using a feed‑forward neural network trained on binding-affinity (BA) and mass-spectrometry eluted-ligand (EL) data.
- **Key Questions:**
  - **How do they split training data?**
     - Positives: EL peptides (observed ligands) or BA peptides (IC₅₀ ≤ 500 nM).
     - Negatives: Random peptides from UniProt, added at 5× excess.
  - **Input/Output:**
     - Input: One-hot encoded peptide sequence and MHC pseudo-sequence.
     - Output: BA (affinity) and EL-likelihood (ligand probability).
  - **Evaluation Metrics:**
     - AUC, AUC₀.₁, Pearson correlation, PPV, and independent benchmarks.

---
**Positive/Negative Training Sets:**

- **BA Instances:**
    - **Positive:** Peptides with IC₅₀ ≤ 500 nM.
    - **Negative:** 100 random peptides from UniProt, target value = 0.01.

- **EL Instances:**
    - **Positive:** Mass-spectrometry eluted ligands.
    - **Negative:** Random peptides (same length distribution) at 5× excess.

- **Source Proteomes:**
    - **Positive:** Experimentally verified peptides.
    - **Negative:** Random sequences from UniProt.

- **Length Filtering:**
    - **Positive:** 8–14 aa (MHC‑I) and 13–21 aa (MHC‑II).
    - **Negative:** Same length ranges for negatives.

---

**Why NetMHCpan Struggles with Antigenicity:**
- Trained on human-derived peptides; bacterial peptides are under-represented.
- Predicts binding/non-binding peptides, not antigenicity.
- Random proteins often match HLA motifs by chance, leading to near-random AUROC (≈0.5).
- Antigenicity depends on additional factors (e.g., proteasomal cleavage, TAP transport, T-cell recognition).

--- 