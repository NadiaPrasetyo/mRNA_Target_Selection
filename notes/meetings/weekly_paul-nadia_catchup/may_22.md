### Candidate Vaccine for *S. aureus*

- Protein-based construct: multiple *S. aureus* proteins linked together into one protein.

#### Clonal Complexes

- Reese Langley: Focus on 6 circulating strains of *S. aureus* that explain 95% of infections globally.
- Covers all of these; focus on these 6 strains:
    - CC1
    - CC5
    - CC8
    - CC22
    - CC30
    - CC93

> **Note:**  
> The 6 Clonal Complexes that cover the majority of *Staphylococcus aureus* infections are CC1, CC5, CC8, CC22, CC30, and CC93.  
>  
> I have attached a Word doc with a table providing information on the genomes I got my SSL sequences from.  
>  
> I have also included the nucleotide and aa sequences of all SSLs (SSL3, SSL7, and SSL11 variants in our vaccine combinations).

I’ve also included some papers on ssl (old name set) evolution.

---

### Preferred Protein Databases

- UniProt
- Pfam (InterPro)

---

### How to Get Sequence Motifs of Antigens

Automate looking for candidate antigens of different pathogens:

- Paul knows absolutely nothing about: find some machine learning tools.
- Look for in previous publications and literature based on experimental validations.
    - Most people use complete proteins (not chopped up into domains): least well-conserved parts of the protein (C and N ends of proteins).

---

### List of Ideas

- If we have poor literature: what are the features that make something a reasonable antigen candidate?
    - Exported: signal peptides (annotated from SignalP)
    - Transmembrane proteins: predict if they are internal and external
    - External proteins
    - Cell-wall anchored proteins
    - Toxins: consider how rapidly the toxin changes (evolution)
    - Things under positive selection (e.g., spike protein is one of the fastest changing molecules: 100+ non-synonymous changes and only 4 synonymous changes compared to Wuhan)
        - Change the vaccine often enough (what we do now with flu and COVID)
        - Cheap for mRNA but expensive testing
- How would I find the different features from each pathogen?
- Could use evolutionary information:
    - jackhammer (self-made alignments, multiple iterations)
    - eggnog (protein full-length alignment database)  
      *Based on JUST the 6 genomes*

#### Pipeline Outline

1. Get the epitopes from IEDB.
2. Map to the genome of the pathogen—find out the motifs and features why these are good epitopes.
3. Find more potential antigens from that information.
4. Get the epitopes/predict the epitopes from each antigen from pathogens.
5. Take the epitopes and see if they conform to any of the features.
6. Use data for training model.

- Use existing patents for training data.
- Vaccine databases for training data.
- Epitope databases for training data.

**=> Train own model**