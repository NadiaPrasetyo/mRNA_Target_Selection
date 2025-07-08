Existing mRNA Vaccines used for:
- COVID-19: BNT162b2, 
- Zika virus
- Influenza

## Pipeline steps for epitope selection, prediction, and evaluation:

1. Choose pathogen of interest:
Universal Flu (influnenza),Staphylococcus aureus,Lupus (mCAR),Chronic HepB

2. Survey existing literature to understand mechanism of actiona and attempted treatments

3. Select target proteins (review common ones from other pathogens to get ideas)

4. Retrieve amino acid sequence from the NCBI protein database

5. Predict epitopes of the protein - use Immune Epitope Database and Analysis Resource (IEDB)
    5.1. predict MHC-1 epitopes - binds with MHC-1 for all cell antigen presentation
    5.2. predict MHC-2 epitopes - binds with MHC-2 of antigen presenting cells to trigger adaptive immune response
    5.3. predict B cell epitopes - binds with B cells to initiate antibody production and immunity memory

6. Evaluate epitopes:
- allergenicity
- antigenicity
- population coverage
- conservancy
- secondary structure prediction
- protein protein interaction with MHC-1 (molecular docking)
- generation of IFN-y
- predict presence and location of signal peptides, mitochondrial, chloroplast, and thylakoid luminal transit peptides

## Tools to use:

- Machine learning:
    - feed-forward neural networks
    - restricted boltzmann machine (RBM)
    - deep learning
    - convolutional neural networks (CNNs)
    - decision trees
    - language models
    - attention mechanism
    - generative models

- NCBI bioinformatic databases (https://www.ncbi.nlm.nih.gov): to get FASTA of proteins/antigens amino acid sequence
    - BLASTp (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) - assess identity of epitopes to human proteins
- protein data bank (https://www.rcsb.org/): to get the target proteins which the antigens interact with: MHC1, MHC2, B-cell, T-cell etc
- IEDB (https://www.iedb.org/):
    - MHC-1 binding prediction tool (http://tools.iedb.org/mhci/, https://nextgen-tools.iedb.org/)
    - MHC-2 binding prediction tool (http://tools.iedb.org/mhcii/, https://nextgen-tools.iedb.org/)
    - B cell epitope prediction - multiple methods
    - T-cell epitope populatoin coverage (http://tools.iedb.org/population/)
    - conservation analysis tool (http://tools.iedb.org/conservancy/)

- Prebuilt tools:
    - NetMHCpan (artificial neural network) - HLA allele binding/interaction prediction
    - [AllerTOP v.2.1](https://www.ddg-pharmfac.net/allertop_test/) - testing epitopes allergenicity *most accurate at 88.7%*
    - [AlgPred](http://crdd.osdd.net/raghava/algpred/) - testing epitopes allergenicity
    - [AllergenFP](https://ddg-pharmfac.net/AllergenFP/) - testing epitopes allergenicity
    - [VaxiJen](https://www.ddg-pharmfac.net/vaxijen3/home/) - testing epitope immunogenicity
    - [IFNepitope](http://crdd.osdd.net/raghava/ifnepitope/) - analyse epitope capacity to generate interferon-gamma (IFN-y) *max accuracy 82.10%*
    - [SignalP](https://services.healthtech.dtu.dk/services/SignalP-5.0/)and [TargetP](https://services.healthtech.dtu.dk/services/TargetP-2.0/) - predict and localize proteins and peptides
    - [PEP-FOLD 3.0](https://bioserv.rpbs.univ-paris-diderot.fr/services/PEP-FOLD3/)
    - [TMHMM 2.0](https://services.healthtech.dtu.dk/services/TMHMM-2.0/)
    - [CLBTope](https://webs.iiitd.edu.in/raghava/clbtope/) - prediction of Linear & Conformational B-cell Epitopes
    - [Ellipro](http://tools.iedb.org/ellipro/) - predicts linear and discontinuous antibody epitopes based on a protein antigen's 3D structure (PDB file)
    - [IL-10Pred](https://webs.iiitd.edu.in/raghava/il10pred/) - Prediction of Interleukin-10 inducing peptides
    - [IL4pred](http://crdd.osdd.net/raghava/il4pred/) - Prediction of Interleukin-4 inducing peptides
    - 

- observe protein-protein interactions between HLA/MHC with epitopes (molecular docking):
    - Autodock Vina
    - Autodock
    - CrankPep (ADCP)
    - HPEPDOCK


## IEDB
> use data from all recent paper publications taken from pubMed

MHC I alleles included: 

`27 Allele panel:Covers 97% of the human population
    HLA-A*01:01, HLA-A*02:01, HLA-A*02:03, HLA-A*02:06, HLA-A*03:01, HLA-A*11:01, HLA-A*23:01, HLA-A*24:02, HLA-A*26:01, HLA-A*30:01, HLA-A*30:02, HLA-A*31:01, HLA-A*32:01, HLA-A*33:01, HLA-A*68:01, HLA-A*68:02, HLA-B*07:02, HLA-B*08:01, HLA-B*15:01, HLA-B*35:01, HLA-B*40:01, HLA-B*44:02, HLA-B*44:03, HLA-B*51:01, HLA-B*53:01, HLA-B*57:01, HLA-B*58:01, 
`

MHC II alleles included:

`7 Allele panel:Representative of the main class II human supertypes
HLA-DRB1*03:01, HLA-DRB1*07:01, HLA-DRB1*15:01, HLA-DRB3*01:01, HLA-DRB3*02:02, HLA-DRB4*01:01, HLA-DRB5*01:01, 
`

`
27 Allele panel:Covers 97% of the human population
HLA-DRB1*01:01, HLA-DRB1*03:01, HLA-DRB1*04:01, HLA-DRB1*04:05, HLA-DRB1*07:01, HLA-DRB1*08:02, HLA-DRB1*09:01, HLA-DRB1*11:01, HLA-DRB1*12:01, HLA-DRB1*13:02, HLA-DRB1*15:01, HLA-DRB3*01:01, HLA-DRB3*02:02, HLA-DRB4*01:01, HLA-DRB5*01:01, HLA-DQA1*05:01/DQB1*02:01, HLA-DQA1*05:01/DQB1*03:01, HLA-DQA1*03:01/DQB1*03:02, HLA-DQA1*04:01/DQB1*04:02, HLA-DQA1*01:01/DQB1*05:01, HLA-DQA1*01:02/DQB1*06:02, HLA-DPA1*02:01/DPB1*01:01, HLA-DPA1*01:03/DPB1*02:01, HLA-DPA1*01:03/DPB1*04:01, HLA-DPA1*03:01/DPB1*04:02, HLA-DPA1*02:01/DPB1*05:01, HLA-DPA1*02:01/DPB1*14:01 
`