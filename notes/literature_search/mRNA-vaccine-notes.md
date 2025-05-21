Existing mRNA Vaccines used for:
- COVID-19: BNT162b2, 
- Zika virus
- Influenza

Pipeline steps for epitope selection, prediction, and evaluation:
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

Tools to use:
- Machine learning:
    - feed-forward neural networks
    - restricted boltzmann machine (RBM)
    - deep learning
    - convolutional neural networks (CNNs)
    - decision trees
    - language models
    - attention mechanism
    - generative models

- NCBI bioinformatic databases (https://www.ncbi.nlm.nih.gov):
    - BLASTp (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) - assess identity of epitopes to human proteins
- protein data bank (https://www.rcsb.org/)
- IEDB (https://www.iedb.org/):
    - MHC-1 binding prediction tool (http://tools.iedb.org/mhci/, https://nextgen-tools.iedb.org/)
    - MHC-2 binding prediction tool (http://tools.iedb.org/mhcii/, https://nextgen-tools.iedb.org/)
    - B cell epitope prediction - multiple methods
    - T-cell epitope populatoin coverage (http://tools.iedb.org/population/)
    - conservation analysis tool (http://tools.iedb.org/conservancy/)

- Prebuilt tools:
    - NetMHCpan (artificial neural network) - HLA allele binding/interaction prediction
    - AllerTOP v.2.1 (https://www.ddg-pharmfac.net/allertop_test/) - testing epitopes allergenicity *most accurate at 88.7%*
    - AlgPred (http://crdd.osdd.net/raghava/algpred/) - testing epitopes allergenicity
    - AllergenFP (https://ddg-pharmfac.net/AllergenFP/) - testing epitopes allergenicity
    - VaxiJen (https://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html, https://www.ddg-pharmfac.net/vaxijen3/home/) - testing epitope immunogenicity
    - IFNepitope (http://crdd.osdd.net/raghava/ifnepitope/) - analyse epitope capacity to generate interferon-gamma (IFN-y) max accuracy 82.10%
    - SignalP and TargetP (https://services.healthtech.dtu.dk/services/SignalP-5.0/, https://services.healthtech.dtu.dk/services/TargetP-2.0/) - predict and localize proteins and peptides
    - PEP-FOLD 3.0 (https://bioserv.rpbs.univ-paris-diderot.fr/services/PEP-FOLD3/)

observe protein-protein interactions between HLA/MHC with epitopes (molecular docking):
    - Autodock Vina
    - Autodock
    - CrankPep (ADCP)
    - HPEPDOCK