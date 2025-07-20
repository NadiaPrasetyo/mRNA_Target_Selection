Done:
- implemented and ran deeplocpro (local) on randoms and positive set
- corrected allergenicity to include non-allergenic peptides/antigens in the results
- fetch PDB sequences of antigens for conformation based predictions
- implement ellipro for Bcell conformational predicted epitopes

Questions:
- IFNepitope, a tool to predict the induction of IFN-gamma - 2 options: 
  - full antigens - using their "scan" function that uses a 15-20 aa slider to predict motifs that induce IFN gamma
  - epitopes - using their "predict" function, calculates the hybrid scores of peptides and predicts whether they are IFN-gamma inducing or not

- alternative tools - compare the auc distances between the negatice and positive sets (would tell me which tool is the most suitable tool)

To do:
- Create  feature matrix: each row is antigen, each column is a label/feature
- create an auc/roc curves for each (between positive or negative) - 
- figure out the correlation of the different features 
- get a cumulative scores, write a random forest distributions
- use coronaviruses, influenza, hepatitis as a test set/ training data; other bacteria as well: pseudomonas, salmonella, clamydia, helicobacter pylori, 
- look into the ancestral sequence accross the strains


other things:
- NZ RNA Symposium 10th November Monday