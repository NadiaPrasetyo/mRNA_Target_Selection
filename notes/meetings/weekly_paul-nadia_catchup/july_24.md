## Done

- Implemented **ifnepitope**
- Implemented **mixmhc2pred**
- Compiled results into AUROC scores
- Generated ROC curves for visualization
- Compiled AUROC scores to categorize features as useful vs. not useful
- Adding documentation on the whole project
- Share Github link with Paul and Zohaib

## To Do

- Further investigate conservation; add more conservation score computations
- Explore tools for detecting regions of positive selection:
    - Perform multiple sequence alignment (e.g., with **MAFFT** or **ClustalO**)
        - Identify highly conserved (frequency 1) and variable (frequency 0) regions
        - Use Kullback-Leibler divergence to measure deviation from expected random population
        - Apply sliding window on the alignment to assess conservation across the protein
    - Tools to consider:
        - [Consurf](https://consurf.tau.ac.il/consurf_index.php)
        - [Rate4Site](https://www.tau.ac.il/~itaymay/cp/rate4site.html)
        - [Scorecons](https://github.com/willvaldar/scorecons)
- Investigate subcellular localization further