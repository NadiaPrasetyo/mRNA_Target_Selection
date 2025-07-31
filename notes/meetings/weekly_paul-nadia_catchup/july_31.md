### Done:
- Implemented MAFFT using Conda package.
- Tried 4 different multiple alignment conservation score tools:
    - 1 did not install.
    - 1 did not compile.
    - 1 is dependent on another tool (i.e., same function).
    - 1 worked (patched source code AND Conda package: Rate4Site).
- Implemented Rate4Site to work post-MAFFT.
- Tested Rate4Site with test sequence and tree.

### To Do:
1. Run Rate4Site with WAG and Maximum Likelihood method.
2. Collect per-site scores and compute a conservation score:
     - Use a sliding window average to get median, max, and min.
3. Compile per-site scores alongside subcellular localization (DeepTMHMM).
4. Integrate per-site scores with subcellular localization data (Outie vs Innie):
     - Use a sliding window to compute only outside peptides.
5. Extend analysis to other bacteria:
     - Research specific strains to include.
     - Consider the number of epitopes available (e.g., from IEDB queries).
     - Account for host type.
6. Explore extending tools to viruses (and potentially cancer).