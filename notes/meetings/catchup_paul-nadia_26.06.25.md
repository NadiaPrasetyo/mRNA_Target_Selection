Have done:
- requested MARIA software access (only 6 months access) - waiting for response
- removed TMHMM analysis
- replaced Cluster with mmseq2 cluster
- got cluster results for randoms
- patching algpred - Allergenicity (figure out)
- collecting epitope evaluation data on the randoms

sort alphabetically: sort -d
minimum epitopes length for cluster and analysis: 5

Order to do things in:
- take the whole proteome/ protein sequences
- cluster through those for the different strains: cluster full protein seq first
- get the epitopes for the proteins
- see if the epitopes are conserved: grep/mmseq2 search to find whether or not each of the cluster have the epitope subsequence
  => may safe from pain later
  1. priority: are the proteins conserved -> split the sets to conserved vs not
  2. are the epitopes conserved