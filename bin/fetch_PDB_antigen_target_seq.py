# 6. Fetch the amino acid sequence for the targets of the antigens:
#     - get the PDB identifier of each protein (MHC1, MHC2, T-cell, B-cell)
#     - use search query: https://search.rcsb.org/rcsbsearch/v2/query
#       - Tutorial: https://search.rcsb.org/index.html#search-api
#         - return type: entry (list of PDB IDs)
#         - query language: full query DSL (domain-specific language) based on JSON
#       - `In GET method, search request should be sent as a URL-encoded query string in json parameter: https://search.rcsb.org/rcsbsearch/v2/query?json={search-request}.`
#     - use: https://www.rcsb.org/docs/programmatic-access/web-apis-overview to fetch aa seq.
#     - find: "MHC I", "MHC II", "T-cell", "B-cell" in the PDB entry name or description
#     - e.g. QUERY: Full Text = "T cell" AND ( Scientific Name of the Source Organism = "Homo sapiens" AND Polymer Entity Type = "Protein" )
#     - batch: QUERY: ( Full Text = "MHC II" OR Full Text = "MHC I" OR Full Text = "T cell" OR Full Text = "B cell" ) AND ( ( Scientific Name of the Source Organism = "Homo sapiens" AND Polymer Entity Type = "Protein" ) AND ( Scientific Name of the Source Organism = "Homo sapiens" AND Polymer Entity Type = "Protein" ) )
#       => 62,561 Structures???
#     - Polymer Entity Type: "Protein"