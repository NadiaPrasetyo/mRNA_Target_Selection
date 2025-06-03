import os
import csv
import sys
import subprocess
from collections import defaultdict

def clean_uniprot_iri(iri):
    """
    Convert parent_source_antigen_iri like 'UNIPROT:O80066' or
    'http://identifiers.org/uniprot/O80066' to 'O80066' format.
    """
    # Strip list brackets or quotes if present
    # For this example, assume you get a string: "['UNIPROT:O80066']"
    # or a direct string like "UNIPROT:O80066"
        # Remove enclosing brackets and quotes if present
    iri = iri.strip("[]'\" ")
    
    # If it's still comma separated (multiple iris?), take first only:
    iri = iri.split(",")[0].strip("[]'\" ")
    
    # Remove prefix if present
    if iri.upper().startswith("UNIPROT:"):
        return iri[8:]  # Remove 'UNIPROT:' prefix
    
    # If URI format
    if iri.startswith("http"):
        return iri.rstrip('/').split('/')[-1]
    
    return iri

def load_antigens(antigen_file):
    antigens = {}
    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            acc = row['uniprot_accession']
            seq = row['sequence']
            antigens[acc] = seq
    return antigens

def load_epitopes(epitope_file):
    epitope_map = defaultdict(list)
    with open(epitope_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ep_seq = row['linear_sequence']
            raw_iris = row['parent_source_antigen_iris']
            # Some IRIs may be stored as a stringified list, e.g. "['UNIPROT:O80066']"
            # You may want to parse that properly or just handle a single one per epitope:
            # Extract the first IRI cleanly
            cleaned_uniprot = clean_uniprot_iri(raw_iris)
            epitope_map[cleaned_uniprot].append(ep_seq)

    return epitope_map

def write_fasta(fasta_path, sequences_dict):
    with open(fasta_path, 'w') as f:
        for name, seq in sequences_dict.items():
            f.write(f">{name}\n{seq}\n")

def write_epitope_fasta(epitope_map, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    fasta_paths = {}
    for antigen_id, epitopes in epitope_map.items():
        epitope_fasta = os.path.join(output_dir, f"{antigen_id}_epitopes.fasta")
        with open(epitope_fasta, 'w') as f:
            for i, seq in enumerate(epitopes):
                f.write(f">{antigen_id}_epi_{i}\n{seq}\n")
        fasta_paths[antigen_id] = epitope_fasta
    return fasta_paths

def run_mmseqs_easy_search(query_fasta, target_fasta, output_file, tmp_dir):
    cmd = [
        'mmseqs', 'easy-search',
        query_fasta, target_fasta, output_file, tmp_dir,
        '--alignment-mode', '3',
        '--format-output', 'query,target,fident,evalue,qaln,taln',
        '-s', '7.0'
    ]
    subprocess.run(cmd, check=True)

def main():
    if len(sys.argv) != 3:
        print("Usage: python bin/sanity_check_antigen_seq.py <pathogen_dir> <pathogen_name>")
        sys.exit(1)

    pathogen_dir = sys.argv[1]
    pathogen_name = sys.argv[2]
    data_dir = os.path.join("data", pathogen_dir)

    if not os.path.isdir(data_dir):
        print(f"Error: Directory '{data_dir}' does not exist.")
        sys.exit(1)

    epitope_file = os.path.join(data_dir, f"{pathogen_name}_IEDB_epitope.csv")
    antigen_file = os.path.join(data_dir, f"{pathogen_name}_compiled_proteins.csv")

    for file in [epitope_file, antigen_file]:
        if not os.path.isfile(file):
            print(f"Error: Required file '{file}' does not exist.")
            sys.exit(1)

    tmp_dir = os.path.join(data_dir, 'mmseqs_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    print("[1] Loading antigen and epitope data...")
    antigens = load_antigens(antigen_file)
    epitope_map = load_epitopes(epitope_file)

    print("[2] Writing antigen FASTA...")
    antigen_fasta = os.path.join(data_dir, f"{pathogen_name}_antigens.fasta")
    write_fasta(antigen_fasta, antigens)

    print("[3] Writing epitope FASTAs...")
    epitope_fastas = write_epitope_fasta(epitope_map, os.path.join(data_dir, "epitope_fastas"))

    print("[4] Running MMseqs2 FASTA-to-FASTA searches...")
    for antigen_id, query_fasta in epitope_fastas.items():
        out_file = os.path.join(data_dir, f"mmseqs_results_{antigen_id}.m8")
        print(f"  ↳ Searching for epitopes of {antigen_id}...")
        try:
            run_mmseqs_easy_search(query_fasta, antigen_fasta, out_file, tmp_dir)
        except subprocess.CalledProcessError:
            print(f"  ❌ Search failed for {antigen_id}")

    print("\n✅ MMseqs2 search complete. Results saved in:")
    for antigen_id in epitope_fastas:
        print(f"  - mmseqs_results_{antigen_id}.m8")

if __name__ == "__main__":
    main()
