import os
import csv
import sys
import subprocess
import tempfile
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
        return iri[8:]
    if iri.startswith("http"):
        return iri.rstrip('/').split('/')[-1]
    
    return iri

def load_antigens(antigen_file):
    antigens = {}
    antigen_names = {}  # Map antigen_id -> antigen_name

    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            acc = row['uniprot_accession']
            seq = row['sequence']
            antigens[acc] = seq

            # Assuming antigen CSV has column 'protein_name'
            antigen_name = row.get('protein_name', '').strip()
            if acc not in antigen_names and antigen_name:
                antigen_names[acc] = antigen_name
    return antigens, antigen_names

def load_epitopes(epitope_file):
    epitope_map = defaultdict(list)
    antigen_names = {}  # Map antigen_id -> antigen_name
    with open(epitope_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ep_seq = row['linear_sequence']
            raw_iri = row['parent_source_antigen_iris']
            cleaned_uniprot = clean_uniprot_iri(raw_iri)

            # Assuming epitope CSV has column 'parent_source_antigen_name'
            antigen_name = row.get('parent_source_antigen_name', '').strip()
            if cleaned_uniprot not in antigen_names and antigen_name:
                antigen_names[cleaned_uniprot] = antigen_name

            epitope_map[cleaned_uniprot].append(ep_seq)
    return epitope_map, antigen_names


def write_fasta(filehandle, sequences_dict):
    for name, seq in sequences_dict.items():
        filehandle.write(f">{name}\n{seq}\n")

def write_epitope_fasta(filehandle, epitope_seqs, antigen_id):
    for i, seq in enumerate(epitope_seqs):
        filehandle.write(f">{antigen_id}_epi_{i}\n{seq}\n")

def run_mmseqs_easy_search(query_fasta, target_fasta, output_file, tmp_dir):
    cmd = [
        'mmseqs', 'easy-search',
        query_fasta, target_fasta, output_file, tmp_dir,
        '--alignment-mode', '3',
        '--format-output', 'query,target,pident,qstart,qend,tstart,tend,qlen,tlen,evalue',
        '-s', '7.0'
    ]
    subprocess.run(cmd, check=True)

def check_mmseqs_installed():
    try:
        subprocess.run(['mmseqs'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print("❌ ERROR: mmseqs2 executable not found. Please install mmseqs2 and ensure it is in your PATH.")
        sys.exit(1)

def parse_mmseqs_results(mmseqs_result_file):
    """
    Parse mmseqs .m8 file to extract matched query IDs (epitope IDs).
    Return a set of matched query sequence IDs.
    """
    matched_queries = set()
    with open(mmseqs_result_file) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    query_id = parts[0]
                    matched_queries.add(query_id)
    return matched_queries

def main():
    if len(sys.argv) != 4:
        print('Usage: python bin/sanity_check_antigen_seq.py <pathogen_dir> "<pathogen_name>" <output_dir>')
        sys.exit(1)

    pathogen_dir = sys.argv[1]
    pathogen_name = sys.argv[2].replace(" ", "_").lower()
    user_output_dir = sys.argv[3]

    check_mmseqs_installed()

    data_dir = os.path.join("data", pathogen_dir)
    output_dir = os.path.join(data_dir, user_output_dir)

    # Check data directory exists
    if not os.path.isdir(data_dir):
        print(f"❌ ERROR: Data directory '{data_dir}' does not exist.")
        sys.exit(1)

    epitope_file = os.path.join(data_dir, f"{pathogen_name}_IEDB_epitope.csv")
    antigen_file = os.path.join(data_dir, f"{pathogen_name}_compiled_proteins.csv")

    # Check files exist
    for f in [epitope_file, antigen_file]:
        if not os.path.isfile(f):
            print(f"❌ ERROR: Required file '{f}' does not exist.")
            sys.exit(1)

    # Check/create output directory
    if not os.path.exists(output_dir):
        print(f"Output directory '{output_dir}' does not exist. Creating it...")
        os.makedirs(output_dir, exist_ok=True)

    print("[1] Loading antigen and epitope data...")
    antigens, protein_names = load_antigens(antigen_file)
    epitope_map, antigen_names = load_epitopes(epitope_file)

    summary_stats = []

    with tempfile.TemporaryDirectory() as tmp_dir:
        print(f"[2] Writing antigen FASTA in temp directory: {tmp_dir}")
        antigen_fasta_path = os.path.join(tmp_dir, f"{pathogen_name}_antigens.fasta")
        with open(antigen_fasta_path, 'w') as f:
            write_fasta(f, antigens)

        print("[3] Writing epitope FASTAs and running MMseqs2 searches...")
        for antigen_id, epitopes in epitope_map.items():
            if antigen_id not in antigens:
                print(f"  ⚠️ Antigen '{antigen_id}' not found in antigen dataset. Skipping...")
                continue

            output_file = os.path.join(output_dir, f"mmseqs_result_{antigen_id}.m8")

            if os.path.isfile(output_file):
                print(f"  ↳ MMseqs results for antigen '{antigen_id}' already exist. Skipping search and analyzing existing results...")
            else:
                query_fasta_path = os.path.join(tmp_dir, f"{antigen_id}_epitopes.fasta")
                with open(query_fasta_path, 'w') as f:
                    write_epitope_fasta(f, epitopes, antigen_id)

                print(f"  ↳ Searching epitopes of antigen '{antigen_id}'...")
                try:
                    run_mmseqs_easy_search(query_fasta_path, antigen_fasta_path, output_file, tmp_dir)
                except subprocess.CalledProcessError:
                    print(f"  ❌ MMseqs search failed for {antigen_id}")
                    continue

                print(f"    Results saved: {output_file}")

            # Analyze results
            matched_queries = parse_mmseqs_results(output_file)

            total_epitopes = len(epitopes)
            mapped_count = 0
            unmapped_epitopes = []

            for i, seq in enumerate(epitopes):
                epitope_id = f"{antigen_id}_epi_{i}"
                if epitope_id in matched_queries:
                    mapped_count += 1
                else:
                    unmapped_epitopes.append(seq)

            not_mapped_count = total_epitopes - mapped_count

            summary_stats.append({
                'antigen_id': antigen_id,
                'antigen_name': antigen_names.get(antigen_id, protein_names.get(antigen_id, '')),
                'expected_epitopes': total_epitopes,
                'mapped_epitopes': mapped_count,
                'not_mapped_epitopes': not_mapped_count,
                'unmapped_epitope_seqs': ";".join(unmapped_epitopes) if unmapped_epitopes else ""
            })

    # Add antigens with no epitopes expected, with zero counts
    proteins_no_epitopes = [acc for acc in antigens.keys() if acc not in epitope_map]
    for acc in proteins_no_epitopes:
        summary_stats.append({
            'antigen_id': acc,
            'antigen_name': protein_names.get(acc, ''),
            'expected_epitopes': 0,
            'mapped_epitopes': 0,
            'not_mapped_epitopes': 0,
            'unmapped_epitope_seqs': ''
        })

    # Write summary CSV
    summary_csv_path = os.path.join(output_dir, f"{pathogen_name}_epitope_mapping_summary.csv")
    print(f"\n[4] Writing epitope mapping summary to: {summary_csv_path}")

    with open(summary_csv_path, 'w', newline='') as csvfile:
        fieldnames = ['antigen_id', 'antigen_name', 'expected_epitopes', 'mapped_epitopes', 'not_mapped_epitopes', 'unmapped_epitope_seqs']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_stats:
            writer.writerow(row)

    print("\n✅ All searches complete. Summary file created. Temporary files cleaned up automatically.")
    
if __name__ == "__main__":
    main()
