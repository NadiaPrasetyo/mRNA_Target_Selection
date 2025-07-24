"""
visualise_cluster.py

Command-line tool to visualize protein clusters across strains from TSV files.

Overview:
    - Parses cluster TSV files and corresponding FASTA files for protein sequences.
    - Generates a network graph of clusters using NetworkX and saves as a PNG image.
    - Optionally splits clusters into separate FASTA files for each connected component.
    - Extracts sequence records and formats node labels for visualization.
    - Supports verbose logging to file.

Arguments:
    pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
    sequence_dir (str): Subdirectory under `data/<pathogen_dir>/` with sequence/cluster files.
    --threads (int, optional): Number of parallel workers (default: 4).
    --output-dir (Path, optional): Output directory for results (default: epitope_outputs).
    --verbose (flag, optional): If set, enables verbose logging to file.
    --split-clusters (flag, optional): If set, writes separate FASTA files for each cluster.

Requirements:
    - Input files in expected formats under:
        data/<pathogen_dir>/<sequence_dir>/*_combined_clu.tsv
        data/<pathogen_dir>/<sequence_dir>/*_combined_clu.fasta
    - Python packages: argparse, logging, pathlib, networkx, matplotlib, Bio (Biopython)

Usage Example:
    python visualise_cluster.py sars_cov_2 protein_clusters --threads 4 --verbose --split-clusters

Outputs:
    data/<pathogen_dir>/<output-dir>/<antigen_id>_clusters.png      # Cluster network image
    data/<pathogen_dir>/<output-dir>/<antigen_id>_cluster_*.fasta   # (optional) FASTA files per cluster
    data/<pathogen_dir>/<output-dir>/cluster_viz.log                # (optional) Verbose log file

Author: Nadia
"""
import argparse
import logging
from pathlib import Path
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
from Bio import SeqIO
import random

def setup_logging(output_dir: Path, verbose: bool):
    """Set up logging configuration.
    Args:
        output_dir (Path): Directory to store log file.
        verbose (bool): If True, enable verbose logging."""
    log_file = output_dir / "cluster_viz.log"
    log_level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=log_level,
        format='[%(asctime)s] %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

    # Silence noisy third-party loggers
    for noisy_logger in ["matplotlib", "PIL"]:
        logging.getLogger(noisy_logger).setLevel(logging.WARNING)

    logging.info("Logging initialized.")

def parse_fasta(fasta_file):
    """Parse a FASTA file and return a dictionary of sequence records.
    Args:
        fasta_file (Path): Path to the FASTA file.
    Returns:
        dict: Dictionary mapping sequence IDs to SeqRecord objects.
    """
    records = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        records[record.id] = record
    logging.info(f"Parsed {len(records)} sequences from {fasta_file}")
    return records


def parse_clusters(tsv_file):
    """Parse a TSV file containing cluster information.
    Args:
        tsv_file (Path): Path to the TSV file.
    Returns:
        dict: Dictionary mapping representative sequence IDs to sets of member IDs.
    """
    clusters = defaultdict(set)
    with open(tsv_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            rep, member = line.strip().split("\t")
            clusters[rep].add(member)
    logging.info(f"Parsed {len(clusters)} clusters from {tsv_file}")
    return clusters


def get_label(seq_id):
    """Generate a label for a sequence ID.
    Args:
        seq_id (str): Sequence ID from the FASTA file.
    Returns:
        str: Formatted label string.
    """
    parts = seq_id.replace('|', '_').split('_')
    strain = parts[0]
    position = parts[1] if len(parts) > 1 else "unknown"
    return f"{strain}:{position}"


def split_fasta_by_cluster(clusters, records, output_dir, antigen_id):
    """Split clusters into separate FASTA files for each connected component.
    Args:
        clusters (dict): Dictionary mapping representative IDs to sets of member IDs.
        records (dict): Dictionary of sequence records from the FASTA file.
        output_dir (Path): Directory to save split FASTA files.
        antigen_id (str): Antigen identifier for naming files.
    """
    for rep, members in clusters.items():
        cluster_file = output_dir / f"{antigen_id}_cluster_{rep}.fasta"
        cluster_seqs = [records[seq_id] for seq_id in members if seq_id in records]
        if rep not in members and rep in records:
            cluster_seqs.insert(0, records[rep])
        SeqIO.write(cluster_seqs, cluster_file, "fasta")
        logging.debug(f"Wrote {len(cluster_seqs)} sequences to {cluster_file}")

def format_node_label(raw_id):
    """Format node label for visualization.
    Args:
        raw_id (str): Raw sequence ID from the FASTA file.
    Returns:
        str: Formatted label string.
    """
    try:
        parts = raw_id.split('|')
        strain_id = parts[3]                  # e.g., BA000018.3
        position = parts[4].replace('tpos:', '')  # e.g., 637462-637595
        return f"{strain_id}\n{position}"
    except IndexError:
        return raw_id  # fallback if unexpected format


def visualize_clusters(clusters, records, output_path: Path, antigen_id: str):
    """Visualize protein clusters as a network graph and save as PNG.
    Args:
        clusters (dict): Dictionary mapping representative IDs to sets of member IDs.
        records (dict): Dictionary of sequence records from the FASTA file.
        output_path (Path): Directory to save the visualization image.
        antigen_id (str): Antigen identifier for naming the output file.
    """
    G = nx.Graph()
    color_map = {}
    cluster_colors = {}

    def get_random_color():
        return "#" + "".join(random.choices("0123456789ABCDEF", k=6))

    for rep, members in clusters.items():
        color = get_random_color()
        cluster_colors[rep] = color
        for member in members:
            label = format_node_label(member)
            G.add_node(member, label=label)
            G.add_edge(rep, member)
            color_map[member] = color
        color_map[rep] = color
        if rep not in members:
            G.add_node(rep, label=get_label(rep))

    pos = nx.spring_layout(G, seed=42)
    node_colors = [color_map[node] for node in G.nodes()]
    node_labels = nx.get_node_attributes(G, "label")

    plt.figure(figsize=(18, 12))
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500, alpha=0.8)
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=20)

    plt.title(f"Protein Clusters for {antigen_id}")
    plt.axis("off")
    plt.tight_layout()
    output_img = output_path / f"{antigen_id}_clusters.png"
    plt.savefig(output_img, dpi=300)
    plt.close()
    logging.info(f"Saved cluster visualization to {output_img}")


def main():
    """Main function to parse arguments and run the visualization."""
    parser = argparse.ArgumentParser(description="Visualize protein clusters across strains from TSVs.")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output-dir", type=Path, default=Path("epitope_outputs"))
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--split-clusters", action="store_true", help="Split clusters into separate files for each connected component.")
    args = parser.parse_args()

    input_dir = Path("data") / args.pathogen_dir / args.sequence_dir
    output_dir = Path("data") / args.pathogen_dir / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(output_dir, args.verbose)

    fasta_files = list(input_dir.glob("*_combined_clu.fasta"))
    tsv_files = list(input_dir.glob("*_combined_clu.tsv"))

    # Build map: antigen_id -> (tsv, fasta)
    antigen_map = {}
    for fasta in fasta_files:
        antigen_id = fasta.stem.replace("_combined_clu", "")
        matching_tsv = input_dir / f"{antigen_id}_combined_clu.tsv"
        if matching_tsv.exists():
            antigen_map[antigen_id] = (matching_tsv, fasta)

    logging.info(f"Found {len(antigen_map)} antigen cluster files.")

    for antigen_id, (tsv_file, fasta_file) in antigen_map.items():
        logging.info(f"Processing antigen: {antigen_id}")
        records = parse_fasta(fasta_file)
        clusters = parse_clusters(tsv_file)

        if args.split_clusters:
            split_fasta_by_cluster(clusters, records, output_dir, antigen_id)

        visualize_clusters(clusters, records, output_dir, antigen_id)

    logging.info("All antigens processed.")


if __name__ == "__main__":
    main()
