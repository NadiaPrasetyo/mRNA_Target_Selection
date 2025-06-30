import argparse
import logging
from pathlib import Path
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
from Bio import SeqIO
import random

def setup_logging(output_dir: Path, verbose: bool):
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
    records = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        records[record.id] = record
    logging.info(f"Parsed {len(records)} sequences from {fasta_file}")
    return records


def parse_clusters(tsv_file):
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
    parts = seq_id.replace('|', '_').split('_')
    strain = parts[0]
    position = parts[1] if len(parts) > 1 else "unknown"
    return f"{strain}:{position}"


def split_fasta_by_cluster(clusters, records, output_dir, antigen_id):
    for rep, members in clusters.items():
        cluster_file = output_dir / f"{antigen_id}_cluster_{rep}.fasta"
        cluster_seqs = [records[seq_id] for seq_id in members if seq_id in records]
        if rep not in members and rep in records:
            cluster_seqs.insert(0, records[rep])
        SeqIO.write(cluster_seqs, cluster_file, "fasta")
        logging.debug(f"Wrote {len(cluster_seqs)} sequences to {cluster_file}")

def format_node_label(raw_id):
    try:
        parts = raw_id.split('|')
        strain_id = parts[3]                  # e.g., BA000018.3
        position = parts[4].replace('tpos:', '')  # e.g., 637462-637595
        return f"{strain_id}\n{position}"
    except IndexError:
        return raw_id  # fallback if unexpected format


def visualize_clusters(clusters, records, output_path: Path, antigen_id: str):
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
