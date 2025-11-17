#!/usr/bin/env python3
"""
filter_proteins_with_auto_human.py

Filter raw CSV, merge with predictions, and remove proteins with human
similarity using BLASTP. If --human-fasta is omitted, the script will attempt
to locate the human proteome in data/human_proteome/. If it does not exist,
the script automatically downloads the human proteome (RefSeq: GCF_000001405.40)
using helpers from bin/fetch_NCBI_strain_genomes.py.

Requirements:
    - BLAST+ installed (makeblastdb, blastp)
    - requests, tqdm (used by the fetch helper)
    - bin/fetch_NCBI_strain_genomes.py must be present and importable

Usage:
    python filter_proteins_with_auto_human.py --input-raw raw.csv \
        --input-pred preds.csv --input-fasta proteins.fasta
    # optionally:
    --human-fasta path/to/human_proteome.fasta
"""

import argparse
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))  # Add parent directory to sys.path


# Attempt to import the NCBI helper functions from the provided script path
# The user indicated the fetch script is at: bin/fetch_NCBI_strain_genomes.py
try:
    from bin.fetch_NCBI_strain_genomes import (
        get_requests_session,
        download_and_extract_zip,
        ensure_dir,
        API_BASE,
    )
except Exception as e:
    # Provide helpful error message if import fails
    print(
        "ERROR: Failed to import required functions from bin/fetch_NCBI_strain_genomes.py.\n"
        "Make sure the file exists at bin/fetch_NCBI_strain_genomes.py and that the current\n"
        "working directory is the project root. Import error details:\n",
        e,
        file=sys.stderr,
    )
    raise


# -----------------------------------------------------------
# Small wrapper: fetch_human_proteome
# (reuses download_and_extract_zip and get_requests_session
#  from bin/fetch_NCBI_strain_genomes)
# -----------------------------------------------------------

def fetch_human_proteome(refseq_id: str, output_dir: str, session=None) -> str:
    """
    Download the human proteome FASTA for the given RefSeq assembly using
    NCBI Datasets v2 API and the project's download helper.

    Parameters:
        refseq_id: e.g. "GCF_000001405.40"
        output_dir: directory to save extracted FASTA (folder will be created)
        session: optional requests.Session (if not provided, get_requests_session() will be used)

    Returns:
        Path (string) to the extracted protein FASTA file, typically:
            {output_dir}/{refseq_id}_proteins.fasta

    Raises:
        FileNotFoundError if the expected FASTA cannot be located after download.
    """
    # Ensure output dir exists
    ensure_dir(output_dir)

    session = session or get_requests_session()

    # Build the download URL in the same style used in the fetch script
    # The fetch script uses: f"{API_BASE}/accession/{acc}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=PROT_FASTA&hydrated=FULLY_HYDRATED"
    # For proteins-only, include PROT_FASTA and hydrated
    url = (
        f"{API_BASE}/accession/{refseq_id}/download?"
        "include_annotation_type=PROT_FASTA&hydrated=FULLY_HYDRATED"
    )

    print(f"📥 Downloading human proteome for {refseq_id} ...")
    # download_and_extract_zip will extract matching files and rename them according to its logic:
    # - .faa -> {accession}_proteins.fasta
    # - .fna -> {accession}.fasta
    # so we expect: {refseq_id}_proteins.fasta
    download_and_extract_zip(url, refseq_id, output_dir, session=session)

    expected = os.path.join(output_dir, f"{refseq_id}_proteins.fasta")
    if not os.path.exists(expected):
        # try a couple of plausible alternative names (be robust)
        alt1 = os.path.join(output_dir, f"{refseq_id}_proteins.faa")
        alt2 = os.path.join(output_dir, f"{refseq_id}.faa")
        alt3 = os.path.join(output_dir, f"{refseq_id}.fasta")
        for alt in (alt1, alt2, alt3):
            if os.path.exists(alt):
                print(f"ℹ️  Found alternative extracted file: {alt}")
                return alt
        raise FileNotFoundError(
            f"Expected protein FASTA not found after download: {expected}. "
            f"Check download logs or that the NCBI endpoint is available."
        )

    print(f"✅ Human proteome saved to: {expected}")
    return expected


# -----------------------------------------------------------
# BLASTP Utilities (unchanged style from your original script)
# -----------------------------------------------------------

def run_blastp(seq, human_db_path, tmpdir):
    """Run BLASTP for a single amino acid sequence and return best hit stats."""
    query_fa = Path(tmpdir) / "query.fa"
    with open(query_fa, "w") as f:
        f.write(">query\n" + seq + "\n")

    out_path = Path(tmpdir) / "blast_out.tsv"
    cmd = [
        "blastp",
        "-query", str(query_fa),
        "-db", str(human_db_path),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-max_target_seqs", "1"
    ]
    subprocess.run(cmd, stdout=open(out_path, "w"), stderr=subprocess.DEVNULL)

    if out_path.stat().st_size == 0:
        return None  # no hits

    df = pd.read_csv(
        out_path, sep="\t", header=None,
        names=["qseqid","sseqid","pident","length","mismatch","gapopen",
               "qstart","qend","sstart","send","evalue","bitscore"]
    )
    return df.iloc[0]  # best hit


def ensure_human_db(human_fasta):
    """Build BLAST DB for human proteome if not already built."""
    required = [human_fasta + ext for ext in [".pin", ".phr", ".psq"]]
    if all(Path(f).exists() for f in required):
        return human_fasta  # already exists

    print("⚙️  Building BLAST database for human proteome...")
    subprocess.run(["makeblastdb", "-in", human_fasta, "-dbtype", "prot"], check=True)
    return human_fasta


def load_fasta_to_dict(fasta_path):
    """Load FASTA into dict: accession → sequence."""
    seqs = {}
    with open(fasta_path) as f:
        acc = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(seq_lines)
                acc = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if acc:
            seqs[acc] = "".join(seq_lines)
    return seqs


# -----------------------------------------------------------
# Main
# -----------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter raw CSV, merge with predictions, and remove proteins with human similarity (BLASTP)."
    )
    parser.add_argument("--input-raw", required=True)
    parser.add_argument("--input-pred", required=True)
    parser.add_argument("--input-fasta", required=True,
                        help="FASTA containing sequences corresponding to accessions in the dataset.")
    parser.add_argument("--human-fasta", default=None,
                        help="FASTA of all human proteins. If omitted the script will look in data/human_proteome/ and download GCF_000001405.40 if missing.")
    parser.add_argument("-o", "--output-file", default=None)

    args = parser.parse_args()

    # Paths
    raw_path = Path(args.input_raw)
    pred_path = Path(args.input_pred)
    fasta_path = Path(args.input_fasta)

    output_path = (
        Path(args.output_file)
        if args.output_file
        else Path("results") / f"filtered_{pred_path.stem}.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load data
    df_raw = pd.read_csv(raw_path)
    df_pred = pd.read_csv(pred_path)

    # Load sequences
    seq_dict = load_fasta_to_dict(fasta_path)

    # Determine human proteome FASTA
    if args.human_fasta:
        human_fasta = Path(args.human_fasta)
        if not human_fasta.exists():
            raise FileNotFoundError(f"Provided --human-fasta does not exist: {human_fasta}")
        print(f"📂 Using provided human FASTA: {human_fasta}")
    else:
        cache_dir = Path("data/human_proteome")
        cache_dir.mkdir(parents=True, exist_ok=True)
        expected = cache_dir / "GCF_000001405.40_proteins.fasta"

        if expected.exists():
            human_fasta = expected
            print(f"📂 Using cached human proteome: {human_fasta}")
        else:
            # Download
            print("📥 Human proteome FASTA not found in cache; downloading now...")
            human_fasta_path = fetch_human_proteome("GCF_000001405.40", str(cache_dir))
            human_fasta = Path(human_fasta_path)

    # Build BLAST DB if necessary
    human_db = ensure_human_db(str(human_fasta))

    # BLASTP filtering
    print("🔍 Running BLASTP similarity filtering against human proteome...")
    blast_flags = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for acc in df_raw["accession"]:
            seq = seq_dict.get(acc, None)
            if seq is None:
                blast_flags.append(False)  # cannot evaluate, keep
                continue

            hit = run_blastp(seq, human_db, tmpdir)
            if hit is None:
                blast_flags.append(False)  # no hit → keep
                continue

            # Filter: >30% identity AND E-value < 0.005
            remove_flag = (hit["pident"] > 30.0) and (hit["evalue"] < 0.005)
            blast_flags.append(remove_flag)

    df_raw["remove_human_similarity"] = blast_flags
    removed_human = df_raw[df_raw["remove_human_similarity"] == True]
    df_raw = df_raw[df_raw["remove_human_similarity"] == False]

    print(f"❌ Removed {len(removed_human)} for high similarity to human proteins.")

    # --- Your existing filtering ---
    filtered_raw = df_raw[
        (df_raw["allergenicity_hybrid_score"] < 0.3) &
        (df_raw["mhci_num_strong_binders"] > 0) &
        (df_raw["mhcii_num_strong_binders"] > 0)
    ]

    removed = df_raw[~df_raw.index.isin(filtered_raw.index)]

    # Merge
    merged = pd.merge(
        filtered_raw, df_pred, on="accession", how="inner",
        suffixes=("_raw", "_pred")
    )

    merged_removed = pd.merge(
        removed, df_pred, on="accession", how="inner",
        suffixes=("_raw", "_pred")
    )

    prob_col = "prob_antigen_raw" if "prob_antigen_raw" in merged.columns else "prob_antigen"

    final_df = merged[["accession", prob_col, "pred_label", "protein_names", "gene_names"]]
    removed_df = merged_removed[["accession", prob_col, "pred_label", "protein_names", "gene_names"]]

    # Save
    final_df.to_csv(output_path, index=False)
    removed_output_path = output_path.parent / f"removed_{pred_path.stem}.csv"
    removed_df.to_csv(removed_output_path, index=False)

    # Save also human-similarity-removed list
    removed_human_path = output_path.parent / f"removed_human_similarity_{pred_path.stem}.csv"
    removed_human.to_csv(removed_human_path, index=False)

    print(f"✅ Final filtered output saved to: {output_path}")
    print(f"❌ Human similarity removed saved to: {removed_human_path}")


if __name__ == "__main__":
    main()
