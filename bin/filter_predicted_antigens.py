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
import logging
import re
from collections import defaultdict, deque
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
    logging.exception(
        "ERROR: Failed to import required functions from bin/fetch_NCBI_strain_genomes.py.\n"
        "Make sure the file exists at bin/fetch_NCBI_strain_genomes.py and that the current\n"
        "working directory is the project root. Import error details:"
    )
    raise


def setup_logging(verbose: bool, output_path: Path):
    """Configure logging: always log to console; if verbose also append to a file."""
    logger = logging.getLogger()
    # Clear existing handlers to avoid duplicate logs when reusing in interactive runs
    for h in list(logger.handlers):
        logger.removeHandler(h)
    level = logging.DEBUG if verbose else logging.INFO
    logger.setLevel(level)
    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    # Console handler (always)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    # File handler (only when verbose)
    if verbose:
        log_file = output_path / "filter_predicted_antigens.log"
        fh = logging.FileHandler(log_file, mode="a")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    logger.debug("Logging initialized (verbose=%s)", verbose)

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

    logging.info(f"📥 Downloading human proteome for {refseq_id} ...")
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
                logging.info(f"ℹ️  Found alternative extracted file: {alt}")
                return alt
        raise FileNotFoundError(
            f"Expected protein FASTA not found after download: {expected}. "
            f"Check download logs or that the NCBI endpoint is available."
        )

    logging.info(f"✅ Human proteome saved to: {expected}")
    return expected


# -----------------------------------------------------------
# BLASTP Utilities (unchanged style from your original script)
# -----------------------------------------------------------

def run_blastp(accession, seq, human_db_path, tmpdir):
    """Run BLASTP for a single amino acid sequence and return best hit stats."""
    query_fa = Path(tmpdir) / "query.fa"
    with open(query_fa, "w") as f:
        f.write(f">{accession}\n" + seq + "\n")

    out_path = Path(tmpdir) / f"{accession}_blast_out.tsv"
    cmd = [
        "blastp",
        "-query", str(query_fa),
        "-db", str(human_db_path),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-max_target_seqs", "1"
    ]
    subprocess.run(cmd, stdout=open(out_path, "w"), stderr=subprocess.DEVNULL)

    if out_path.stat().st_size == 0:
        logging.info(f"No BLASTP hits for {accession}.")
        out_path.unlink()
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

    logging.info("⚙️  Building BLAST database for human proteome...")
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
                    print(f"adding seq to {acc}")
                acc = line[1:].split('|')[1]
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
    parser.add_argument("--verbose", action="store_true")

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

    setup_logging(args.verbose, output_path.parent)

    # Load data
    df_raw = pd.read_csv(raw_path)
    df_pred = pd.read_csv(pred_path)

    # -----------------------------------------------------------
    # Collapse homologous proteins using BOTH protein_names and gene_names
    # -----------------------------------------------------------

    def extract_codes(name_string):
        """
        Extract any 4-letter alphanumeric codes from protein_names or gene_names.
        E.g.: 'sasC', 'ebh1', 'mrpA', etc.
        """
        if pd.isna(name_string):
            return []
        parts = re.split(r"[ ,;/()]+", name_string)
        codes = [p for p in parts if re.fullmatch(r"[A-Za-z0-9]{4}", p)]
        return [c.lower() for c in codes]


    # 1. Build code → accession sets for BOTH protein_names AND gene_names
    code_groups = defaultdict(set)

    for idx, row in df_pred.iterrows():
        acc = row["accession"]

        # protein-name codes
        prot_codes = extract_codes(row.get("protein_names", ""))

        # gene-name codes (NEW)
        gene_codes = extract_codes(row.get("gene_names", ""))

        for c in prot_codes + gene_codes:
            code_groups[c].add(acc)


    # 2. Merge overlapping sets (Union-Find–like merging)
    groups = []
    for _, accs in code_groups.items():
        accs = set(accs)
        if len(accs) <= 1:
            continue

        merged = False
        for g in groups:
            if g & accs:       # overlap → merge
                g |= accs
                merged = True
                break

        if not merged:
            groups.append(accs)


    # 3. Create collapse map: accession → representative
    collapse_map = {}     # accession → rep
    alt_map = {}          # rep → [alts]


    for cluster in groups:
        cluster = list(cluster)

        # Choose representative = highest prob_antigen
        sub = df_pred[df_pred["accession"].isin(cluster)]
        rep = sub.loc[sub["prob_antigen"].idxmax()]["accession"]

        others = [a for a in cluster if a != rep]
        alt_map[rep] = others

        for o in others:
            collapse_map[o] = rep

    # Load FASTA sequences
    logging.info("📥 Loading input FASTA sequences...")
    seq_dict = load_fasta_to_dict(fasta_path)

    # 4. Apply collapse to df_pred, df_raw and integrate names/means
    if collapse_map:
        logging.info(f"Collapsing homologous proteins (protein+gene based): {len(collapse_map)} removed.")

        # Save removed homologs list
        homolog_removed = pd.DataFrame({
            "accession": list(collapse_map.keys()),
            "representative": [collapse_map[a] for a in collapse_map.keys()]
        })
        homolog_removed_path = output_path.parent / f"removed_homologs_{pred_path.stem}.csv"
        homolog_removed.to_csv(homolog_removed_path, index=False)
        logging.info(f"Homolog-removed proteins saved to: {homolog_removed_path}")

        # Remove collapsed accessions
        df_pred = df_pred[~df_pred["accession"].isin(collapse_map.keys())].copy()
        df_raw  = df_raw[~df_raw["accession"].isin(collapse_map.keys())].copy()

        # For representative accessions, we add alt list using comma instead of semicolon
        df_pred["alternative_accessions"] = df_pred["accession"].map(
            lambda a: ",".join(alt_map.get(a, []))
        )

        # Load full original predictions to recompute merged names and means
        orig_pred = pd.read_csv(pred_path)

        split_pattern = r"[,;/]+"  # normalize splitting on anything previously used

        for rep, alts in alt_map.items():
            accs = [rep] + alts

            # Merge protein_names + gene_names from all members
            prot_names = []
            gene_names = []

            for a in accs:
                row_p = orig_pred[orig_pred["accession"] == a]
                if row_p.empty:
                    continue

                pn = row_p.iloc[0].get("protein_names", "")
                gn = row_p.iloc[0].get("gene_names", "")

                if isinstance(pn, str) and pn.strip():
                    prot_names.extend(
                        x.strip()
                        for x in re.split(split_pattern, pn)
                        if x.strip()
                    )

                if isinstance(gn, str) and gn.strip():
                    gene_names.extend(
                        x.strip()
                        for x in re.split(split_pattern, gn)
                        if x.strip()
                    )

            # Deduplicate but preserve order
            prot_names = list(dict.fromkeys(prot_names))
            gene_names = list(dict.fromkeys(gene_names))

            # Update representative rows using comma join instead of semicolon
            df_pred.loc[df_pred["accession"] == rep, "protein_names"] = ", ".join(prot_names)
            df_pred.loc[df_pred["accession"] == rep, "gene_names"] = ", ".join(gene_names)

            # --- mean prob_antigen ---
            mean_prob = orig_pred[orig_pred["accession"].isin(accs)]["prob_antigen"].mean()
            df_pred.loc[df_pred["accession"] == rep, "prob_antigen"] = mean_prob

            # --- mean allergenicity ---
            mean_allergen = df_raw[df_raw["accession"].isin(accs)]["allergenicity_hybrid_score"].mean()
            df_raw.loc[df_raw["accession"] == rep, "allergenicity_hybrid_score"] = mean_allergen

        # Remove FASTA sequences for collapsed accessions
        for alt in collapse_map.keys():
            if alt in seq_dict:
                del seq_dict[alt]


    # Determine human proteome FASTA
    if args.human_fasta:
        human_fasta = Path(args.human_fasta)
        if not human_fasta.exists():
            raise FileNotFoundError(f"Provided --human-fasta does not exist: {human_fasta}")
        logging.info(f"📂 Using provided human FASTA: {human_fasta}")
    else:
        cache_dir = Path("data/human_proteome")
        cache_dir.mkdir(parents=True, exist_ok=True)
        expected = cache_dir / "GCF_000001405.40_proteins.fasta"

        if expected.exists():
            human_fasta = expected
            logging.info(f"📂 Using cached human proteome: {human_fasta}")
        else:
            # Download
            logging.info("📥 Human proteome FASTA not found in cache; downloading now...")
            human_fasta_path = fetch_human_proteome("GCF_000001405.40", str(cache_dir))
            human_fasta = Path(human_fasta_path)

    # Build BLAST DB if necessary
    human_db = ensure_human_db(str(human_fasta))

    # BLASTP filtering
    logging.info("🔍 Running BLASTP similarity filtering against human proteome...")
    blast_flags = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for acc in df_pred["accession"]:
            seq = seq_dict.get(acc, None)
            if seq is None:
                blast_flags.append(False)  # cannot evaluate, keep
                continue

            if args.verbose:
                logging.info(f"Running BLASTP for accession: {acc}")
                temp_dir = Path("data/human_blast_tmp")
                temp_dir.mkdir(parents=True, exist_ok=True)
            else:
                temp_dir = tmpdir

            hit = run_blastp(acc, seq, human_db, temp_dir)
            if hit is None:
                blast_flags.append(False)  # no hit → keep
                continue

            # Filter: >30% identity AND E-value < 0.005 or E-value < 1e-6
            remove_flag = (hit["pident"] > 30.0) and (hit["evalue"] < 0.005) or (hit["evalue"] < 1e-6)
            blast_flags.append(remove_flag)

    df_pred["remove_human_similarity"] = blast_flags
    removed_human = df_pred[df_pred["remove_human_similarity"] == True]
    df_pred = df_pred[df_pred["remove_human_similarity"] == False]

    logging.info(f"Removed {len(removed_human)} for high similarity to human proteins.")
    
    # --- Allergenicity and MHC filtering with per-filter removal counts ---
    removed_allergenicity = df_raw[df_raw["allergenicity_hybrid_score"] >= 0.3]
    removed_mhci = df_raw[
        (df_raw["allergenicity_hybrid_score"] < 0.3) &
        (df_raw["mhci_num_strong_binders"] == 0)
    ]
    removed_mhcii = df_raw[
        (df_raw["allergenicity_hybrid_score"] < 0.3) &
        (df_raw["mhci_num_strong_binders"] > 0) &
        (df_raw["mhcii_num_strong_binders"] == 0)
    ]
 
    logging.info(f"🧪 Removed {len(removed_allergenicity)} proteins with allergenicity_hybrid_score >= 0.3")
    # save removed_allergenic rows to csv
    removed_allergenicity.to_csv(args.output_dir + "removed_allergenicity.csv", index=False)
    
    logging.info(f"🧬 Removed {len(removed_mhci)} proteins with mhci_num_strong_binders == 0 (after allergenicity filter)")
    logging.info(f"🧬 Removed {len(removed_mhcii)} proteins with mhcii_num_strong_binders == 0 (after allergenicity + mhci filters)")

    # --- Allergenicity and mhc filtering ---
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

    # include alternative accessions if present
    columns = ["accession", prob_col, "pred_label", "protein_names", "gene_names"]
    if "alternative_accessions" in merged.columns:
        columns.append("alternative_accessions")

    final_df = merged[columns]

    columns_removed = ["accession", prob_col, "pred_label", "protein_names", "gene_names"]
    if "alternative_accessions" in merged_removed.columns:
        columns_removed.append("alternative_accessions")

    removed_df = merged_removed[columns_removed]

    # Save
    final_df.to_csv(output_path, index=False)
    removed_output_path = output_path.parent / f"removed_{pred_path.stem}.csv"
    removed_df.to_csv(removed_output_path, index=False)

    # Save also human-similarity-removed list
    removed_human_path = output_path.parent / f"removed_human_similarity_{pred_path.stem}.csv"
    removed_human.to_csv(removed_human_path, index=False)

    logging.info(f"✅ Final filtered output saved to: {output_path}")
    logging.info(f"❌ Human similarity removed saved to: {removed_human_path}")


if __name__ == "__main__":
    main()
