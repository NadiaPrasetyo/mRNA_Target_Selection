#!/usr/bin/env python3
import argparse
import requests
import pandas as pd
import random
from pathlib import Path
from time import sleep
import sys

SEARCH_URL = "https://rest.uniprot.org/uniparc/search"
UNIPARC_URL = "https://rest.uniprot.org/uniparc"
DELAY = 0.3  # seconds between API calls


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Fetch one active UniProtKB accession per input accession from UniParc API. "
            "If multiple exist, prefer the one matching the given organism scientific name."
        )
    )
    parser.add_argument(
        "-i", "--input_dir", required=True,
        help="Directory containing input CSVs with 'accession' column."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output file (.csv) or directory."
    )
    parser.add_argument(
        "-org", "--organism", required=True,
        help="Target organism scientific name (exact match)."
    )
    return parser.parse_args()


def resolve_output_path(output_arg: str) -> Path:
    """Create directories if needed and return output file path."""
    out_path = Path(output_arg)
    if out_path.suffix.lower() != ".csv":
        out_path.mkdir(parents=True, exist_ok=True)
        out_path = out_path / "uniprot_accessions.csv"
    else:
        out_path.parent.mkdir(parents=True, exist_ok=True)
    return out_path


def load_unique_accessions(input_dir: str):
    """Read all CSVs and collect unique accessions."""
    all_accessions = set()
    for file in Path(input_dir).glob("*.csv"):
        try:
            df = pd.read_csv(file, usecols=["accession"])
            accs = df["accession"].dropna().unique()
            all_accessions.update(accs)
            print(f"📄 {file.name}: {len(accs)} accessions")
        except ValueError:
            print(f"⚠️ Skipping {file.name}: no 'accession' column found")
    print(f"\n✅ Found {len(all_accessions)} unique accessions total.\n")
    return sorted(all_accessions)


def get_uniparc_upi(accession):
    """Look up UPI for a given accession using UniParc search endpoint."""
    headers = {"accept": "application/json"}
    params = {"query": accession, "fields": "upi", "size": 1}

    response = requests.get(SEARCH_URL, headers=headers, params=params, timeout=30)
    if not response.ok:
        print(f"⚠️ Failed to get UPI for {accession}: HTTP {response.status_code}")
        return None

    data = response.json()
    results = data.get("results", [])
    if not results:
        print(f"⚠️ No UPI found for {accession}")
        return None

    return results[0].get("uniParcId")


def get_best_active_uniprot_accession(upi, query_accession, organism_name):
    """Fetch all active UniProtKB accessions and select one based on organism name."""
    headers = {"accept": "application/json"}
    url = f"{UNIPARC_URL}/{upi}"
    response = requests.get(url, headers=headers, timeout=30)
    if not response.ok:
        print(f"⚠️ API error for UPI {upi}: HTTP {response.status_code}")
        return None

    data = response.json()
    xrefs = data.get("uniParcCrossReferences", [])
    active_xrefs = [
        x for x in xrefs
        if x.get("database") in {"UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL"} and x.get("active") is True
    ]

    if not active_xrefs:
        return None

    # Prefer exact organism match if available
    matched = [
        x for x in active_xrefs
        if x.get("organism", {}).get("scientificName", "").strip().lower() == organism_name.strip().lower()
    ]

    chosen = matched[0] if matched else random.choice(active_xrefs)
    acc_id = chosen.get("id", "").split(".")[0]
    org_name = chosen.get("organism", {}).get("scientificName", "N/A")

    print(f"   ✅ Selected {acc_id} (organism: {org_name})")

    return {
        "query_accession": query_accession,
        "uniparc_upi": upi,
        "Uniprot_ID": acc_id,
        "organism": org_name
    }


def main():
    args = parse_args()
    output_path = resolve_output_path(args.output)
    organism_name = args.organism
    accessions = load_unique_accessions(args.input_dir)

    all_results = []

    for i, acc in enumerate(accessions, 1):
        print(f"\n🔎 [{i}/{len(accessions)}] Processing {acc} ...")

        upi = get_uniparc_upi(acc)
        if not upi:
            continue

        best_acc = get_best_active_uniprot_accession(upi, acc, organism_name)
        if best_acc:
            all_results.append(best_acc)

        pd.DataFrame(all_results).to_csv(output_path, index=False)
        sleep(DELAY)

    print(f"\n✅ Finished! Saved {len(all_results)} accessions to {output_path}")


if __name__ == "__main__":
    main()

