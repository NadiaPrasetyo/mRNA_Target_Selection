"""
IEDB_fetch.py
Command-line tool to fetch antigen and epitope data for a given organism from the IEDB API.

Overview:
    - Queries the IEDB (Immune Epitope Database) API for antigen and epitope records associated with a specified source organism.
    - Saves the retrieved data as structured CSV files in a user-specified output directory.

Arguments:
    output_folder (str): Subdirectory under `data/` where results will be saved.
    source_organism (str): Scientific name of the target organism to query (e.g., "SARS-CoV-2").

Requirements:
    - Python packages: requests, csv, argparse, os.

Usage Example:
    python IEDB_fetch.py sars_cov_2 "SARS-CoV-2"

Outputs:
    data/<output_folder>/<organism_tag>_IEDB_antigens.csv   # Antigen records for the organism
    data/<output_folder>/<organism_tag>_IEDB_epitope.csv    # Epitope records for the organism

Author: Nadia
"""

import requests
import csv
import os
import argparse

# IEDB API endpoints
antigen_url = "https://query-api.iedb.org/antigen_search"
epitope_url = "https://query-api.iedb.org/epitope_search"

def fetch_and_save(url, params, out_path):
    """
    Fetch data from IEDB API and save to CSV file. Queries the specified URL with given parameters,
    and writes the results to a CSV file at the specified output path. Handles errors and logs progress.
    Args:
        url (str): The IEDB API endpoint to query.
        params (dict): Dictionary of query parameters.
        out_path (str): Path to the output CSV file.
    Returns:
        None
    """
    headers = {
        "Accept": "application/json",
        "Prefer": "count=exact"
    }
    try:
        response = requests.get(url, params=params, headers=headers)
    except requests.exceptions.RequestException as e:
        print(f"❌ Error during request: {e}")
        return

    if response.status_code == 200:
        data = response.json()
        if data:
            fieldnames = data[0].keys()
            with open(out_path, mode="w", newline='', encoding='utf-8') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(data)
            print(f"✅ Data saved to '{out_path}' with {len(data)} records.")
        else:
            print(f"⚠️ No data returned for {out_path}.")
        content_range = response.headers.get("Content-Range")
        if content_range:
            total_count = content_range.split("/")[-1]
            print(f"🔢 Total matching records: {total_count}")
    else:
        print(f"❌ Failed to fetch data from {url}. Status code: {response.status_code}")
        print("Response content:", response.text)

def main():
    """Main function to parse command-line arguments and fetch IEDB data."""
    # Argument parser setup
    parser = argparse.ArgumentParser(
        description="Fetches antigen and epitope data for a given organism from the IEDB API.",
        usage="python fetch_iedb_data.py <output_folder> <source_organism>"
    )
    parser.add_argument("output_folder", help="Directory name under data/")
    parser.add_argument("source_organism", help='Full organism name (e.g., "SARS-CoV-2")')
    args = parser.parse_args()

    output_folder = f"data/{args.output_folder}"
    source_organism = args.source_organism.lower()
    organism_tag = source_organism.replace(" ", "_")

    os.makedirs(output_folder, exist_ok=True)

    # Antigen data fetch
    antigen_params = {
        "host_organism_name": "ilike.*human*",
        "source_organism_name": f"ilike.*{source_organism}*",
        "select": "parent_source_antigen_iri,parent_source_antigen_names,source_organism_names,host_organism_names",
    }
    antigen_out = os.path.join(output_folder, f"{organism_tag}_IEDB_antigens.csv")
    fetch_and_save(antigen_url, antigen_params, antigen_out)

    # Epitope data fetch
    epitope_params = {
        "host_organism_name": "ilike.*human*",
        "source_organism_name": f"ilike.*{source_organism}*",
        "select": "structure_iri,pdb_ids,qualitative_measures,linear_sequence,linear_sequence_length,parent_source_antigen_iris,parent_source_antigen_names,source_organism_names,host_organism_names",
    }
    epitope_out = os.path.join(output_folder, f"{organism_tag}_IEDB_epitope.csv")
    fetch_and_save(epitope_url, epitope_params, epitope_out)

if __name__ == "__main__":
    main()
