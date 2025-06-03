"""
/**
 * @file compile_antigens.py
 * @brief Compiles antigen data from IEDB and literature sources for a given organism.
 *
 * This script aggregates antigen information from IEDB CSV files and literature/patent Excel files
 * for a specified organism. It standardizes the data format and outputs a single compiled CSV file.
 *
 * General Flow:
 *   1. Loads IEDB antigen data from a CSV file.
 *   2. Loads literature/patent antigen data from Excel files.
 *   3. Standardizes and merges the data from all sources.
 *   4. Outputs a compiled CSV file containing all antigen information.
 *
 * Parameters:
 *   short_name (str): Short identifier for the organism (used as a directory name).
 *   long_name (str): Full organism name (used for file naming and metadata).
 *
 * Usage:
 *   python compile_antigens.py <short_name> <long_name>
 *
 * Example:
 *   python compile_antigens.py sars_cov_2 "SARS-CoV-2"
 *
 * @author Nadia
 */
"""

import pandas as pd
import os
from glob import glob
import re
import argparse

"""
/**
* @brief Extracts the UniProt ID from an IRI string.
*
* This function parses a given IRI (Internationalized Resource Identifier) and extracts
* the UniProt accession if present.
*
* @param iri (str): The IRI string potentially containing a UniProt ID.
* @return (str or None): The extracted UniProt ID, or None if not found or input is NaN.
*/
"""
def extract_uniprot_id(iri):
    if pd.isna(iri):
        return None
    match = re.search(r"(?:UNIPROT:)?([A-Z0-9]+)$", iri)
    if match:
        return f"{match.group(1)}"
    return None

"""
/**
* @brief Loads and standardizes antigen data from an IEDB CSV file.
*
* Reads a CSV file containing IEDB antigen data, renames columns to a standard format,
* extracts UniProt IDs, and adds metadata columns.
*
* @param file_path (str): Path to the IEDB antigen CSV file.
* @return (pd.DataFrame): Standardized DataFrame with antigen information.
*/
"""
def load_iedb_antigens(file_path):
    if not os.path.exists(file_path):
        print(f"Warning: IEDB file not found: {file_path}. Skipping.")
        return pd.DataFrame()
    df = pd.read_csv(file_path)

    df_out = df.rename(columns={
        "source_organism_names": "source_organism",
        "host_organism_names": "host_organisms",
        "parent_source_antigen_names": "antigen_name"
    })

    # Extract and format Uniprot_ID
    if "parent_source_antigen_iri" in df.columns:
        df_out["Uniprot_ID"] = df["parent_source_antigen_iri"].apply(extract_uniprot_id)
    else:
        df_out["Uniprot_ID"] = None
    df_out["gene_name"] = None  # IEDB does not provide gene names
    df_out["source"] = "IEDB"
    return df_out[["source_organism", "host_organisms", "antigen_name", "gene_name", "Uniprot_ID", "source"]]

"""
/**
* @brief Loads and standardizes antigen data from a literature/patent Excel file.
*
* Reads an Excel file containing antigen data from literature or patents, standardizes
* the columns, and adds metadata.
*
* @param file_path (str): Path to the literature/patent Excel file.
* @param source_organism (str): Name of the source organism.
* @return (pd.DataFrame): Standardized DataFrame with antigen information.
*/
"""
def load_literature_antigens(file_path, source_organism):
    if not os.path.exists(file_path):
        print(f"Warning: Literature file not found: {file_path}. Skipping.")
        return pd.DataFrame()
    df = pd.read_excel(file_path)
    df_out = pd.DataFrame()
    df_out["antigen_name"] = df["Name"]
    df_out["gene_name"] = df["Gene"]
    df_out["source_organism"] = source_organism
    df_out["host_organisms"] = "Homo sapiens"
    df_out["Uniprot_ID"] = None
    df_out["source"] = "literature"
    return df_out[["source_organism", "host_organisms", "antigen_name", "gene_name", "Uniprot_ID", "source"]]

"""
/**
* @brief Main pipeline for compiling antigen data from all sources for a given organism.
*
* Loads IEDB and literature antigen data, standardizes and merges them, and writes the
* compiled data to a CSV file.
*
* @param short_name (str): Short identifier for the organism.
* @param long_name (str): Full organism name.
* @return: None
*/
    """
def main(short_name, long_name):
    organism_tag = str(long_name).replace(" ", "_").lower()
    base_path = f"data/{short_name}"

    # Load IEDB
    iedb_file = os.path.join(base_path, f"{organism_tag}_IEDB_antigens.csv")
    iedb_df = load_iedb_antigens(iedb_file)

    # Load literature/patents (XLSX files)
    literature_dfs = []
    for xlsx_file in glob(os.path.join(base_path, "*.xlsx")):
        lit_df = load_literature_antigens(xlsx_file, long_name)
        if not lit_df.empty:
            literature_dfs.append(lit_df)
    literature_df = pd.concat(literature_dfs, ignore_index=True) if literature_dfs else pd.DataFrame()

    # Combine all sources
    combined_df = pd.concat([iedb_df, literature_df], ignore_index=True)
    if combined_df.empty:
        print("No antigen data found. Exiting.")
        return

    # Output
    output_file = os.path.join(base_path, f"{organism_tag}_compiled_antigens.csv")
    combined_df.to_csv(output_file, index=False)
    print(f"Compiled antigen data saved to: {output_file}")

"""
/**
* @brief Entry point for the script. Parses arguments and runs the main pipeline.
*
* Expects two command-line arguments: short_name and long_name.
*/
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compile antigen data from IEDB and literature sources for a given organism.",
        usage="python compile_antigens.py <pathogen_directory> <pathogen_name>"
    )
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("pathogen_name", help='Prefix used in filenames (e.g., "staphylococcus aureus")')
    args = parser.parse_args()

    if not os.path.exists(f"data/{args.pathogen_directory}"):
        os.makedirs(f"data/{args.pathogen_directory}")
        print(f"Created directory: data/{args.pathogen_directory}")
    main(args.pathogen_directory, args.pathogen_name)
