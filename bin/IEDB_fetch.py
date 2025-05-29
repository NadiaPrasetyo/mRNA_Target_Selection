import requests
import csv
import sys
import os

antigen_url = "https://query-api.iedb.org/antigen_search"
epitope_url = "https://query-api.iedb.org/epitope_search"

if len(sys.argv) > 2:
    output_folder = f'data/{sys.argv[1]}'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"📁 Output folder '{output_folder}' will be created.")

    source_organism = sys.argv[2].lower()
    organism_tag = source_organism.replace(" ", "_")
else:
    print("❌ Please provide the output folder and source organism name as arguments.\n"
          "Usage: python script.py <output_folder> <source_organism>")
    exit(1)

os.makedirs(output_folder, exist_ok=True)

def fetch_and_save(url, params, out_path):
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

# Antigen fetch
antigen_params = {
    "host_organism_name": "ilike.*human*",
    "source_organism_name": f"ilike.*{source_organism}*",
    "select": "parent_source_antigen_iri,parent_source_antigen_names,source_organism_names,host_organism_names",
}
antigen_out = os.path.join(output_folder, f"{organism_tag}_IEDB_antigens.csv")
fetch_and_save(antigen_url, antigen_params, antigen_out)

# Epitope fetch
epitope_params = {
    "host_organism_name": "ilike.*human*",
    "source_organism_name": f"ilike.*{source_organism}*",
    "select": "structure_iri,pdb_ids,qualitative_measures,linear_sequence,linear_sequence_length,parent_source_antigen_iris,parent_source_antigen_names,source_organism_names,host_organism_names",
}
epitope_out = os.path.join(output_folder, f"{organism_tag}_IEDB_epitope.csv")
fetch_and_save(epitope_url, epitope_params, epitope_out)
