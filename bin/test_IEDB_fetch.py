import requests
import csv

# API URL with proper query syntax
url = "https://query-api.iedb.org/antigen_search"
params = {
    "host_organism_name": "ilike.*human*", 
    "source_organism_name": "ilike.*staphylococcus aureus*",
    "select": "parent_source_antigen_iri,parent_source_antigen_names,source_organism_names,host_organism_names",
}

# Headers with Accept and Prefer
headers = {
    "Accept": "application/json",   # Requesting JSON
    "Prefer": "count=exact" # Requesting exact count of the total results
}

try:
    # Make the request
    response = requests.get(url, params=params, headers=headers)
except requests.exceptions.RequestException as e:
    print(f"❌ An error occurred while making the request: {e}")
    exit(1)

# Check for success
if response.status_code == 200:
    data = response.json()

    # Save to CSV
    output_file = "data/S.aureus/IEDB_antigens.csv"
    if data:
        fieldnames = data[0].keys()
        with open(output_file, mode="w", newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
        print(f"✅ Data saved to '{output_file}' with {len(data)} records.")
    else:
        print("⚠️ No data returned by the API.")

    # Total count from header if available
    content_range = response.headers.get("Content-Range")
    if content_range:
        total_count = content_range.split("/")[-1]
        print(f"🔢 Total matching records: {total_count}")

else:
    print(f"❌ Failed to fetch data. Status code: {response.status_code}")
    print("Response content:", response.text)
