# This script fetches, displays and saves all supported organisms from the PANTHER database.

import requests
import json
import os

URL = "https://pantherdb.org/services/oai/pantherdb/supportedgenomes"

try:
    # Send GET request to the PANTHER API
    response = requests.get(URL, timeout=30)
    response.raise_for_status()  # raise an error for HTTP problems

    # Parse the JSON response
    data = response.json()
    genomes = data["search"]["output"]["genomes"]["genome"]

    # Display formatted results in console
    print(f"Found {len(genomes)} supported genomes:\n")
    for g in genomes:
        name = g.get("name", "N/A")
        long_name = g.get("long_name", "N/A")
        taxon_id = g.get("taxon_id", "N/A")
        version = g.get("version", "N/A")

        print(f"{name:20s} | {long_name:35s} | id = NCBITaxon:{taxon_id:<8} | {version}")

    # --- Restructure into a simplified JSON array ---
    json_array = [
        {
            "name": g.get("name", "N/A"),
            "long_name": g.get("long_name", "N/A"),
            "taxon_id": f"NCBITaxon:{g.get('taxon_id', 'N/A')}",
            "version": g.get("version", "N/A")
        }
        for g in genomes
    ]

    # Ensure output directory exists
    # os.makedirs("data_files", exist_ok=True)
    # output_path = os.path.join("data_files", "supported_genomes.json")
    output_path = "supported_organisms.json"

    # Write to JSON file
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(json_array, f, ensure_ascii=False, indent=4)

    print(f"\n✅ Saved structured JSON to: {output_path}")

except requests.exceptions.RequestException as e:
    print(f"❌ Network or HTTP error: {e}")
except (KeyError, ValueError, json.JSONDecodeError) as e:
    print(f"❌ JSON structure error: {e}")