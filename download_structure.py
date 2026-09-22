import requests
import urllib.request
import argparse
import sys

# 1. COMMAND LINE ARGUMENTS
parser = argparse.ArgumentParser(description='Dynamically download the latest AlphaFold CIF structure.')
parser.add_argument('--uniprot', required=True, help='UniProt Accession ID (e.g., P01130)')
parser.add_argument('--output', default=None, help='Output filename (optional)')
args = parser.parse_args()

# Default the filename to the UniProt ID if not provided
out_file = args.output if args.output else f"{args.uniprot}_structure.cif"

print(f"Querying AlphaFold API for the current {args.uniprot} structure...")
url = f"https://alphafold.ebi.ac.uk/api/prediction/{args.uniprot}"

try:
    response = requests.get(url, timeout=10)
    
    if response.ok:
        data = response.json()
        if not data:
            print(f"No AlphaFold structure found for {args.uniprot}.")
            sys.exit(1)
            
        latest_cif = data[0]['cifUrl'] 
        print(f"Resolved active URL: {latest_cif}")
        
        print(f"Downloading to {out_file}...")
        urllib.request.urlretrieve(latest_cif, out_file)
        print(f"Success! Saved as {out_file}")
    else:
        print(f"API Error: {response.status_code} - {response.reason}")
        
except Exception as e:
    print(f"Connection failed: {e}")
