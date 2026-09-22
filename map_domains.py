import requests
import pandas as pd
import argparse
import sys

# 1. COMMAND LINE ARGUMENTS (No more hardcoding)
parser = argparse.ArgumentParser(description='Map genomic variants to protein domains.')
parser.add_argument('--input', required=True, help='Input CSV file (e.g., top_ldlr_candidates.csv)')
parser.add_argument('--chrom', required=True, help='Chromosome number to filter (e.g., 19)')
args = parser.parse_args()

print(f"Loading {args.input} and isolating Chromosome {args.chrom}...")
try:
    df = pd.read_csv(args.input)
except FileNotFoundError:
    print(f"Error: Could not find {args.input}")
    sys.exit(1)

# Filter for the target chromosome and take the top 5 candidates
df['Chromosome'] = df['Chromosome'].astype(str)
true_targets = df[df['Chromosome'] == str(args.chrom)].head(5)

server = "https://rest.ensembl.org"
headers = {"Content-Type": "application/json", "Accept": "application/json"}
results = []

print("Querying Ensembl API for Canonical Transcripts, HGVS nomenclature, and Domains...\n")

for index, row in true_targets.iterrows():
    ref, alt = row['Mutation'].split('>')
    
    # Added hgvs=1 to pull exact nomenclature
    ext = f"/vep/human/region/{row['Chromosome']}:{row['Position']}-{row['Position']}:1/{alt}?Domains=1&hgvs=1"
    
    try:
        response = requests.get(server+ext, headers=headers)
        if not response.ok:
            continue
            
        decoded = response.json()
        consequences = decoded[0].get('transcript_consequences', [])
        
        # 2. ISOLATE THE CANONICAL PROTEIN TRANSCRIPT
        best_transcript = {}
        for t in consequences:
            if 'hgvsp' in t and t.get('biotype') == 'protein_coding':
                best_transcript = t
                if t.get('canonical') == 1:
                    break # Perfect match found
        
        hgvsp = best_transcript.get('hgvsp', 'Unknown').split(':')[-1] 
        amino_acids = best_transcript.get('amino_acids', 'Unknown')
        
        # 3. EXTRACT STRUCTURAL DOMAINS
        domains = []
        if 'domains' in best_transcript:
            for d in best_transcript['domains']:
                # Filter for recognized biological structural databases
                if d.get('db') in ['InterPro', 'Pfam', 'PROSITE']:
                    domains.append(f"{d.get('db')}:{d.get('name')}")
        
        # Remove duplicates and format cleanly
        domain_str = " | ".join(list(dict.fromkeys(domains))[:2]) if domains else "No defined domain"
        
        results.append({
            'Genomic_Pos': row['Position'],
            'Protein_Change': hgvsp,
            'Amino_Acids': amino_acids,
            'Domain_Location': domain_str
        })
        
    except Exception as e:
        print(f"Error on {row['Position']}: {e}")

final_df = pd.DataFrame(results)
print("--- FINAL BIOLOGICAL INTERPRETATION TABLE ---")
print(final_df.to_string(index=False))
