import pandas as pd
import requests
import json
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Top candidates CSV')
parser.add_argument('--uniprot', required=True, help='UniProt ID for domain boundaries')
args = parser.parse_args()

df = pd.read_csv(args.input)

# 1. Fetch gold-standard domain boundaries directly from UniProt
uniprot_url = f"https://rest.uniprot.org/uniprotkb/{args.uniprot}.json"
response = requests.get(uniprot_url)
domains = []

if response.ok:
    features = response.json().get('features', [])
    for f in features:
        if f.get('type') in ['Domain', 'Repeat', 'Region', 'Topological domain']:
            start = f.get('location', {}).get('start', {}).get('value')
            end = f.get('location', {}).get('end', {}).get('value')
            desc = f.get('description', f.get('type'))
            if start and end:
                domains.append({'start': int(start), 'end': int(end), 'desc': desc})
else:
    print("Failed to reach UniProt API.")
    exit(1)

# 2. Query Ensembl ONLY to translate Genomic POS to Protein POS
server = "https://rest.ensembl.org"
ext = "/vep/human/region"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

print(f"{'Genomic_Pos':<15} {'Protein_Change':<18} {'Biological_Domain'}")
print("-" * 75)

# Take the top 5 candidates for the report
top_5 = df.head(5)
queries = [f"{row['Chromosome']}:{row['Position']}-{row['Position']}:1/{row['Mutation'].split('>')[1]}" for _, row in top_5.iterrows()]
payload = {"variants": queries}

try:
    res = requests.post(server+ext, headers=headers, data=json.dumps(payload), timeout=15)
    decoded = res.json()
    
    for var_data in decoded:
        input_str = var_data.get('input', '')
        genomic_pos = input_str.split(':')[1].split('-')[0]
        
        transcript = var_data.get('transcript_consequences', [{}])[0]
        protein_start = transcript.get('protein_start')
        
        # Smart HGVS extraction/fallback
        hgvsp = transcript.get('hgvsp', '')
        if hgvsp:
            prot_change = hgvsp.split(':')[-1]
        else:
            aa = transcript.get('amino_acids', '')
            if aa and protein_start and '/' in aa:
                ref_aa, alt_aa = aa.split('/')
                prot_change = f"p.{ref_aa}{protein_start}{alt_aa}"
            else:
                prot_change = "Unknown"
        
        domain_name = "Intergenic / No Domain"
        if protein_start:
            p_loc = int(protein_start)
            # Check which UniProt boundary this amino acid falls inside
            matches = [d['desc'] for d in domains if d['start'] <= p_loc <= d['end']]
            if matches:
                domain_name = " | ".join(matches)
                
        print(f"{genomic_pos:<15} {prot_change:<18} {domain_name}")
except Exception as e:
    print(f"Error mapping variants: {e}")
