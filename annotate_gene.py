import requests
import json
import pandas as pd
import time
import argparse
import sys

parser = argparse.ArgumentParser(description='Batch API Annotation for VUS')
parser.add_argument('--input', required=True, help='Input target file (e.g., targets.txt)')
parser.add_argument('--output', required=True, help='Output CSV file (e.g., annotated_vus.csv)')
args = parser.parse_args()

print(f"Loading {args.input}...")
try:
    variants = pd.read_csv(args.input, sep='\t', header=None, names=['CHROM', 'POS', 'REF', 'ALT'], dtype=str)
except FileNotFoundError:
    print(f"Error: Could not find {args.input}")
    sys.exit(1)

# Filter for SNVs
variants = variants[(variants['REF'].str.len() == 1) & (variants['ALT'].str.len() == 1)]
total = len(variants)
print(f"Total valid SNVs: {total}. Initiating Ensembl REST API bulk requests...\n")

server = "https://rest.ensembl.org"
ext = "/vep/human/region"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

query_map = {}
variant_queries = []

for _, row in variants.iterrows():
    query_str = f"{row['CHROM']}:{row['POS']}-{row['POS']}:1/{row['ALT']}"
    variant_queries.append(query_str)
    query_map[query_str] = {
        'Chromosome': row['CHROM'], 'Position': row['POS'], 'REF': row['REF'], 'ALT': row['ALT']
    }

results = []
batch_size = 50 

for i in range(0, len(variant_queries), batch_size):
    batch = variant_queries[i : i + batch_size]
    payload = {"variants": batch}
    
    # Robust Retry Logic for API Timeouts
    max_retries = 3
    response = None
    
    for attempt in range(max_retries):
        try:
            response = requests.post(server+ext, headers=headers, data=json.dumps(payload), timeout=45)
            if response.ok:
                break
            else:
                print(f"[Warning] Batch {i//batch_size + 1} hiccuped (HTTP {response.status_code}). Retrying ({attempt+1}/{max_retries})...")
                time.sleep(5)
        except requests.exceptions.RequestException as e:
            print(f"[Warning] Batch {i//batch_size + 1} network timeout. Retrying ({attempt+1}/{max_retries})...")
            time.sleep(5)
            
    if response is None or not response.ok:
        print(f"[Error] Batch {i//batch_size + 1} failed after {max_retries} attempts. Skipping.")
        continue
        
    try:
        decoded = response.json()
        
        for var_data in decoded:
            input_str = var_data.get('input', '')
            orig = query_map.get(input_str)
            if not orig: continue
            
            transcript = var_data.get('transcript_consequences', [{}])[0] if var_data.get('transcript_consequences') else {}
            sift_score = transcript.get('sift_score', 'N/A')
            polyphen_score = transcript.get('polyphen_score', 'N/A')
            
            gnomad_af = 'N/A'
            if 'colocated_variants' in var_data:
                for cv in var_data['colocated_variants']:
                    if 'frequencies' in cv and orig['ALT'] in cv['frequencies']:
                        gnomad_af = cv['frequencies'][orig['ALT']].get('gnomad', 'N/A')
                        break
                        
            results.append({
                'Chromosome': orig['Chromosome'], 'Position': orig['Position'],
                'Mutation': f"{orig['REF']}>{orig['ALT']}",
                'SIFT_Score': sift_score, 'PolyPhen_Score': polyphen_score, 'gnomAD_Freq': gnomad_af
            })
            
        print(f"--> Successfully processed batch {i//batch_size + 1}")
        time.sleep(1) 
        
    except Exception as e:
         print(f"[Error] Batch {i//batch_size + 1} JSON parsing failed: {e}")

df = pd.DataFrame(results)
df.to_csv(args.output, index=False)
print(f"\nSuccess! Saved annotations to {args.output}")
