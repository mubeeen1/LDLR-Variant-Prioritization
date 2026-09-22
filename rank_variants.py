import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(description='Strict Pathogenicity Filtering')
parser.add_argument('--input', required=True, help='Raw annotated CSV')
parser.add_argument('--output', required=True, help='Top candidates CSV')
args = parser.parse_args()

print(f"Filtering {args.input} for pathogenic thresholds...")
try:
    df = pd.read_csv(args.input)
except FileNotFoundError:
    print(f"Error: Could not find {args.input}")
    sys.exit(1)

df['gnomAD_Freq'] = pd.to_numeric(df['gnomAD_Freq'], errors='coerce').fillna(0) 
df['SIFT_Score'] = pd.to_numeric(df['SIFT_Score'], errors='coerce')
df['PolyPhen_Score'] = pd.to_numeric(df['PolyPhen_Score'], errors='coerce')

rare_mask = df['gnomAD_Freq'] < 0.0001 
sift_mask = df['SIFT_Score'] < 0.05 
polyphen_mask = df['PolyPhen_Score'] > 0.90 

top_candidates = df[rare_mask & sift_mask & polyphen_mask].copy()
top_candidates = top_candidates.sort_values(by='PolyPhen_Score', ascending=False)
top_candidates.to_csv(args.output, index=False)

print(f"Isolated {len(top_candidates)} highly pathogenic variants. Saved to {args.output}")
