#!/bin/bash

# Ensure correct usage
if [ "$#" -ne 4 ]; then
    echo "Usage: ./run_pipeline.sh <clinvar.vcf.gz> <GENE_NAME> <CHROMOSOME> <UNIPROT_ID>"
    echo "Example: ./run_pipeline.sh clinvar.vcf.gz LDLR 19 P01130"
    exit 1
fi

VCF=$1
GENE=$2
CHROM=$3
UNIPROT=$4

echo "=========================================================="
echo " Starting End-to-End VUS Pipeline for $GENE (Chr $CHROM)"
echo "=========================================================="

# STEP 1: Genomic Extraction
echo -e "\n[1/5] Extracting VUS from ClinVar..."
bcftools view $VCF | bcftools filter -i "INFO/GENEINFO ~ '(^|\|)${GENE}:' && INFO/CLNSIG='Uncertain_significance'" > ${GENE}_vus_only.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${GENE}_vus_only.vcf > ${GENE}_targets.txt
echo "Extracted strictly anchored raw targets to ${GENE}_targets.txt"

echo -e "\n[2/5] Fetching API Annotations (SIFT, PolyPhen, gnomAD)..."
python annotate_gene.py --input ${GENE}_targets.txt --output ${GENE}_annotated_vus.csv

echo -e "\n[3/5] Applying Clinical Filters..."
python rank_variants.py --input ${GENE}_annotated_vus.csv --output top_${GENE}_candidates.csv

echo -e "\n[4/5] Mapping Biological Domains..."
python map_domains.py --input top_${GENE}_candidates.csv --uniprot $UNIPROT > ${GENE}_final_report.txt

echo -e "\n[5/5] Downloading AlphaFold Structure..."
python download_structure.py --uniprot $UNIPROT --output ${GENE}_structure.cif

echo -e "\n=========================================================="
echo " PIPELINE COMPLETE."
echo " 1. Review biological context in: ${GENE}_final_report.txt"
echo " 2. Load ${GENE}_structure.cif into PyMOL for visualization."
echo "=========================================================="
