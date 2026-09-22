# LDLR Variant Prioritization Pipeline

### High-Throughput Computational Reclassification of Variants of Uncertain Significance (VUS) in Familial Hypercholesterolemia

---

## Overview

Familial Hypercholesterolemia (FH) is a common, severe autosomal dominant genetic disorder caused primarily by mutations in the low-density lipoprotein receptor (*LDLR*) gene. Clinical sequencing frequently uncovers **Variants of Uncertain Significance (VUS)**, leaving patients without actionable diagnostic insights.

This repository houses a fully programmatic, reproducible command-line interface (CLI) pipeline designed to bypass manual web-based databases. By integrating genomic extraction, programmatic batch API annotation, multi-parameter clinical filtering, and 3D structural biophysics, this pipeline systematically reclassifies *LDLR* VUS into high-confidence pathogenic candidates.

---

## Pipeline Architecture

The workflow is managed via a fault-tolerant Bash orchestrator (`run_pipeline.sh`) executing five core programmatic modules:

1. **Genomic Data Extraction** — Isolates single-nucleotide VUS from raw ClinVar VCF files using `bcftools`.
2. **High-Throughput Batch Annotation** — Interfaces with the Ensembl REST Protein–Ensembl API (`/vep/human/region`) via Python (`requests`) with robust retry logic to extract evolutionary (SIFT), biophysical (PolyPhen-2), and population frequency (gnomAD) metrics.
3. **Strict Pathogenicity Filtering** — Applies mathematical cutoffs (`gnomAD < 0.0001`, `SIFT < 0.05`, `PolyPhen-2 > 0.90`) via `pandas` to isolate top-tier candidates.
4. **Functional Domain Mapping** — Extracts canonical transcript nomenclature (HGVS.p) and structural database classifications (InterPro/Pfam).
5. **Structural Validation** — Dynamically queries the AlphaFold API to retrieve the latest structural model (.cif) for PyMOL visualization.

---

## Repository Structure

```text
├── run_pipeline.sh          # Master Bash orchestrator script
├── annotate_gene.py         # Batch Ensembl REST API annotation module
├── rank_variants.py         # Multi-parameter clinical filtering module
├── map_domains.py           # Canonical transcript & HGVS domain mapping module
├── download_structure.py    # Dynamic AlphaFold CIF fetcher module
├── top_LDLR_candidates.csv  # Mathematically prioritized variant dataset
├── figurr1_ldlr_final.png   # 3D PyMOL structural visualization figure
└── README.md                # Project documentation
```

---

## Prerequisites & Dependencies

Ensure your environment (WSL Ubuntu or Linux) satisfies the following requirements:

* **Core Tools:** `bash`, `git`, `bcftools`
* **Python Environment (v3.10+):** `pandas`, `requests`, `argparse`
* **Visualization:** PyMOL (v3.x)

---

## Step-by-Step Execution Guide

### 1. Clone the Repository

```bash
git clone https://github.com/<YOUR_USERNAME>/LDLR-Variant-Prioritization.git
cd LDLR-Variant-Prioritization
```

### 2. Prepare the Input Data

Download the latest human ClinVar VCF and its index file into your root working directory:

```bash
# Example paths for ClinVar GRCh38 VCF
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

### 3. Run the Master Orchestrator

Execute the complete end-to-end workflow with a single command:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh clinvar.vcf.gz LDLR 19 P01130
```

---

## Results & Biological Interpretation

Running the pipeline on the *LDLR* gene (Chr 19, UniProt: P01130) yielded the following structural and evolutionary insights:

* **Data Funnel** — Out of **1,319** raw ClinVar SNV VUS targets, **242** high-confidence pathogenic candidates were isolated following strict frequency and score thresholds.
* **Structural Impact** — Top variants (such as `p.Cys46Trp`, `p.Asp57His`, `p.Cys63Phe/Tyr`, and `p.Val295Gly`) cluster aggressively inside critical functional zones.
* **Disulfide Network Disruption** — Mutations like `Cys63` systematically destroy the rigid cysteine residues required to maintain the disulfide bond scaffolding of the ligand-binding domain, explaining the structural mechanism of pathogenicity.

### 3D Structural Validation (PyMOL)

To render the publication-grade structural model displaying the disrupted cysteine–disulfide bonding network:

```python
load LDLR_structure.cif, ldlr
bg_color white
hide everything
show cartoon
color gray70, ldlr
select pathogenic_vus, resi 46+57+63+295
show sticks, pathogenic_vus
color red, pathogenic_vus
select disulfides, resname CYS and resi 40-80
show sticks, disulfides
color yellow, disulfides
color red, pathogenic_vus
zoom pathogenic_vus, 15
ray 1600, 1200
png figurr1_ldlr_final.png, dpi=300
```

---

## License & Author

Developed as an open-source, reproducible bioinformatics pipeline suitable for academic research and clinical genomics variant interpretation portfolios.

* **Author:** [Mobeen Nasir / mubeennasir117@gmail.com]
* **Target Condition:** Familial Hypercholesterolemia (CDC Tier-1 Genomic Application)
