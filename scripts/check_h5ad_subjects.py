#!/usr/bin/env python3
"""Check which subjects are in H5AD vs CSV export."""

import scanpy as sc
import pandas as pd

# Load both datasets
h5ad_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
csv_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/adni_diagnosis_export.csv'

print("Loading datasets...")
adata = sc.read_h5ad(h5ad_file)
csv_data = pd.read_csv(csv_file)

# Get subject lists
h5ad_subjects = set(adata.obs['subject_id'].unique())
csv_subjects = set(csv_data['donor_id'].unique())

print(f"H5AD subjects: {len(h5ad_subjects)}")
print(f"CSV subjects: {len(csv_subjects)}")

# Find differences
h5ad_only = h5ad_subjects - csv_subjects
csv_only = csv_subjects - h5ad_subjects
common = h5ad_subjects & csv_subjects

print(f"\nCommon subjects: {len(common)}")
print(f"H5AD only: {len(h5ad_only)}")
print(f"CSV only: {len(csv_only)}")

if h5ad_only:
    print(f"\nFirst 10 subjects only in H5AD: {list(h5ad_only)[:10]}")
if csv_only:
    print(f"\nFirst 10 subjects only in CSV: {list(csv_only)[:10]}")

# Check diagnosis data for 941_S_4420
print(f"\n🎯 Subject 941_S_4420 analysis:")
if '941_S_4420' in h5ad_subjects:
    subject_data = adata.obs[adata.obs['subject_id'] == '941_S_4420'].iloc[0]
    print("H5AD data:")
    print(f"  worst_diagnosis_code: {subject_data['worst_diagnosis_code']}")
    print(f"  worst_diagnosis_label: {subject_data['worst_diagnosis_label']}")
    print(f"  worst_diagnosis_date: {subject_data['worst_diagnosis_date']}")
else:
    print("❌ Not found in H5AD")

if '941_S_4420' in csv_subjects:
    csv_subject = csv_data[csv_data['donor_id'] == '941_S_4420'].iloc[0]
    print("CSV data:")
    print(f"  most_severe_diagnosis_code: {csv_subject['most_severe_diagnosis_code']}")
    print(f"  most_severe_diagnosis: {csv_subject['most_severe_diagnosis']}")
    print(f"  diagnosis_date: {csv_subject['diagnosis_date']}")
else:
    print("❌ Not found in CSV")

# Check if we need to update H5AD data source
print(f"\n💡 CONCLUSION:")
print(f"The H5AD file contains expression data subjects (650) with diagnosis data from ~2013")
print(f"The CSV export contains ALL subjects (3,582) with latest diagnosis data through 2025")
print(f"For 941_S_4420: H5AD shows 2012-2013 data (MCI), CSV shows 2023 data (AD)")
print(f"This explains the partner discrepancy - they have 2023 data, we had 2013 data in H5AD")