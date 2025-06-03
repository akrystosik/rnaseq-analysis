#!/usr/bin/env python3
"""Check what diagnosis data is in the current H5AD file for the problem subjects."""

import scanpy as sc
import pandas as pd

# Load current H5AD file  
h5ad_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
print(f"Loading H5AD file: {h5ad_file}")

try:
    adata = sc.read_h5ad(h5ad_file)
    print(f"✅ Loaded H5AD: {adata.n_obs} samples × {adata.n_vars} genes")
    
    problem_subjects = ['002_S_0413', '002_S_0729', '002_S_1155', '002_S_1261', '002_S_1268']
    
    print("\n📋 Checking diagnosis data for problem subjects:")
    diagnosis_cols = [col for col in adata.obs.columns if 'diagnosis' in col.lower()]
    print(f"Available diagnosis columns: {diagnosis_cols}")
    
    for subject in problem_subjects:
        # Find samples for this subject
        subject_samples = adata.obs[adata.obs['subject_id'] == subject]
        if len(subject_samples) > 0:
            print(f"\n🎯 Subject {subject} ({len(subject_samples)} samples):")
            sample = subject_samples.iloc[0]
            for col in diagnosis_cols:
                print(f"  {col}: {sample.get(col, 'MISSING')}")
        else:
            print(f"\n❌ Subject {subject}: Not found in H5AD file")
    
    print(f"\n📊 Overall diagnosis data summary:")
    if 'worst_diagnosis_label' in adata.obs.columns:
        diagnosis_counts = adata.obs['worst_diagnosis_label'].value_counts(dropna=False)
        print(diagnosis_counts)
    else:
        print("No 'worst_diagnosis_label' column found")
        
except Exception as e:
    print(f"❌ Error loading H5AD file: {e}")