#!/usr/bin/env python3
"""
Debug ADNI sex data issue
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def debug_adni_sex_data():
    """Debug what's happening to ADNI sex data during loading"""
    
    print("="*60)
    print("DEBUGGING ADNI SEX DATA ISSUE")
    print("="*60)
    
    # Try different ADNI dataset files
    adni_paths = [
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/adni_standardized_preprocessed.h5ad",
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/adni_standardized_v2.h5ad",
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/latest_v2.2/adni_standardized_v2.h5ad",
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/adni_standardized.h5ad"
    ]
    
    for adni_path in adni_paths:
        if Path(adni_path).exists():
            print(f"\n{'='*50}")
            print(f"Checking: {adni_path}")
            print(f"{'='*50}")
            
            try:
                adata = sc.read_h5ad(adni_path)
                print(f"Dataset shape: {adata.shape}")
                
                if 'sex' in adata.obs.columns:
                    print("Sex column found!")
                    print("Sex value counts:")
                    print(adata.obs['sex'].value_counts())
                    print("Sex unique values:")
                    print(f"  {list(adata.obs['sex'].unique())}")
                    print("Sex data type:")
                    print(f"  {adata.obs['sex'].dtype}")
                    
                    # Check for null values
                    null_count = adata.obs['sex'].isna().sum()
                    print(f"Null values in sex: {null_count}")
                    
                    # Sample some values
                    print("Sample sex values:")
                    print(adata.obs['sex'].head(10).tolist())
                    
                else:
                    print("❌ No 'sex' column found!")
                    print("Available columns:")
                    for col in adata.obs.columns:
                        if 'sex' in col.lower():
                            print(f"  - {col}: {adata.obs[col].dtype}")
                
            except Exception as e:
                print(f"❌ Error loading {adni_path}: {str(e)}")
                continue
    
    print("\n" + "="*60)
    print("CHECKING COMBINED DATASET")
    print("="*60)
    
    # Also check the combined dataset to see if ADNI samples have sex data there
    combined_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/latest_v2.2/combined_dataset_all_genes_sparse.h5ad"
    
    if Path(combined_path).exists():
        print(f"Loading combined dataset: {combined_path}")
        adata_combined = sc.read_h5ad(combined_path)
        
        # Filter to ADNI samples
        adni_mask = adata_combined.obs['dataset'] == 'adni'
        adni_samples = adata_combined[adni_mask]
        
        print(f"ADNI samples in combined dataset: {adni_samples.n_obs}")
        
        if 'sex' in adni_samples.obs.columns:
            print("Sex data in ADNI samples from combined dataset:")
            print(adni_samples.obs['sex'].value_counts())
            print("Unique values:")
            print(f"  {list(adni_samples.obs['sex'].unique())}")
        else:
            print("❌ No sex column in combined dataset")

if __name__ == "__main__":
    debug_adni_sex_data()