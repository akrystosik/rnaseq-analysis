#!/usr/bin/env python3
"""
Create minimal test datasets for rapid testing
"""
import os
import scanpy as sc
import pandas as pd

# Base paths
BASE_DIR = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
OUTPUT_DIR = f"{BASE_DIR}/standardized_data/minimal_test"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Creating minimal test datasets...")

# Create minimal ENCODE dataset
try:
    # Load original ENCODE dataset
    encode_file = f"{BASE_DIR}/standardized_data/latest/encode_standardized_v1.h5ad"
    if os.path.exists(encode_file):
        adata = sc.read_h5ad(encode_file)
        
        # Select just 10-20 samples and 1000 genes
        mini_encode = adata[:20, :1000].copy()
        
        # Save minimal dataset
        mini_encode.write_h5ad(f"{OUTPUT_DIR}/encode_standardized_v1.h5ad")
        print(f"Created minimal ENCODE dataset: {mini_encode.shape}")
except Exception as e:
    print(f"Error creating minimal ENCODE dataset: {e}")

# Create minimal MAGE dataset
try:
    # Load original MAGE dataset
    mage_file = f"{BASE_DIR}/standardized_data/latest/mage_standardized_v1.h5ad"
    if os.path.exists(mage_file):
        adata = sc.read_h5ad(mage_file)
        
        # Select just 10-20 samples and 1000 genes
        mini_mage = adata[:20, :1000].copy()
        
        # Save minimal dataset
        mini_mage.write_h5ad(f"{OUTPUT_DIR}/mage_standardized_v1.h5ad")
        print(f"Created minimal MAGE dataset: {mini_mage.shape}")
except Exception as e:
    print(f"Error creating minimal MAGE dataset: {e}")

# Create minimal ADNI dataset
try:
    # Load original ADNI dataset
    adni_file = f"{BASE_DIR}/standardized_data/latest/adni_standardized_v1.h5ad"
    if os.path.exists(adni_file):
        adata = sc.read_h5ad(adni_file)
        
        # Select just 10-20 samples and 1000 genes
        mini_adni = adata[:20, :1000].copy()
        
        # Save minimal dataset
        mini_adni.write_h5ad(f"{OUTPUT_DIR}/adni_standardized_v1.h5ad")
        print(f"Created minimal ADNI dataset: {mini_adni.shape}")
except Exception as e:
    print(f"Error creating minimal ADNI dataset: {e}")

print("Minimal test datasets created in", OUTPUT_DIR)
