#!/usr/bin/env python
# Script to investigate missing datasets and placeholder gene IDs

import os
import scanpy as sc
import pandas as pd
import numpy as np

# Base directory
BASE_DIR = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
LATEST_RUN = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250428_234225"

# Check for individual dataset files
print("=== CHECKING INDIVIDUAL DATASET FILES ===")
for dataset in ['adni', 'encode', 'entex', 'gtex', 'mage']:
    # Check both standardized and preprocessed versions
    std_path = os.path.join(LATEST_RUN, f"{dataset}_standardized.h5ad")
    prep_path = os.path.join(BASE_DIR, f"preprocessed_data/run_20250428_234225/{dataset}_standardized_preprocessed.h5ad")
    
    # Check standardized file
    if os.path.exists(std_path):
        try:
            std_data = sc.read_h5ad(std_path)
            print(f"{dataset} standardized file exists: {std_path}")
            print(f"  - Shape: {std_data.shape[0]} samples × {std_data.shape[1]} genes")
        except Exception as e:
            print(f"Error reading {dataset} standardized file: {str(e)}")
    else:
        print(f"{dataset} standardized file NOT FOUND: {std_path}")
    
    # Check preprocessed file
    if os.path.exists(prep_path):
        try:
            prep_data = sc.read_h5ad(prep_path)
            print(f"{dataset} preprocessed file exists: {prep_path}")
            print(f"  - Shape: {prep_data.shape[0]} samples × {prep_data.shape[1]} genes")
        except Exception as e:
            print(f"Error reading {dataset} preprocessed file: {str(e)}")
    else:
        print(f"{dataset} preprocessed file NOT FOUND: {prep_path}")

# Load the combined dataset to analyze placeholder genes
print("\n=== ANALYZING PLACEHOLDER GENE IDs ===")
combined_path = os.path.join(LATEST_RUN, "combined_all_genes_sparse_standardized.h5ad")
adata = sc.read_h5ad(combined_path)

# Find placeholder genes and analyze their metadata
placeholder_mask = [True if 'placeholder' in str(x).lower() else False for x in adata.var_names]
placeholder_count = sum(placeholder_mask)
print(f"Total placeholder genes: {placeholder_count}")

if placeholder_count > 0:
    placeholder_data = adata.var[placeholder_mask]
    
    # 1. Check which datasets contain these placeholder genes
    if 'present_in_datasets' in placeholder_data.columns:
        dataset_counts = {}
        for dataset in ['adni', 'encode', 'entex', 'gtex', 'mage']:
            count = sum(1 for x in placeholder_data['present_in_datasets'] if dataset in str(x))
            percentage = count / len(placeholder_data) * 100
            dataset_counts[dataset] = count
            print(f"Placeholder genes in {dataset}: {count}/{placeholder_count} ({percentage:.2f}%)")
    
    # 2. Check if there's a mapping source pattern
    if 'mapping_source' in placeholder_data.columns:
        source_counts = placeholder_data['mapping_source'].value_counts()
        print("\nPlaceholder gene mapping sources:")
        for source, count in source_counts.items():
            percentage = count / len(placeholder_data) * 100
            print(f"  - {source}: {count} ({percentage:.2f}%)")
    
    # 3. Check for any original IDs that might give clues
    if 'original_ids' in placeholder_data.columns:
        print("\nSample of original IDs for placeholder genes:")
        sample_size = min(10, len(placeholder_data))
        for i, (idx, row) in enumerate(placeholder_data.head(sample_size).iterrows()):
            print(f"  {i+1}. {idx}: {row.get('original_ids', 'No original ID')}")

# Check the gene ID mapping reference file
print("\n=== CHECKING GENE ID REFERENCE MAPPING ===")
mapping_file = f"{BASE_DIR}/metadata/json/gene_id_reference_mapping.csv"
if os.path.exists(mapping_file):
    try:
        mapping_df = pd.read_csv(mapping_file)
        print(f"Gene ID reference mapping file exists: {mapping_file}")
        print(f"  - Total mappings: {len(mapping_df)}")
        
        # Check for columns
        print(f"  - Columns: {list(mapping_df.columns)}")
        
        # Look for problematic mappings (missing Ensembl IDs)
        if 'ensembl_id' in mapping_df.columns:
            missing_ensembl = mapping_df['ensembl_id'].isna().sum()
            missing_pct = missing_ensembl / len(mapping_df) * 100
            print(f"  - Mappings missing Ensembl ID: {missing_ensembl}/{len(mapping_df)} ({missing_pct:.2f}%)")
        
            # Sample problematic mappings
            if missing_ensembl > 0:
                print("\nSample of entries missing Ensembl IDs:")
                problematic = mapping_df[mapping_df['ensembl_id'].isna()]
                for i, (_, row) in enumerate(problematic.head(5).iterrows()):
                    print(f"  {i+1}. {row.to_dict()}")
    except Exception as e:
        print(f"Error reading gene ID reference mapping: {str(e)}")
else:
    print(f"Gene ID reference mapping file NOT FOUND: {mapping_file}")

# Check the ADNI and ENTEX raw data directories to confirm they exist
print("\n=== CHECKING RAW DATA DIRECTORIES ===")
adni_dir = f"{BASE_DIR}/adni_microarray"
entex_dir = f"{BASE_DIR}/encode/entex"

if os.path.exists(adni_dir):
    print(f"ADNI directory exists: {adni_dir}")
    print(f"  - Contents: {os.listdir(adni_dir)[:5]} ...")
else:
    print(f"ADNI directory NOT FOUND: {adni_dir}")

if os.path.exists(entex_dir):
    print(f"ENTEX directory exists: {entex_dir}")
    print(f"  - Contents: {os.listdir(entex_dir)[:5]} ...")
else:
    print(f"ENTEX directory NOT FOUND: {entex_dir}")

# Examine log file for errors
print("\n=== EXAMINING PIPELINE LOG ===")
log_file = f"{BASE_DIR}/logs/pipeline_20250428_234225.log"
if os.path.exists(log_file):
    print(f"Log file exists: {log_file}")
    
    # Look for error messages
    with open(log_file, 'r') as f:
        log_content = f.read()
    
    error_lines = [line for line in log_content.split('\n') if 'error' in line.lower() or 'failed' in line.lower() or 'warning' in line.lower()]
    if error_lines:
        print("Found potential error messages in log:")
        for i, line in enumerate(error_lines[:10]):  # Show first 10 errors
            print(f"  {i+1}. {line}")
        if len(error_lines) > 10:
            print(f"  ... and {len(error_lines) - 10} more error/warning messages")
    else:
        print("No obvious error messages found in log file")
else:
    print(f"Log file NOT FOUND: {log_file}")