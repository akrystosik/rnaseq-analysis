#!/usr/bin/env python
# Script to diagnose and fix critical issues in the RNA-seq pipeline

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import shutil
import json

# Base paths
BASE_DIR = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
LATEST_RUN = f"{BASE_DIR}/standardized_data/run_20250428_234225"
SCRIPTS_DIR = f"{BASE_DIR}/scripts"
METADATA_DIR = f"{BASE_DIR}/metadata/json"

print("=== RNA-seq Pipeline Diagnostic and Fix Script ===")

# PART 1: Examine the ENCODE dataset to understand placeholder genes
print("\n=== EXAMINING ENCODE DATASET ===")
encode_v2_path = f"{LATEST_RUN}/encode_standardized_v2.h5ad"

if os.path.exists(encode_v2_path):
    print(f"Loading ENCODE dataset: {encode_v2_path}")
    encode_data = sc.read_h5ad(encode_v2_path)
    print(f"ENCODE dataset shape: {encode_data.shape}")
    
    # Check gene IDs
    print("\nGene ID format distribution:")
    ensembl_count = sum(1 for x in encode_data.var_names if str(x).startswith('ENSG'))
    placeholder_count = sum(1 for x in encode_data.var_names if 'placeholder' in str(x).lower())
    other_count = len(encode_data.var) - ensembl_count - placeholder_count
    
    print(f"  - Ensembl IDs: {ensembl_count}/{len(encode_data.var)} ({ensembl_count/len(encode_data.var)*100:.2f}%)")
    print(f"  - Placeholder IDs: {placeholder_count}/{len(encode_data.var)} ({placeholder_count/len(encode_data.var)*100:.2f}%)")
    print(f"  - Other format IDs: {other_count}/{len(encode_data.var)} ({other_count/len(encode_data.var)*100:.2f}%)")
    
    # Sample placeholder gene info
    if placeholder_count > 0:
        placeholder_mask = [True if 'placeholder' in str(x).lower() else False for x in encode_data.var_names]
        placeholder_genes = encode_data.var[placeholder_mask]
        
        print("\nSample placeholder genes from ENCODE:")
        for i, (idx, row) in enumerate(placeholder_genes.head(5).iterrows()):
            gene_info = {col: row[col] for col in encode_data.var.columns if col in row}
            print(f"  {i+1}. {idx}: {gene_info}")
        
        # Check if there's any additional metadata for these genes
        if 'original_ids' in encode_data.var.columns:
            orig_ids = encode_data.var.loc[placeholder_mask, 'original_ids']
            has_real_ids = sum(1 for x in orig_ids if x != x.index)
            print(f"\nPlaceholder genes with real original IDs: {has_real_ids}/{placeholder_count}")
            
            # Sample ones with different original IDs
            different_ids = orig_ids[orig_ids != orig_ids.index]
            if len(different_ids) > 0:
                print("Sample of placeholder genes with real original IDs:")
                for i, (idx, orig_id) in enumerate(different_ids.head(5).items()):
                    print(f"  {i+1}. Placeholder: {idx}, Original ID: {orig_id}")
else:
    print(f"ENCODE dataset file not found: {encode_v2_path}")

# PART 2: Examine the gene ID mapping reference
print("\n=== EXAMINING GENE ID MAPPING REFERENCE ===")
mapping_file = f"{METADATA_DIR}/gene_id_reference_mapping.csv"

if os.path.exists(mapping_file):
    print(f"Loading gene ID mapping reference: {mapping_file}")
    mapping_df = pd.read_csv(mapping_file, low_memory=False)
    print(f"Mapping file shape: {mapping_df.shape}")
    
    # Check for placeholder entries in the mapping
    placeholder_in_mapping = sum(1 for x in mapping_df['gene_id'] if isinstance(x, str) and 'placeholder' in x.lower())
    print(f"Placeholder IDs in mapping file: {placeholder_in_mapping}")
    
    # Check mapping confidence distribution
    if 'mapping_confidence' in mapping_df.columns:
        confidence_counts = mapping_df['mapping_confidence'].value_counts()
        print("\nMapping confidence distribution:")
        for conf, count in confidence_counts.items():
            pct = count / len(mapping_df) * 100
            print(f"  - {conf}: {count} ({pct:.2f}%)")
    
    # Examine potentially problematic mappings
    low_conf_mask = mapping_df['mapping_confidence'] == 'low' if 'mapping_confidence' in mapping_df.columns else None
    if low_conf_mask is not None and low_conf_mask.sum() > 0:
        low_conf = mapping_df[low_conf_mask]
        print(f"\nLow confidence mappings: {len(low_conf)}")
        print("Sample of low confidence mappings:")
        for i, (_, row) in enumerate(low_conf.head(5).iterrows()):
            print(f"  {i+1}. {row['gene_id']}: {row.to_dict()}")
else:
    print(f"Gene ID mapping reference file not found: {mapping_file}")

# PART 3: Examine the ENCODE mapping script
print("\n=== EXAMINING ENCODE MAPPING PROCESS ===")
mapping_script = f"{SCRIPTS_DIR}/generate_encode_mapping.py"

if os.path.exists(mapping_script):
    print(f"ENCODE mapping script exists: {mapping_script}")
    
    # Check the script content
    with open(mapping_script, 'r') as f:
        script_content = f.read()
    
    # Look for placeholder generation logic
    placeholder_lines = [line for line in script_content.split('\n') if 'placeholder' in line.lower()]
    if placeholder_lines:
        print("\nLines mentioning placeholder in generate_encode_mapping.py:")
        for i, line in enumerate(placeholder_lines[:10]):
            print(f"  {i+1}. {line.strip()}")
else:
    print(f"ENCODE mapping script not found: {mapping_script}")

# PART 4: Examine the issue with missing files
print("\n=== EXAMINING FILE NAMING ISSUE ===")
dataset_files = os.listdir(LATEST_RUN)
print(f"Files in the latest run directory:")
for file in sorted(dataset_files):
    if file.endswith('.h5ad'):
        print(f"  - {file}")

# Check how the combined dataset script references the individual files
combined_script = f"{SCRIPTS_DIR}/create_combined_dataset_all_genes_sparse.py"
if os.path.exists(combined_script):
    print(f"\nChecking combined dataset script: {combined_script}")
    with open(combined_script, 'r') as f:
        combined_content = f.read()
    
    # Look for file path patterns
    file_pattern_lines = [line for line in combined_content.split('\n') 
                         if ('file_pattern' in line.lower() or 'h5ad' in line.lower()) 
                         and not line.strip().startswith('#')]
    if file_pattern_lines:
        print("Lines related to file patterns:")
        for i, line in enumerate(file_pattern_lines[:15]):
            print(f"  {i+1}. {line.strip()}")
else:
    print(f"Combined dataset script not found: {combined_script}")

# PART 5: Check the ADNI dataset issue
print("\n=== EXAMINING ADNI DATASET ISSUE ===")
adni_dir = f"{BASE_DIR}/adni_microarray"
if os.path.exists(adni_dir):
    print(f"ADNI directory exists: {adni_dir}")
    adni_files = os.listdir(adni_dir)
    print(f"Total files/dirs in ADNI directory: {len(adni_files)}")
    
    # Check if the standardize_datasets.py script handles ADNI
    standardize_script = f"{SCRIPTS_DIR}/standardize_datasets.py"
    if os.path.exists(standardize_script):
        with open(standardize_script, 'r') as f:
            std_content = f.read()
        
        adni_lines = [line for line in std_content.split('\n') 
                     if 'adni' in line.lower() and not line.strip().startswith('#')]
        if adni_lines:
            print("\nLines handling ADNI in standardize_datasets.py:")
            for i, line in enumerate(adni_lines[:15]):
                print(f"  {i+1}. {line.strip()}")
else:
    print(f"ADNI directory not found: {adni_dir}")

# PART 6: Suggested fixes
print("\n=== SUGGESTED FIXES ===")

# 1. Fix for file naming issues (wrong path references)
print("\n1. Create symbolic links with expected filenames:")
print("""
# Create symbolic links for the dataset files with the expected names
ln -s $LATEST_RUN/encode_standardized_v2.h5ad $LATEST_RUN/encode_standardized.h5ad
ln -s $LATEST_RUN/gtex_standardized_v2.h5ad $LATEST_RUN/gtex_standardized.h5ad
ln -s $LATEST_RUN/mage_standardized_v2.h5ad $LATEST_RUN/mage_standardized.h5ad
""")

# 2. Fix for ADNI processing
print("\n2. Check ADNI processing in standardize_datasets.py:")
print("""
# Add debug info to the ADNI processing section in standardize_datasets.py
# Look for lines like:
#   if args.adni_dir:
#       process_adni(args.adni_dir, args.output_dir)
# Add print statements to verify ADNI is being processed
""")

# 3. Fix for placeholder gene IDs in ENCODE
print("\n3. Improve gene ID mapping for ENCODE:")
print("""
# Modify the gene_id_mapping_reference.py script to better handle ENCODE gene IDs
# Look for placeholder generation logic and ensure it's using all available mapping sources
# Run the following command to regenerate the mapping with verbose output:
python $SCRIPTS_DIR/gene_id_mapping_reference.py \\
    --encode-dir "$BASE_DIR/encode/raw_data" \\
    --entex-dir "$BASE_DIR/encode/entex" \\
    --entrez-mapping "$METADATA_DIR/entrez_to_ensembl_mapping.csv" \\
    --output "$METADATA_DIR/gene_id_reference_mapping.csv" \\
    --verbose \\
    --force
""")

# 4. Cell line issues
print("\n4. Add cell line information to the metadata:")
print("""
# Modify the standardize_metadata.py script to ensure cell_line information is preserved
# Look for cell_line handling in ENCODE processing and ensure it's being included in the obs dataframe
""")

# 5. Fix the combined dataset generation
print("\n5. Fix combined dataset generation:")
print("""
# Modify create_combined_dataset_all_genes_sparse.py to properly handle file paths
# Ensure it's using the correct file patterns (should look for *_standardized_v2.h5ad or create symlinks)
# Add debug logging to show which files are being included
""")

# Provide a simplified script to create the symbolic links
print("\n=== FIX IMPLEMENTATION: CREATE SYMBOLIC LINKS ===")
fix_script = """#!/bin/bash
# Script to create symbolic links with correct filenames

# Set base directory
LATEST_RUN="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250428_234225"

# Create symbolic links
ln -sf "$LATEST_RUN/encode_standardized_v2.h5ad" "$LATEST_RUN/encode_standardized.h5ad"
ln -sf "$LATEST_RUN/gtex_standardized_v2.h5ad" "$LATEST_RUN/gtex_standardized.h5ad"
ln -sf "$LATEST_RUN/mage_standardized_v2.h5ad" "$LATEST_RUN/mage_standardized.h5ad"

echo "Created symbolic links with expected filenames"
echo "Next step: re-run the combined dataset generation script"
"""

# Write the fix script to a file
fix_script_path = f"{SCRIPTS_DIR}/debug/create_symlinks.sh"
os.makedirs(os.path.dirname(fix_script_path), exist_ok=True)
with open(fix_script_path, 'w') as f:
    f.write(fix_script)
os.chmod(fix_script_path, 0o755)  # Make executable

print(f"Created fix script: {fix_script_path}")
print("Run this script to create the symbolic links, then re-run the combined dataset step")