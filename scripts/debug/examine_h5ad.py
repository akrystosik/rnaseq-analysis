#!/usr/bin/env python
# Diagnostic script to examine the combined h5ad dataset

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the combined dataset
file_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250428_234225/combined_all_genes_sparse_standardized.h5ad"
print(f"Loading dataset from: {file_path}")
adata = sc.read_h5ad(file_path)

# Basic dataset stats
print("\n=== BASIC DATASET STATS ===")
print(f"Dataset shape: {adata.shape[0]} samples × {adata.shape[1]} genes")

# Check observation metadata (samples)
print("\n=== SAMPLE METADATA OVERVIEW ===")
print("Observation (sample) keys:", list(adata.obs.columns))
print("\nDataset distribution:")
dataset_counts = adata.obs['dataset'].value_counts()
for dataset, count in dataset_counts.items():
    percentage = count / len(adata.obs) * 100
    print(f"  - {dataset}: {count} ({percentage:.2f}%)")

# Check tissue annotation status
print("\n=== TISSUE ANNOTATION STATUS ===")
if 'tissue' in adata.obs.columns:
    missing_tissue = adata.obs['tissue'].isna().sum()
    missing_pct = (missing_tissue / len(adata.obs)) * 100
    print(f"Samples missing tissue annotation: {missing_tissue}/{len(adata.obs)} ({missing_pct:.2f}%)")
    
    print("\nTop 10 tissues:")
    top_tissues = adata.obs['tissue'].value_counts().head(10)
    for tissue, count in top_tissues.items():
        print(f"  - {tissue}: {count}")
else:
    print("ERROR: 'tissue' column not found in sample metadata!")

# Check cell line annotation status
print("\n=== CELL LINE ANNOTATION STATUS ===")
if 'cell_line' in adata.obs.columns:
    has_cell_line = (~adata.obs['cell_line'].isna()).sum()
    has_cell_line_pct = (has_cell_line / len(adata.obs)) * 100
    print(f"Samples with cell line annotation: {has_cell_line}/{len(adata.obs)} ({has_cell_line_pct:.2f}%)")
    
    if has_cell_line > 0:
        print("\nTop 10 cell lines:")
        top_cell_lines = adata.obs['cell_line'].value_counts().head(10)
        for cell_line, count in top_cell_lines.items():
            if pd.notna(cell_line):  # Filter out NaN values
                print(f"  - {cell_line}: {count}")
else:
    print("INFO: 'cell_line' column not found in sample metadata")

# Check gene metadata
print("\n=== GENE METADATA OVERVIEW ===")
print("Variable (gene) keys:", list(adata.var.columns))

# Check gene ID distribution
print("\n=== GENE ID FORMAT DISTRIBUTION ===")
# Using var_names (index of var dataframe)
ensembl_count = sum(1 for x in adata.var_names if str(x).startswith('ENSG'))
spike_in_count = sum(1 for x in adata.var_names if str(x).startswith('gSpikein'))
placeholder_count = sum(1 for x in adata.var_names if 'placeholder' in str(x).lower())
other_count = len(adata.var) - ensembl_count - spike_in_count - placeholder_count

print(f"Total genes: {len(adata.var)}")
print(f"  - Ensembl IDs: {ensembl_count} ({ensembl_count/len(adata.var)*100:.2f}%)")
print(f"  - Spike-in controls: {spike_in_count} ({spike_in_count/len(adata.var)*100:.2f}%)")
print(f"  - Placeholder IDs: {placeholder_count} ({placeholder_count/len(adata.var)*100:.2f}%)")
print(f"  - Other IDs: {other_count} ({other_count/len(adata.var)*100:.2f}%)")

# Sample some problematic gene IDs
if placeholder_count > 0:
    print("\nSample of placeholder gene IDs:")
    placeholder_ids = [x for x in adata.var_names if 'placeholder' in str(x).lower()]
    for i, pid in enumerate(placeholder_ids[:10]):
        print(f"  {i+1}. {pid}")

# Check gene presence by dataset
print("\n=== GENE PRESENCE BY DATASET ===")
if 'present_in_datasets' in adata.var.columns:
    dataset_gene_counts = {}
    for dataset in ['adni', 'encode', 'entex', 'gtex', 'mage']:
        count = sum(1 for x in adata.var['present_in_datasets'] if dataset in str(x))
        percentage = count / len(adata.var) * 100
        print(f"  - {dataset}: {count}/{len(adata.var)} ({percentage:.2f}%)")
        dataset_gene_counts[dataset] = count
        
    # Check dataset-unique genes
    print("\nDataset-unique genes:")
    for dataset in ['adni', 'encode', 'entex', 'gtex', 'mage']:
        unique_count = sum(1 for x in adata.var['present_in_datasets'] 
                          if str(x) == f"['{dataset}']" or str(x) == f"[\"{dataset}\"]")
        unique_pct = unique_count / dataset_gene_counts[dataset] * 100 if dataset_gene_counts[dataset] > 0 else 0
        print(f"  - {dataset}: {unique_count} ({unique_pct:.2f}% of dataset genes)")
else:
    print("ERROR: 'present_in_datasets' column not found in gene metadata!")

# Check if gene_id and ensembl_id columns exist
print("\n=== GENE ID MAPPING STATUS ===")
if 'gene_id' in adata.var.columns:
    missing_gene_id = adata.var['gene_id'].isna().sum()
    missing_pct = (missing_gene_id / len(adata.var)) * 100
    print(f"Genes missing gene_id: {missing_gene_id}/{len(adata.var)} ({missing_pct:.2f}%)")
else:
    print("WARNING: 'gene_id' column not found in gene metadata!")

if 'ensembl_id' in adata.var.columns:
    missing_ensembl = adata.var['ensembl_id'].isna().sum()
    missing_pct = (missing_ensembl / len(adata.var)) * 100
    print(f"Genes missing ensembl_id: {missing_ensembl}/{len(adata.var)} ({missing_pct:.2f}%)")
else:
    print("INFO: 'ensembl_id' column not found in gene metadata")

# Check gene name information
if 'gene_name' in adata.var.columns:
    missing_gene_name = adata.var['gene_name'].isna().sum()
    missing_pct = (missing_gene_name / len(adata.var)) * 100
    print(f"Genes missing gene_name: {missing_gene_name}/{len(adata.var)} ({missing_pct:.2f}%)")
else:
    print("INFO: 'gene_name' column not found in gene metadata")

# Check data sparsity
print("\n=== DATA SPARSITY ===")
sparsity = 1 - (adata.X.nnz / (adata.shape[0] * adata.shape[1]))
print(f"Matrix sparsity: {sparsity:.4f} ({sparsity*100:.2f}%)")

# Check for zeros
zero_count = (adata.X == 0).sum()
if hasattr(zero_count, 'compute'):  # If dask array
    zero_count = zero_count.compute()
zero_percent = zero_count / (adata.shape[0] * adata.shape[1]) * 100
print(f"Zero values: {zero_count}/{adata.shape[0] * adata.shape[1]} ({zero_percent:.2f}%)")

# Check individual dataset files
print("\n=== INDIVIDUAL DATASET STATUS ===")
import os
dataset_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250428_234225"
for dataset_name in ['encode', 'gtex', 'mage', 'adni', 'entex']:
    filename = f"{dataset_name}_standardized.h5ad"
    fullpath = os.path.join(dataset_path, filename)
    if os.path.exists(fullpath):
        try:
            ds_adata = sc.read_h5ad(fullpath)
            print(f"{dataset_name}: {ds_adata.shape[0]} samples × {ds_adata.shape[1]} genes")
            if 'tissue' in ds_adata.obs.columns:
                missing = ds_adata.obs['tissue'].isna().sum()
                print(f"  - Missing tissue annotations: {missing}/{ds_adata.shape[0]} ({missing/ds_adata.shape[0]*100:.2f}%)")
            else:
                print(f"  - No 'tissue' column found!")
                
            # Check gene ID format
            ensembl_count = sum(1 for x in ds_adata.var_names if str(x).startswith('ENSG'))
            print(f"  - Genes with Ensembl IDs: {ensembl_count}/{ds_adata.shape[1]} ({ensembl_count/ds_adata.shape[1]*100:.2f}%)")
        except Exception as e:
            print(f"Error loading {dataset_name}: {str(e)}")
    else:
        print(f"{dataset_name}: File not found at {fullpath}")