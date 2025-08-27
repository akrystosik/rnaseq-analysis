import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

# Define paths
data_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data"
datasets = {
    "adni": f"{data_dir}/adni_standardized_v2.h5ad",
    "encode": f"{data_dir}/encode_standardized_v2.h5ad",
    "entex": f"{data_dir}/entex_standardized_v2.h5ad",
    "gtex": f"{data_dir}/gtex_standardized_v2.h5ad",
    "mage": f"{data_dir}/mage_standardized_v2.h5ad",
    "combined": f"{data_dir}/combined_all_genes_standardized.h5ad"
}

# Function to analyze gene IDs
def analyze_gene_ids(dataset_name, file_path):
    print(f"\n=== Analyzing {dataset_name} dataset ===")
    
    try:
        # Load the dataset
        adata = sc.read_h5ad(file_path)
        print(f"Dataset shape: {adata.shape}")
        
        # Get gene IDs
        gene_ids = adata.var_names.tolist()
        
        # Check if gene IDs are in var DataFrame
        if 'gene_id' in adata.var.columns:
            gene_ids_in_var = adata.var['gene_id'].tolist()
            print(f"Gene IDs in var['gene_id']: {len(gene_ids_in_var)}")
            print(f"First 5 gene IDs in var['gene_id']: {gene_ids_in_var[:5]}")
            
            # Check for mismatches between var_names and var['gene_id']
            mismatches = sum(1 for i, gene in enumerate(gene_ids) if str(gene) != str(gene_ids_in_var[i]))
            print(f"Mismatches between var_names and var['gene_id']: {mismatches}")
        
        # Analyze gene ID patterns
        print(f"First 5 gene IDs: {gene_ids[:5]}")
        
        # Check for ENSG pattern (Ensembl)
        ensembl_count = sum(1 for gene in gene_ids if str(gene).startswith('ENSG'))
        print(f"Ensembl IDs (ENSG): {ensembl_count} ({ensembl_count/len(gene_ids)*100:.2f}%)")
        
        # Check for numeric IDs
        numeric_count = sum(1 for gene in gene_ids if str(gene).isdigit())
        print(f"Numeric IDs: {numeric_count} ({numeric_count/len(gene_ids)*100:.2f}%)")
        
        # Check for other patterns
        other_count = len(gene_ids) - ensembl_count - numeric_count
        print(f"Other ID formats: {other_count} ({other_count/len(gene_ids)*100:.2f}%)")
        
        if other_count > 0:
            other_examples = [gene for gene in gene_ids if not str(gene).startswith('ENSG') and not str(gene).isdigit()][:5]
            print(f"Examples of other formats: {other_examples}")
        
        # Check var columns
        print(f"Var columns: {list(adata.var.columns)}")
        
        # Check for gene metadata
        has_gene_name = 'gene_name' in adata.var.columns
        has_chromosome = 'chromosome' in adata.var.columns
        has_gene_type = 'gene_type' in adata.var.columns
        
        print(f"Has gene_name: {has_gene_name}")
        print(f"Has chromosome: {has_chromosome}")
        print(f"Has gene_type: {has_gene_type}")
        
        if has_gene_name:
            empty_names = adata.var['gene_name'].isna().sum() + sum(1 for name in adata.var['gene_name'] if name == '')
            print(f"Empty gene_name: {empty_names} ({empty_names/len(adata.var)*100:.2f}%)")
        
        return adata
    
    except Exception as e:
        print(f"Error analyzing {dataset_name}: {e}")
        return None

# Analyze each dataset
dataset_adatas = {}
for name, path in datasets.items():
    if os.path.exists(path):
        adata = analyze_gene_ids(name, path)
        if adata is not None:
            dataset_adatas[name] = adata
    else:
        print(f"File not found: {path}")

# Analyze gene overlap between datasets
if len(dataset_adatas) >= 2 and 'combined' in dataset_adatas:
    print("\n=== Analyzing Gene Overlap ===")
    combined = dataset_adatas['combined']
    
    # Get gene IDs from each dataset
    gene_sets = {}
    for name, adata in dataset_adatas.items():
        if name != 'combined':
            gene_sets[name] = set(adata.var_names.astype(str))
    
    # Calculate overlap with combined dataset
    for name, genes in gene_sets.items():
        combined_genes = set(combined.var_names.astype(str))
        overlap = genes.intersection(combined_genes)
        print(f"Overlap {name} with combined: {len(overlap)}/{len(genes)} ({len(overlap)/len(genes)*100:.2f}%)")
        
        # Check genes not in combined dataset
        missing = genes - combined_genes
        if missing:
            print(f"Sample of genes in {name} but not in combined: {list(missing)[:5]}")