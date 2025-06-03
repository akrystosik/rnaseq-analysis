import scanpy as sc
import pandas as pd

# Load the combined dataset
adata = sc.read_h5ad('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_all_genes_standardized.h5ad')

# Examine the var dataframe (gene annotations)
print(f"Shape of adata.var: {adata.var.shape}")
print(f"Columns in adata.var: {list(adata.var.columns)}")

# Look at the first few rows
print("\nFirst 10 rows of adata.var:")
print(adata.var.head(10))

# Check for potential issues
print("\nChecking for potential issues:")
print(f"Missing values in gene_id: {adata.var['gene_id'].isna().sum() if 'gene_id' in adata.var.columns else 'Column not found'}")

# Check how many genes come from each dataset
if 'present_in_datasets' in adata.var.columns:
    print("\nDataset origin counts:")
    for dataset in ['encode', 'entex', 'gtex', 'adni', 'mage']:
        count = adata.var['present_in_datasets'].str.contains(dataset).sum()
        print(f"  - {dataset}: {count} genes")

# Check for unmapped genes
if 'mapping_source' in adata.var.columns:
    unmapped = adata.var[adata.var['mapping_source'] == 'unmapped']
    print(f"\nUnmapped genes: {len(unmapped)} ({len(unmapped)/len(adata.var)*100:.2f}%)")
    if len(unmapped) > 0:
        print("Sample of unmapped genes:")
        print(unmapped.head())

# Check if all genes have a gene_type
if 'gene_type' in adata.var.columns:
    missing_type = adata.var[adata.var['gene_type'].isna()]
    print(f"\nGenes missing gene_type: {len(missing_type)} ({len(missing_type)/len(adata.var)*100:.2f}%)")