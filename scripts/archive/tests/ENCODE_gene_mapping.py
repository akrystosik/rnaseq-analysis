import scanpy as sc
import pandas as pd
import numpy as np

# Load the preprocessed ENCODE dataset
encode_data = sc.read_h5ad('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/encode_standardized_preprocessed.h5ad')

# Load the combined dataset
combined_data = sc.read_h5ad('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_all_genes_standardized.h5ad')

# Load the entrez to ensembl mapping
entrez_mapping = pd.read_csv('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/entrez_to_ensembl_mapping.csv')

print(f"ENCODE dataset shape: {encode_data.shape}")
print(f"Combined dataset shape: {combined_data.shape}")
print(f"Entrez to Ensembl mapping entries: {len(entrez_mapping)}")

# Check ENCODE var columns
print("\nENCODE var columns:")
print(encode_data.var.columns)

# Check ENCODE var_names format
encode_var_names = list(encode_data.var_names[:5])
print(f"\nENCODE var_names sample: {encode_var_names}")

# Check ensembl_id column in ENCODE var
if 'ensembl_id' in encode_data.var.columns:
    # Count non-empty ensembl_ids
    non_empty = sum(1 for x in encode_data.var['ensembl_id'] if x and str(x).strip() != '')
    percentage = non_empty / len(encode_data.var) * 100
    print(f"\nENCODE genes with ensembl_id: {non_empty}/{len(encode_data.var)} ({percentage:.2f}%)")
    
    # Sample of ensembl_ids
    sample_ensembl_ids = encode_data.var['ensembl_id'].dropna().iloc[:5].tolist()
    print(f"Sample ENCODE ensembl_ids: {sample_ensembl_ids}")

# Check how many ENCODE genes are in the combined dataset
encode_genes_in_combined = sum(1 for x in combined_data.var['present_in_datasets'] if 'encode' in str(x))
print(f"\nENCODE genes in combined dataset: {encode_genes_in_combined}/{len(combined_data.var)} ({encode_genes_in_combined/len(combined_data.var)*100:.2f}%)")

# Check numeric IDs in entrez mapping
if 'entrez_id' in entrez_mapping.columns:
    # Count how many ENCODE numeric IDs are in the entrez mapping
    encode_numeric_ids = set(encode_data.var_names.astype(str))
    entrez_ids = set(entrez_mapping['entrez_id'].astype(str))
    
    overlap = encode_numeric_ids.intersection(entrez_ids)
    print(f"\nENCODE numeric IDs found in entrez mapping: {len(overlap)}/{len(encode_numeric_ids)} ({len(overlap)/len(encode_numeric_ids)*100:.2f}%)")