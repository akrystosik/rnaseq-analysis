import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the combined dataset
adata = sc.read_h5ad('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_all_genes_standardized.h5ad')

print(f"Dataset shape: {adata.shape}")

# Check gene metadata completeness
print("\nGene metadata completeness:")
for col in ['gene_name', 'gene_type', 'chromosome']:
    non_empty = sum(1 for x in adata.var[col] if x and str(x).strip() != '')
    percentage = non_empty / len(adata.var) * 100
    print(f"  - {col}: {non_empty}/{len(adata.var)} ({percentage:.2f}%)")

# Check mapping source distribution
print("\nMapping source distribution:")
mapping_counts = adata.var['mapping_source'].value_counts()
for source, count in mapping_counts.items():
    percentage = count / len(adata.var) * 100
    print(f"  - {source}: {count} ({percentage:.2f}%)")

# Check dataset sample distribution
print("\nSample distribution by dataset:")
dataset_counts = adata.obs['dataset'].value_counts()
for dataset, count in dataset_counts.items():
    percentage = count / len(adata.obs) * 100
    print(f"  - {dataset}: {count} ({percentage:.2f}%)")

# Check gene presence by dataset
print("\nGene presence by dataset:")
for dataset in ['adni', 'encode', 'entex', 'gtex', 'mage']:
    count = sum(1 for x in adata.var['present_in_datasets'] if dataset in str(x))
    percentage = count / len(adata.var) * 100
    print(f"  - {dataset}: {count}/{len(adata.var)} ({percentage:.2f}%)")

# Check ensembl ID vs. other ID format distribution
ensembl_count = sum(1 for x in adata.var_names if str(x).startswith('ENSG'))
spike_in_count = sum(1 for x in adata.var_names if str(x).startswith('gSpikein'))
other_count = len(adata.var) - ensembl_count - spike_in_count

print("\nGene ID format distribution:")
print(f"  - Ensembl IDs: {ensembl_count}/{len(adata.var)} ({ensembl_count/len(adata.var)*100:.2f}%)")
print(f"  - Spike-in controls: {spike_in_count}/{len(adata.var)} ({spike_in_count/len(adata.var)*100:.2f}%)")
print(f"  - Other IDs: {other_count}/{len(adata.var)} ({other_count/len(adata.var)*100:.2f}%)")

# Check dataset-specific metrics
print("\nData type distribution:")
data_type_counts = adata.obs['data_type'].value_counts()
for data_type, count in data_type_counts.items():
    percentage = count / len(adata.obs) * 100
    print(f"  - {data_type}: {count} ({percentage:.2f}%)")

# Check for missing values in key metadata
print("\nMissing values in key metadata:")
for col in ['subject_id', 'tissue', 'sex', 'age']:
    missing = adata.obs[col].isna().sum() + sum(1 for x in adata.obs[col] if str(x).strip() == '')
    percentage = missing / len(adata.obs) * 100
    print(f"  - {col}: {missing}/{len(adata.obs)} ({percentage:.2f}%)")