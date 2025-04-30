# Multi-Dataset Gene Expression Analysis Pipeline

A comprehensive pipeline for integrating, standardizing, and querying gene expression data from multiple sources including ENCODE, GTEx, MAGE, ADNI, and ENTEx.

## Overview

This project provides tools for:

1. **Downloading and organizing** gene expression data from multiple sources
2. **Standardizing** data into a consistent format for cross-dataset analysis
3. **Querying and analyzing** gene expression across tissues, donors, and datasets
4. **Optimized data access** for model training and large-scale analysis

## Features

- Standardized data representation using AnnData format
- Cross-dataset gene expression querying
- High-performance data access optimized for repeated queries
- Tissue and donor-specific expression analysis
- Efficient caching for rapid access during model training
- Comprehensive metadata and annotation support
- **NEW: Combined dataset with ALL genes** using sparse matrix representation

## Directory Structure

```
rnaseq/
├── scripts/
│   ├── standardize_datasets.py                   # Data standardization script
│   ├── create_combined_dataset.py                # Create combined dataset (common genes)
│   ├── create_combined_dataset_all_genes_sparse.py # Create combined dataset (all genes, sparse)
│   ├── expression-query.py                       # Command-line query tool
│   ├── entex_pipeline.sh                         # ENTEx integration pipeline
│   ├── download_entex_files.py                   # ENTEx file downloader
│   ├── generate_entex_metadata.py                # ENTEx metadata generator
│   ├── optimized_gene_expression.py              # Optimized query functions
│   └── setup_workspace.sh                        # Environment setup script
├── encode/
│   ├── metadata/
│   │   ├── entex_metadata.json                   # ENTEx metadata
│   │   └── entex_manual_metadata.json            # Manual ENTEx annotations
│   └── raw_data/
│       ├── entex/                                # ENTEx raw data files
│       └── cell_lines/                           # ENCODE cell line data
├── gtex/
│   ├── metadata/                                 # GTEx metadata
│   └── raw_data/                                 # GTEx raw data
├── mage/
│   └── raw_data/                                 # MAGE raw data
├── adni/
│   └── raw_data/                                 # ADNI raw data
└── standardized_data/                            # Standardized AnnData files
    ├── encode_standardized.h5ad
    ├── entex_standardized.h5ad
    ├── gtex_standardized.h5ad
    ├── mage_standardized.h5ad
    ├── adni_standardized.h5ad
    ├── combined_standardized.h5ad                # Combined dataset (common genes only)
    └── combined_all_genes_standardized.h5ad      # Combined dataset (all genes, sparse matrix)
```

## Installation

### Requirements

- Python 3.6+
- anndata
- numpy
- pandas
- matplotlib
- seaborn
- requests
- scipy
- scanpy

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/gene-expression-pipeline.git
cd gene-expression-pipeline

# Set up the environment
bash rnaseq/scripts/setup_workspace.sh
```

## Usage

### Data Standardization

```bash
# Standardize ENTEx data
python rnaseq/scripts/standardize_datasets.py \
  --metadata /path/to/entex_metadata.json \
  --dataset-type entex \
  --output-dir /path/to/standardized_data
```

```bash
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py \
  --encode-dir /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines \
  --encode-entex-dir /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/entex \
  --gtex-file /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz \
  --mage-dir /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage \
  --adni-dir /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni_microarray \
  --output-dir /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data
```

### Creating Combined Datasets

Two types of combined datasets are now available:

#### 1. Combined Dataset with Common Genes Only (Original Approach)

This is the original approach that only includes genes found in all datasets:

```bash
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/create_combined_dataset.py
```

#### 2. Combined Dataset with ALL Genes (Sparse Matrix Approach)

This new approach includes all genes from all datasets using a sparse matrix representation:

```bash
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/create_combined_dataset_all_genes_sparse.py
```

The all-genes dataset has the following advantages:
- Preserves dataset-specific genes that might be important for certain analyses
- Provides a complete gene expression landscape across all datasets
- Uses sparse matrices to efficiently represent missing genes
- **Important**: Missing genes are properly represented as null/NaN values, not zeros
- Tracks which datasets each gene appears in through metadata

### Command-line Queries

```bash
# Query gene expression
python rnaseq/scripts/expression-query.py --symbol TP53 --tissue lung

# List available tissues
python rnaseq/scripts/expression-query.py --list-tissues

# Search for genes
python rnaseq/scripts/expression-query.py --search-gene BRCA
```

### Programmatic Access

```python
from rnaseq.scripts.optimized_gene_expression import (
    load_expression_data, get_gene_expression
)

# Load data (only needs to be done once)
loader = load_expression_data()

# Query gene expression
expr = get_gene_expression("entex", "TP53", tissue="lung")
print(f"TP53 expression in lung: Mean TPM = {expr['mean_expression']}")
```

### High-Performance Queries for Model Training

```python
from rnaseq.scripts.optimized_gene_expression import (
    get_gene_expression_matrix
)

# Get expression matrix for multiple genes across tissues
genes = ("TP53", "BRCA1", "BRCA2")
tissues = ("lung", "colon", "thyroid")
matrix, gene_ids, sample_ids = get_gene_expression_matrix("entex", genes, tissues=tissues)

# Use the expression matrix for model training
# ...
```

### Working with the All-Genes Combined Dataset

The sparse matrix version of the combined dataset properly distinguishes between:
- Genes that were measured with zero expression (actual zeros in the matrix)
- Genes that weren't measured in a dataset (null/NaN values)

This distinction is important for proper statistical analysis. Here's how to work with this dataset:

```python
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Load the all-genes combined dataset
adata = sc.read_h5ad('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_all_genes_standardized.h5ad')

# The data is stored in a sparse matrix
print(f"Matrix is sparse: {sp.issparse(adata.X)}")
print(f"Matrix shape: {adata.X.shape}")
print(f"Sparsity: {1.0 - (adata.X.count_nonzero() / (adata.X.shape[0] * adata.X.shape[1])):.2%}")

# See which datasets each gene appears in
gene_presence = adata.var['present_in_datasets']

# Find genes unique to specific datasets
gtex_specific_genes = [gene for gene in adata.var_names if 'gtex' in adata.var.loc[gene, 'present_in_datasets'] and 
                      adata.var.loc[gene, 'present_in_datasets'].count(',') == 0]

# Get all samples from a specific dataset
encode_samples = adata[adata.obs['source_dataset'] == 'encode', :]

# Working with sparse matrices - converting to dense only when needed
# For small operations, you can convert to dense
dense_subset = encode_samples[:, gtex_specific_genes[:10]].X.toarray()

# For large operations, keep sparse format and use sparse-compatible operations
mean_expr = encode_samples.X.mean(axis=0)

# Handle missing values appropriately (sparse matrices typically don't show NaNs)
def count_missing_genes(adata, dataset):
    """Count genes that are missing (not measured) in a specific dataset."""
    # Get the dataset's samples
    samples = adata[adata.obs['source_dataset'] == dataset, :]
    missing_count = 0
    
    # A gene is "missing" if it has no non-zero values across all samples in the dataset
    for gene_idx in range(samples.n_vars):
        gene_slice = samples.X[:, gene_idx]
        if gene_slice.count_nonzero() == 0:
            missing_count += 1
            
    return missing_count

# Example usage
missing_in_encode = count_missing_genes(adata, 'encode')
print(f"Genes not measured in ENCODE: {missing_in_encode}")
```

## ENTEx Data Integration

The ENTEx data integration pipeline includes:

1. Downloading gene quantification files
2. Generating standardized metadata
3. Processing and standardizing data
4. Verification of successful integration

To run the complete pipeline:

```bash
cd rnaseq/scripts
./entex_pipeline.sh
```

## Performance Optimization

The optimized gene expression module provides:

- Efficient indexing of genes, tissues, and donors
- LRU caching for repeated queries
- Hashmap-based lookups for O(1) access
- Memory-efficient sparse matrix handling
- Bulk query capabilities for matrix operations

## Technical Implementation Details

### Sparse Matrix Representation

The combined all-genes dataset uses sparse matrices to efficiently represent the expression data:

- **Missing genes** are represented as null/NaN values (not zeros)
- **Measured genes with zero expression** are represented as actual zeros
- **Non-zero expression values** are stored efficiently

This approach:
1. Maintains the biological distinction between "not measured" and "measured but not expressed"
2. Dramatically reduces memory requirements (often >90% reduction)
3. Works efficiently with scipy.sparse and anndata operations