#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/create_combined_dataset.py
#!/usr/bin/env python3
"""
Create Combined RNA-seq Dataset - Fixed Version

This script creates a combined AnnData object from individual standardized datasets,
excluding GTEx due to its large size.
"""

import os
import sys
import logging
import time
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('combine_datasets')

# Define paths
DATA_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data")
OUTPUT_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_standardized.h5ad")

def load_dataset(file_path):
    """Load a single dataset from an h5ad file."""
    logger.info(f"Loading dataset from {file_path}")
    start_time = time.time()
    try:
        adata = ad.read_h5ad(file_path)
        load_time = time.time() - start_time
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes in {load_time:.2f}s")
        return adata
    except Exception as e:
        logger.error(f"Error loading dataset: {e}")
        return None

def find_common_genes(datasets):
    """Find genes that are common to all datasets."""
    if not datasets:
        return set()
    
    # Start with genes from the first dataset
    common_genes = set(datasets[0].var_names)
    
    # Intersect with genes from each additional dataset
    for adata in datasets[1:]:
        common_genes = common_genes.intersection(set(adata.var_names))
    
    logger.info(f"Found {len(common_genes)} genes common to all datasets")
    return common_genes

def combine_datasets(datasets, dataset_names):
    """Combine multiple datasets into a single AnnData object."""
    if not datasets:
        logger.error("No datasets to combine")
        return None
    
    # Find common genes
    common_genes = find_common_genes(datasets)
    if not common_genes:
        logger.error("No common genes found across datasets")
        return None
    
    # Convert common genes set to sorted list for consistent ordering
    common_genes = sorted(list(common_genes))
    
    # For each dataset, filter to common genes and add dataset identifier
    filtered_datasets = []
    
    for adata, name in zip(datasets, dataset_names):
        # Subset to common genes
        logger.info(f"Filtering {name} to {len(common_genes)} common genes")
        adata_filtered = adata[:, common_genes].copy()
        
        # Add dataset identifier to obs
        adata_filtered.obs['source_dataset'] = name
        
        # Add dataset to list
        filtered_datasets.append(adata_filtered)
    
    # Concatenate datasets
    logger.info("Concatenating datasets")
    combined = ad.concat(
        filtered_datasets,
        join='inner',  # Inner join on genes (should already be filtered to common genes)
        label='source_dataset',  # Use dataset name as batch key
        keys=dataset_names,  # Dataset names as keys
        index_unique='-'  # Add hyphen to ensure unique sample IDs
    )
    
    # Merge var DataFrames to get comprehensive gene metadata
    logger.info("Merging gene metadata")
    
    # Use first dataset's var as base and convert categorical columns to string
    base_var = datasets[0].var.loc[common_genes].copy()
    
    # Convert categorical columns to string
    for col in base_var.columns:
        if pd.api.types.is_categorical_dtype(base_var[col]):
            logger.info(f"Converting categorical column {col} to string")
            base_var[col] = base_var[col].astype(str)
    
    # Merge metadata from other datasets
    for adata in datasets[1:]:
        var_df = adata.var.loc[common_genes].copy()
        
        # Go through each column and merge
        for col in var_df.columns:
            # Skip columns already in base_var
            if col in base_var.columns:
                # For categorical columns, convert to string first
                if pd.api.types.is_categorical_dtype(var_df[col]):
                    var_df[col] = var_df[col].astype(str)
                
                # Fill NAs in base_var with values from var_df
                base_var[col] = base_var[col].fillna(var_df[col])
            else:
                # For categorical columns, convert to string first
                if pd.api.types.is_categorical_dtype(var_df[col]):
                    var_df[col] = var_df[col].astype(str)
                
                # Add new column to base_var
                base_var[col] = var_df[col]
    
    # Replace var DataFrame in combined AnnData
    for col in base_var.columns:
        combined.var[col] = base_var[col]
    
    # Add additional metadata
    combined.uns['dataset_info'] = {
        'source': 'combined',
        'included_datasets': dataset_names,
        'sample_counts': {name: adata.n_obs for name, adata in zip(dataset_names, datasets)},
        'creation_date': pd.Timestamp.now().strftime('%Y-%m-%d'),
        'common_genes': len(common_genes)
    }
    
    return combined

def main():
    # Start timing
    start_time = time.time()
    
    # Find datasets
    dataset_files = {
        'adni': DATA_DIR / 'adni_standardized.h5ad',
        'encode': DATA_DIR / 'encode_standardized.h5ad',
        'entex': DATA_DIR / 'entex_standardized.h5ad',
        'mage': DATA_DIR / 'mage_standardized.h5ad'
    }
    
    # Load each dataset
    datasets = []
    dataset_names = []
    
    for name, file_path in dataset_files.items():
        if file_path.exists():
            adata = load_dataset(file_path)
            if adata is not None:
                datasets.append(adata)
                dataset_names.append(name)
        else:
            logger.warning(f"Dataset file not found: {file_path}")
    
    if not datasets:
        logger.error("No datasets loaded")
        return 1
    
    # Combine datasets
    logger.info(f"Combining {len(datasets)} datasets")
    combined = combine_datasets(datasets, dataset_names)
    
    if combined is None:
        logger.error("Failed to combine datasets")
        return 1
    
    # Save combined dataset
    logger.info(f"Saving combined dataset with {combined.n_obs} samples and {combined.n_vars} genes to {OUTPUT_FILE}")
    try:
        combined.write_h5ad(OUTPUT_FILE)
        logger.info("Successfully saved combined dataset")
    except Exception as e:
        logger.error(f"Error saving combined dataset: {e}")
        return 1
    
    # Log total processing time
    total_time = time.time() - start_time
    logger.info(f"Total processing time: {total_time:.2f} seconds")
    
    # Print dataset summary
    logger.info("Combined dataset summary:")
    logger.info(f"  Total samples: {combined.n_obs}")
    logger.info(f"  Common genes: {combined.n_vars}")
    logger.info("  Samples per dataset:")
    for name, count in combined.uns['dataset_info']['sample_counts'].items():
        logger.info(f"    {name}: {count}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())