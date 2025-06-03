#!/usr/bin/env python3
"""
Check Ensembl Gene ID Uniqueness in Standardized Datasets

This script checks whether Ensembl gene IDs are unique within each standardized dataset file
and reports any duplicates that are found.
"""

import os
import sys
import logging
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
from collections import Counter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gene_id_checker')

# Define paths
DATA_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data")

def check_dataset_gene_ids(file_path):
    """Check if gene IDs are unique within a dataset."""
    try:
        # Load dataset
        logger.info(f"Loading dataset from {file_path}")
        adata = sc.read_h5ad(file_path)
        
        # Get basic stats
        n_genes = adata.n_vars
        logger.info(f"Dataset has {n_genes} genes")
        
        # Check if var_names (gene IDs) are unique
        unique_genes = len(set(adata.var_names))
        logger.info(f"Number of unique gene IDs: {unique_genes}")
        
        if unique_genes < n_genes:
            # Find duplicates
            gene_counts = Counter(adata.var_names)
            duplicates = {gene: count for gene, count in gene_counts.items() if count > 1}
            
            logger.info(f"Found {len(duplicates)} genes with duplicates")
            logger.info(f"Duplicate examples:")
            
            # Show some examples of duplicated genes
            for i, (gene, count) in enumerate(sorted(duplicates.items(), key=lambda x: x[1], reverse=True)):
                if i >= 5:  # Show top 5 duplicates
                    break
                
                # Check if there are differences in the gene metadata for duplicates
                gene_indices = [i for i, name in enumerate(adata.var_names) if name == gene]
                metadata_diff = {}
                
                for field in adata.var.columns:
                    values = [adata.var[field].iloc[idx] for idx in gene_indices]
                    if len(set(values)) > 1:
                        metadata_diff[field] = values
                
                logger.info(f"  Gene {gene} appears {count} times")
                if metadata_diff:
                    logger.info(f"  Differences in metadata for {gene}:")
                    for field, values in metadata_diff.items():
                        logger.info(f"    {field}: {values}")
            
            return False, duplicates
        else:
            logger.info("All gene IDs are unique in this dataset")
            return True, {}
    
    except Exception as e:
        logger.error(f"Error analyzing dataset {file_path}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None, {}

def main():
    # Find all dataset files
    dataset_files = {
        'adni': DATA_DIR / 'adni_standardized_v2.h5ad',
        'encode': DATA_DIR / 'encode_standardized_v2.h5ad',
        'entex': DATA_DIR / 'entex_standardized_v2.h5ad',
        'mage': DATA_DIR / 'mage_standardized_v2.h5ad',
        'gtex': DATA_DIR / 'gtex_standardized_v2.h5ad',
        'combined': DATA_DIR / 'combined_standardized.h5ad',
        'combined_all_genes': DATA_DIR / 'combined_all_genes_standardized.h5ad'
    }
    
    # Check each dataset
    results = {}
    
    for name, file_path in dataset_files.items():
        if file_path.exists():
            logger.info(f"=== Checking {name} dataset ===")
            unique, duplicates = check_dataset_gene_ids(file_path)
            results[name] = {
                'all_unique': unique,
                'duplicate_count': len(duplicates) if duplicates else 0,
                'duplicates': duplicates
            }
            logger.info(f"=== {name} check complete ===\n")
        else:
            logger.warning(f"Dataset file not found: {file_path}")
    
    # Summarize results
    logger.info("=== Summary of Gene ID Uniqueness ===")
    for name, result in results.items():
        if result['all_unique']:
            logger.info(f"{name}: All gene IDs are unique")
        else:
            logger.info(f"{name}: Found {result['duplicate_count']} genes with duplicates")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())