#!/usr/bin/env python3
"""
Gene ID Mapping Analysis Script

This script analyzes the results of gene ID preprocessing to check
how well the mapping performed and whether version numbers are preserved.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze gene ID mapping results')
    parser.add_argument('--input-dir', type=str, required=True,
                        help='Directory containing preprocessed datasets')
    
    return parser.parse_args()

def analyze_dataset(file_path, dataset_name):
    """Analyze gene ID mapping for a single dataset."""
    print(f"\n=== Analyzing {dataset_name} dataset ===")
    
    # Load dataset
    adata = sc.read_h5ad(file_path)
    print(f"Dataset shape: {adata.shape[0]} samples × {adata.shape[1]} genes")
    
    # Check if required columns exist
    required_cols = ['gene_id', 'original_gene_id', 'ensembl_id']
    missing_cols = [col for col in required_cols if col not in adata.var.columns]
    
    if missing_cols:
        print(f"WARNING: Missing columns: {', '.join(missing_cols)}")
        return
    
    # Check preserved version numbers
    ensembl_with_version = sum(1 for gene_id in adata.var['original_gene_id'] 
                             if str(gene_id).startswith('ENSG') and '.' in str(gene_id))
    ensembl_without_version = sum(1 for gene_id in adata.var['original_gene_id'] 
                                if str(gene_id).startswith('ENSG') and '.' not in str(gene_id))
    
    print(f"Original Ensembl IDs with version numbers: {ensembl_with_version}")
    print(f"Original Ensembl IDs without version numbers: {ensembl_without_version}")
    
    # Check mapping statistics
    mapped_genes = sum(1 for gene in adata.var['ensembl_id'] if str(gene).strip() != '')
    mapping_percentage = mapped_genes / adata.shape[1] * 100
    
    print(f"Mapped genes: {mapped_genes}/{adata.shape[1]} ({mapping_percentage:.2f}%)")
    
    # Check mapping sources
    if 'mapping_source' in adata.var.columns:
        source_counts = adata.var['mapping_source'].value_counts()
        print("\nMapping sources:")
        for source, count in source_counts.items():
            source_percentage = count / adata.shape[1] * 100
            print(f"  - {source}: {count} ({source_percentage:.2f}%)")
    
    # Check for placeholders
    placeholders = sum(1 for gene in adata.var['ensembl_id'] 
                     if str(gene).startswith('PLACEHOLDER_'))
    if placeholders > 0:
        placeholder_percentage = placeholders / adata.shape[1] * 100
        print(f"\nPlaceholder IDs: {placeholders}/{adata.shape[1]} ({placeholder_percentage:.2f}%)")
        
        # Show sample placeholders
        placeholder_sample = adata.var[adata.var['ensembl_id'].astype(str).str.startswith('PLACEHOLDER_')].head(5)
        print("\nSample placeholder entries:")
        print(placeholder_sample[['gene_id', 'original_gene_id', 'ensembl_id']])

def main():
    """Main function to analyze all preprocessed datasets."""
    args = parse_arguments()
    
    # Find all preprocessed datasets
    datasets = []
    for file in os.listdir(args.input_dir):
        if file.endswith("_preprocessed.h5ad"):
            dataset_name = file.split("_")[0]
            file_path = os.path.join(args.input_dir, file)
            datasets.append((dataset_name, file_path))
    
    print(f"Found {len(datasets)} preprocessed datasets in {args.input_dir}")
    
    # Analyze each dataset
    for dataset_name, file_path in datasets:
        analyze_dataset(file_path, dataset_name)

if __name__ == "__main__":
    main()