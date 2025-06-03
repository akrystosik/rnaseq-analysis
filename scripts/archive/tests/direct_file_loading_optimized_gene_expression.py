#!/usr/bin/env python3
"""
Simplified test for direct file loading in the optimized_gene_expression module
"""

import sys
import os

# Add the parent directory to the Python path so we can import the module
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)

# Import just the loader function
from optimized_gene_expression import ExpressionDataLoader

def test_direct_file_loading():
    print("=== Testing Direct File Loading ===")
    
    # Define the file path for the sparse combined file
    combined_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_all_genes_standardized.h5ad"
    
    print(f"Loading file: {combined_file}")
    
    # Create a loader instance
    loader = ExpressionDataLoader()
    
    # Use the load_file_directly method
    adata = loader.load_file_directly("combined_direct", combined_file)
    
    if adata is not None:
        print(f"Successfully loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Check dataset distribution if available
        if 'dataset' in adata.obs.columns:
            dataset_counts = adata.obs['dataset'].value_counts()
            print("\nDataset distribution:")
            for dataset, count in dataset_counts.items():
                percent = (count / adata.n_obs) * 100
                print(f"  {dataset}: {count} samples ({percent:.2f}%)")
        
        # Look at some basic gene information
        if 'gene_name' in adata.var.columns:
            print("\nSample genes in the dataset:")
            # Get 5 random gene names
            import random
            sample_indices = random.sample(range(adata.n_vars), 5)
            for idx in sample_indices:
                gene_id = adata.var_names[idx]
                gene_name = adata.var['gene_name'][idx] if 'gene_name' in adata.var.columns else "Unknown"
                print(f"  {gene_id}: {gene_name}")
        
        return True
    else:
        print("Failed to load the dataset directly")
        return False

if __name__ == "__main__":
    success = test_direct_file_loading()
    print(f"\nTest {'succeeded' if success else 'failed'}")