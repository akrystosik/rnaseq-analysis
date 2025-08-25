#!/usr/bin/env python3
"""
Examine h5ad file structures to understand dataset schemas for Croissant metadata generation.

This script analyzes the structure and metadata of h5ad files to inform
the creation of Croissant JSON-LD metadata files.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import json
from pathlib import Path
import sys

def examine_h5ad_structure(file_path):
    """Examine the structure of an h5ad file and return metadata summary."""
    
    print(f"\n{'='*60}")
    print(f"Examining: {Path(file_path).name}")
    print(f"{'='*60}")
    
    try:
        # Load the h5ad file
        adata = sc.read_h5ad(file_path)
        
        # Basic structure
        structure_info = {
            'file_path': str(file_path),
            'shape': adata.shape,
            'n_vars': adata.n_vars,
            'n_obs': adata.n_obs,
            'obs_columns': adata.obs.columns.tolist(),
            'var_columns': adata.var.columns.tolist(),
            'uns_keys': list(adata.uns.keys()) if adata.uns else [],
            'layers': list(adata.layers.keys()) if adata.layers else [],
            'obsm_keys': list(adata.obsm.keys()) if adata.obsm else [],
            'varm_keys': list(adata.varm.keys()) if adata.varm else []
        }
        
        print(f"Shape: {adata.shape} (observations × variables)")
        print(f"Genes: {adata.n_vars:,}")
        print(f"Samples: {adata.n_obs:,}")
        
        print(f"\nObservation metadata columns ({len(adata.obs.columns)}):")
        for col in adata.obs.columns:
            unique_vals = adata.obs[col].nunique() if hasattr(adata.obs[col], 'nunique') else 'N/A'
            dtype = str(adata.obs[col].dtype)
            print(f"  - {col} ({dtype}): {unique_vals} unique values")
        
        print(f"\nVariable metadata columns ({len(adata.var.columns)}):")
        for col in adata.var.columns:
            unique_vals = adata.var[col].nunique() if hasattr(adata.var[col], 'nunique') else 'N/A'
            dtype = str(adata.var[col].dtype)
            print(f"  - {col} ({dtype}): {unique_vals} unique values")
        
        if adata.uns:
            print(f"\nUnstructured metadata (.uns) keys:")
            for key in adata.uns.keys():
                value_type = type(adata.uns[key]).__name__
                print(f"  - {key} ({value_type})")
        
        if adata.layers:
            print(f"\nData layers:")
            for layer in adata.layers.keys():
                print(f"  - {layer}")
        
        # Sample some key metadata for understanding
        print(f"\nSample observation metadata (first 3 rows):")
        display_cols = adata.obs.columns[:min(5, len(adata.obs.columns))]
        print(adata.obs[display_cols].head(3).to_string())
        
        # Check for specific important columns that would be relevant for Croissant
        important_columns = [
            'tissue_type', 'tissue_ontology_term_id', 'cell_type', 
            'disease', 'sex', 'age', 'ethnicity', 'ancestry_prediction',
            'dataset', 'assay', 'organism'
        ]
        
        print(f"\nImportant columns for Croissant metadata:")
        for col in important_columns:
            if col in adata.obs.columns:
                unique_count = adata.obs[col].nunique()
                sample_values = adata.obs[col].value_counts().head(3)
                print(f"  ✓ {col}: {unique_count} unique values")
                print(f"    Top values: {dict(sample_values)}")
            else:
                print(f"  ✗ {col}: not found")
        
        return structure_info
        
    except Exception as e:
        print(f"Error examining {file_path}: {str(e)}")
        return None

def main():
    """Main function to examine all h5ad files in current directory."""
    
    # Find all h5ad files
    current_dir = Path('.')
    h5ad_files = list(current_dir.glob('*.h5ad'))
    
    if not h5ad_files:
        print("No h5ad files found in current directory")
        return
    
    print(f"Found {len(h5ad_files)} h5ad files to examine")
    
    all_structures = {}
    
    for h5ad_file in sorted(h5ad_files):
        structure_info = examine_h5ad_structure(h5ad_file)
        if structure_info:
            dataset_name = h5ad_file.stem.replace('_standardized_preprocessed', '')
            all_structures[dataset_name] = structure_info
    
    # Save summary to JSON
    summary_file = 'h5ad_structure_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(all_structures, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"Structure analysis complete!")
    print(f"Summary saved to: {summary_file}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()