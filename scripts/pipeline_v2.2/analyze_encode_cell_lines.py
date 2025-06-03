#!/usr/bin/env python3
"""
Analyze ENCODE cell line mapping from H5AD file
"""

import anndata as ad
import pandas as pd
import json

def analyze_encode_h5ad():
    """Analyze ENCODE H5AD file for cell line information"""
    
    # Load the ENCODE standardized H5AD file
    encode_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250528_235853/encode_standardized_v2.h5ad"
    
    try:
        adata = ad.read_h5ad(encode_file)
        
        print("=== ENCODE H5AD File Structure ===")
        print(f"Shape: {adata.shape}")
        print(f"Observations (samples): {adata.n_obs}")
        print(f"Variables (genes): {adata.n_vars}")
        
        print("\n=== Observation Metadata Columns ===")
        obs_columns = adata.obs.columns.tolist()
        for i, col in enumerate(obs_columns):
            print(f"{i+1:2d}. {col}")
        
        print("\n=== Sample Index (first 10) ===")
        print(adata.obs.index[:10].tolist())
        
        print("\n=== Sample of observation metadata ===")
        print(adata.obs.head(10))
        
        # Look for cell line related columns
        cell_line_cols = [col for col in obs_columns if 'cell' in col.lower() or 'line' in col.lower()]
        if cell_line_cols:
            print(f"\n=== Cell Line Related Columns ===")
            for col in cell_line_cols:
                print(f"Column: {col}")
                print(f"Unique values: {adata.obs[col].unique()[:20]}")  # First 20 unique values
                print()
        
        # Check for ENCFF patterns in index or any column
        print("\n=== ENCFF ID Analysis ===")
        encff_samples = [idx for idx in adata.obs.index if 'ENCFF' in str(idx)]
        print(f"Number of samples with ENCFF in index: {len(encff_samples)}")
        if encff_samples:
            print("Sample ENCFF IDs:", encff_samples[:10])
        
        # Look for any column that might contain ENCFF IDs
        for col in obs_columns:
            encff_in_col = adata.obs[col].astype(str).str.contains('ENCFF', na=False).sum()
            if encff_in_col > 0:
                print(f"Column '{col}' contains {encff_in_col} ENCFF entries")
        
        # Look for common cell line names
        common_cell_lines = ['A549', 'K562', 'HepG2', 'GM12878', 'MCF-7', 'H1-hESC', 'HeLa']
        print(f"\n=== Common Cell Line Detection ===")
        for cell_line in common_cell_lines:
            for col in obs_columns:
                matches = adata.obs[col].astype(str).str.contains(cell_line, case=False, na=False).sum()
                if matches > 0:
                    print(f"Cell line '{cell_line}' found {matches} times in column '{col}'")
        
        return adata.obs
        
    except Exception as e:
        print(f"Error loading H5AD file: {e}")
        return None

def analyze_encode_metadata():
    """Analyze ENCODE metadata JSON file"""
    
    metadata_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/encode_metadata.json"
    
    try:
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        
        print("\n" + "="*50)
        print("=== ENCODE Metadata JSON Analysis ===")
        print("="*50)
        
        print(f"Top-level keys: {list(metadata.keys())}")
        
        # Analyze structure
        for key, value in metadata.items():
            print(f"\nKey: {key}")
            print(f"Type: {type(value)}")
            if isinstance(value, dict):
                print(f"Sub-keys: {list(value.keys())[:10]}")  # First 10 keys
                if len(value) <= 20:  # If small enough, show all
                    for subkey, subval in value.items():
                        print(f"  {subkey}: {subval}")
            elif isinstance(value, list):
                print(f"List length: {len(value)}")
                if value and len(value) <= 10:  # Show small lists
                    print(f"Items: {value}")
        
        return metadata
        
    except Exception as e:
        print(f"Error loading metadata JSON: {e}")
        return None

if __name__ == "__main__":
    print("Analyzing ENCODE cell line mapping...")
    
    # Analyze H5AD file
    obs_data = analyze_encode_h5ad()
    
    # Analyze metadata JSON
    metadata = analyze_encode_metadata()
    
    print("\n" + "="*50)
    print("Analysis complete.")