#!/usr/bin/env python3
"""
Wrapper script to save AnnData objects with proper string conversion.
Usage: python anndata_save_wrapper.py input_file output_file
"""

import sys
import os
import scanpy as sc
import pandas as pd
import numpy as np
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('anndata_save_wrapper')

def convert_to_serializable(obj):
    """Convert dict values to serializable types."""
    if isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif pd.isna(obj):
        return ""
    elif isinstance(obj, (pd.Series, pd.DataFrame, np.ndarray)):
        if hasattr(obj, "tolist"):
            return obj.tolist()
        else:
            return str(obj)
    else:
        return str(obj)

def save_adata_safely(input_file, output_file):
    """Load AnnData and save it with proper string conversion."""
    try:
        logger.info(f"Loading AnnData from {input_file}")
        adata = sc.read_h5ad(input_file)
        
        logger.info(f"Loaded AnnData with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Make a copy of the original uns
        original_uns = adata.uns.copy()
        
        # Convert all uns values to serializable types
        logger.info("Converting uns dictionary to serializable types")
        adata.uns = convert_to_serializable(original_uns)
        
        # Save the AnnData object
        logger.info(f"Saving AnnData to {output_file}")
        adata.write(output_file)
        
        # Verify the save worked
        logger.info("Verifying saved file")
        test_adata = sc.read_h5ad(output_file)
        logger.info(f"Verification successful: {test_adata.n_obs} samples, {test_adata.n_vars} genes")
        
        return True
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python anndata_save_wrapper.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    success = save_adata_safely(input_file, output_file)
    
    if success:
        print(f"Successfully saved AnnData to {output_file}")
        sys.exit(0)
    else:
        print("Failed to save AnnData")
        sys.exit(1)
