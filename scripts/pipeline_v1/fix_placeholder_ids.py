#!/usr/bin/env python3
"""
Fix placeholder IDs in preprocessed datasets by converting them to proper Entrez IDs.
"""
import os
import sys
import re
import logging
import scanpy as sc
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('placeholder_id_fixer')

def fix_placeholder_ids(h5ad_file, output_file):
    """
    Fix placeholder IDs in the specified h5ad file.
    
    Parameters:
    -----------
    h5ad_file : str
        Path to the h5ad file to fix
    output_file : str
        Path to save the fixed h5ad file
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Load the dataset
        logger.info(f"Loading dataset from {h5ad_file}")
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Check if there are placeholder IDs
        placeholder_pattern = re.compile(r'^PLACEHOLDER_(.+)$')
        
        # Find placeholder IDs in gene_id and ensembl_id columns
        placeholders = []
        for col in ['gene_id', 'ensembl_id']:
            if col in adata.var.columns:
                # Count placeholders
                placeholder_count = sum(1 for g in adata.var[col] if isinstance(g, str) and g.startswith('PLACEHOLDER_'))
                if placeholder_count > 0:
                    logger.info(f"Found {placeholder_count} placeholder IDs in {col} column")
                    placeholders.append(col)
        
        if not placeholders:
            logger.info("No placeholder IDs found, no fix needed")
            return True
        
        # Make a copy of var DataFrame
        var_df = adata.var.copy()
        
        # Fix categorical columns to prevent comparison issues
        for col in var_df.columns:
            if pd.api.types.is_categorical_dtype(var_df[col]):
                logger.info(f"Converting categorical column {col} to string")
                var_df[col] = var_df[col].astype(str)
        
        # Replace placeholder IDs with Entrez IDs
        entrez_count = 0
        other_count = 0
        
        for col in placeholders:
            for idx in var_df.index:
                value = var_df.loc[idx, col]
                if isinstance(value, str) and value.startswith('PLACEHOLDER_'):
                    # Extract the ID from the placeholder
                    match = placeholder_pattern.match(value)
                    if match:
                        id_part = match.group(1)
                        
                        # If it's a numeric ID, use Entrez prefix
                        if id_part.isdigit():
                            var_df.loc[idx, col] = f"ENTREZ:{id_part}"
                            entrez_count += 1
                        else:
                            # For non-numeric IDs, keep as is but remove PLACEHOLDER_
                            var_df.loc[idx, col] = id_part
                            other_count += 1
        
        # Replace the var DataFrame
        adata.var = var_df
        
        # Log results
        logger.info(f"Fix results: {entrez_count + other_count} placeholder IDs converted")
        logger.info(f"Total Entrez IDs: {entrez_count}")
        logger.info(f"Total other IDs: {other_count}")
        
        # Save the fixed dataset
        logger.info(f"Saving fixed dataset to {output_file}")
        adata.write_h5ad(output_file)
        
        logger.info("Fix completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error fixing placeholder IDs: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fix_placeholder_ids.py <input_h5ad_file> [output_h5ad_file]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        # Use same file with .fixed suffix
        output_file = input_file + ".fixed"
    
    success = fix_placeholder_ids(input_file, output_file)
    
    if success:
        print(f"Successfully fixed placeholder IDs and saved to {output_file}")
        sys.exit(0)
    else:
        print(f"Failed to fix placeholder IDs")
        sys.exit(1)
