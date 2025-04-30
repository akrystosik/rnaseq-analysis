#!/usr/bin/env python3
"""
Fix Placeholder Gene IDs

This script fixes placeholder gene IDs in preprocessed datasets by converting
them to proper Entrez identifiers, properly handling categorical data types.

Usage:
    python fix_placeholder_ids.py input_file output_file
"""

import sys
import scanpy as sc
import pandas as pd
import numpy as np
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('placeholder_id_fixer')

def fix_placeholder_ids(input_file, output_file):
    """Fix placeholder gene IDs in a preprocessed dataset."""
    try:
        logger.info(f"Loading dataset from {input_file}")
        adata = sc.read_h5ad(input_file)
        
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Count placeholder IDs before fix
        placeholder_count_before = sum(1 for id in adata.var['ensembl_id'] 
                                      if str(id).startswith('PLACEHOLDER_'))
        
        if placeholder_count_before > 0:
            logger.info(f"Found {placeholder_count_before} placeholder IDs to fix")
            
            # First handle categorical columns by converting to string
            var_df = adata.var.copy()
            
            # Convert categorical columns to strings
            for col in var_df.columns:
                if pd.api.types.is_categorical_dtype(var_df[col]):
                    logger.info(f"Converting categorical column {col} to string")
                    var_df[col] = var_df[col].astype(str)
            
            # Create a mapping from placeholder IDs to Entrez IDs
            placeholder_to_entrez = {}
            placeholder_indices = []
            
            # Identify placeholder IDs and create mapping
            for idx, row in var_df.iterrows():
                ensembl_id = row['ensembl_id']
                if str(ensembl_id).startswith('PLACEHOLDER_'):
                    # Extract the numeric part from the placeholder ID
                    placeholder_num = str(ensembl_id).replace('PLACEHOLDER_', '')
                    # Use this as the Entrez ID with ENTREZ: prefix
                    entrez_id = f"ENTREZ:{placeholder_num}"
                    placeholder_to_entrez[ensembl_id] = entrez_id
                    placeholder_indices.append(idx)
            
            # Fix the placeholder IDs
            logger.info(f"Replacing {len(placeholder_indices)} placeholder IDs with Entrez IDs")
            for idx in placeholder_indices:
                old_id = var_df.at[idx, 'ensembl_id']
                var_df.at[idx, 'ensembl_id'] = placeholder_to_entrez[old_id]
                var_df.at[idx, 'mapping_source'] = 'entrez_id'
            
            # Replace adata.var with the fixed DataFrame
            adata.var = var_df
            
            # Count IDs after fix
            entrez_count = sum(1 for id in adata.var['ensembl_id'] 
                              if str(id).startswith('ENTREZ:'))
            
            logger.info(f"Fix results: {placeholder_count_before} placeholder IDs converted to Entrez IDs")
            logger.info(f"Total Entrez IDs: {entrez_count}")
            
            # Save the fixed dataset
            logger.info(f"Saving fixed dataset to {output_file}")
            adata.write(output_file)
            
            logger.info("Fix completed successfully")
            return True
        else:
            logger.info("No placeholder IDs found, no fix needed")
            # Save the dataset as is
            adata.write(output_file)
            return True
    
    except Exception as e:
        logger.error(f"Error fixing placeholder IDs: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_placeholder_ids.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    success = fix_placeholder_ids(input_file, output_file)
    
    if success:
        print(f"Successfully fixed placeholder IDs and saved to {output_file}")
        sys.exit(0)
    else:
        print("Failed to fix placeholder IDs")
        sys.exit(1)
