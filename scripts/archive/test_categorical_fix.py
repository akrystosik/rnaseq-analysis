#!/usr/bin/env python3
"""
Test script for categorical data comparison fix
"""
import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('categorical_fix_test')

def fix_categorical_columns(adata):
    """
    Fix categorical columns to ensure they use the same categories.
    This prevents the "Categoricals can only be compared if 'categories' are the same" error.
    """
    logger.info("Checking for categorical columns that need fixing...")
    
    # Convert specified columns from categorical to string
    for col in adata.var.columns:
        if pd.api.types.is_categorical_dtype(adata.var[col]):
            logger.info(f"Converting categorical column {col} to string")
            adata.var[col] = adata.var[col].astype(str)
    
    logger.info("Categorical columns fixed")
    return adata

def test_with_example_data():
    """Create a test dataset with categorical columns that would cause comparison issues."""
    # Create a small AnnData object with categorical columns
    obs = pd.DataFrame({
        'sample_id': ['s1', 's2', 's3'],
        'group': pd.Categorical(['A', 'B', 'A'])
    })
    
    var = pd.DataFrame({
        'gene_id': pd.Categorical(['g1', 'g2', 'g3']),
        'gene_type': pd.Categorical(['protein_coding', 'lncRNA', 'protein_coding'])
    }, index=['g1', 'g2', 'g3'])
    
    X = np.random.rand(3, 3)
    
    adata = sc.AnnData(X=X, obs=obs, var=var)
    
    # Try a comparison that would fail without the fix
    try:
        logger.info("Testing comparison before fix...")
        result = adata.var['gene_id'] == 'g1'
        logger.info("Comparison succeeded unexpectedly")
    except TypeError as e:
        logger.info(f"Expected error occurred before fix: {e}")
    
    # Apply the fix
    adata = fix_categorical_columns(adata)
    
    # Try the comparison again
    try:
        logger.info("Testing comparison after fix...")
        result = adata.var['gene_id'] == 'g1'
        logger.info(f"Comparison succeeded after fix: {sum(result)} matches found")
        return True
    except TypeError as e:
        logger.error(f"Fix failed - error still occurs: {e}")
        return False

def test_with_real_data():
    """Test the fix with real data from a preprocessed dataset if available."""
    # Look for a preprocessed dataset file
    data_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/latest/preprocessed_data"
    
    if not os.path.exists(data_dir):
        logger.warning(f"Directory not found: {data_dir}")
        data_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data"
    
    # Look for any preprocessed h5ad file
    for file in os.listdir(data_dir):
        if file.endswith("_preprocessed.h5ad"):
            file_path = os.path.join(data_dir, file)
            logger.info(f"Found preprocessed file: {file_path}")
            
            try:
                # Load the dataset
                adata = sc.read_h5ad(file_path)
                logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
                
                # Check for categorical columns
                cat_columns = []
                for col in adata.var.columns:
                    if pd.api.types.is_categorical_dtype(adata.var[col]):
                        cat_columns.append(col)
                
                if not cat_columns:
                    logger.warning("No categorical columns found in var")
                    continue
                
                logger.info(f"Found categorical columns: {', '.join(cat_columns)}")
                
                # Try a comparison that might fail
                try:
                    logger.info(f"Testing comparison on {cat_columns[0]} before fix...")
                    test_value = adata.var[cat_columns[0]].iloc[0]
                    result = adata.var[cat_columns[0]] == test_value
                    logger.info("Comparison succeeded (no fix needed)")
                except TypeError as e:
                    logger.info(f"Expected error occurred before fix: {e}")
                    
                    # Apply the fix
                    adata = fix_categorical_columns(adata)
                    
                    # Try the comparison again
                    try:
                        logger.info(f"Testing comparison on {cat_columns[0]} after fix...")
                        result = adata.var[cat_columns[0]] == test_value
                        logger.info(f"Comparison succeeded after fix: {sum(result)} matches found")
                        return True
                    except TypeError as e:
                        logger.error(f"Fix failed - error still occurs: {e}")
                        return False
                
                return True
            
            except Exception as e:
                logger.error(f"Error loading or processing {file_path}: {e}")
    
    logger.warning("No suitable preprocessed file found for testing")
    return None

if __name__ == "__main__":
    logger.info("=== Testing Categorical Column Comparison Fix ===")
    
    # Test with example data
    logger.info("Testing with example data...")
    example_success = test_with_example_data()
    
    if example_success:
        logger.info("Example data test PASSED")
    else:
        logger.error("Example data test FAILED")
    
    # Test with real data if available
    logger.info("Testing with real data if available...")
    real_success = test_with_real_data()
    
    if real_success is True:
        logger.info("Real data test PASSED")
    elif real_success is False:
        logger.error("Real data test FAILED")
    else:
        logger.warning("Real data test SKIPPED (no suitable data found)")
    
    # Overall result
    if example_success:
        logger.info("Categorical comparison fix works!")
    else:
        logger.error("Categorical comparison fix needs more work")
