#!/usr/bin/env python3
"""
Test ADNI file parsing with fixed tab delimiter handling
"""
import os
import pandas as pd
import logging
import sys

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('adni_parser_test')

# Import the necessary functions from standardize_datasets.py
sys.path.append('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/')
from standardize_datasets import preprocess_adni_file, identify_columns

def test_adni_file(file_path):
    """Test preprocessing and parsing an ADNI file"""
    logger.info(f"Testing file: {file_path}")
    
    # Test the preprocessor
    try:
        fixed_path = preprocess_adni_file(file_path)
        logger.info(f"Preprocessed path: {fixed_path}")
        
        # Try to read with tab delimiter
        df = pd.read_csv(fixed_path, sep='\t')
        logger.info(f"SUCCESS: Read with tab delimiter. Shape: {df.shape}")
        
        # Identify columns
        gene_id_col, tpm_col = identify_columns(df)
        logger.info(f"Identified columns: gene_id={gene_id_col}, tpm={tpm_col}")
        
        # Check values
        if gene_id_col and tpm_col:
            logger.info(f"First few genes and values:")
            for i, (gene, val) in enumerate(zip(df[gene_id_col], df[tpm_col])):
                logger.info(f"  {gene}: {val}")
                if i >= 4:  # Just show first 5
                    break
            return True
        else:
            logger.error("Failed to identify columns")
            return False
            
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

def main():
    # Test with a few specific files
    adni_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni_microarray'
    test_files = []
    
    # Find a few subject directories
    subject_dirs = []
    for item in os.listdir(adni_dir):
        if os.path.isdir(os.path.join(adni_dir, item)) and '_S_' in item:
            subject_dirs.append(os.path.join(adni_dir, item))
            if len(subject_dirs) >= 3:  # Test with 3 subjects
                break
    
    # Get one file from each subject directory
    for subject_dir in subject_dirs:
        for file in os.listdir(subject_dir):
            if file.endswith('.csv'):
                test_files.append(os.path.join(subject_dir, file))
                break  # Just get one file per directory
    
    # Run the tests
    success_count = 0
    for file_path in test_files:
        if test_adni_file(file_path):
            success_count += 1
    
    logger.info(f"Test summary: {success_count}/{len(test_files)} files processed successfully")
    return 0 if success_count == len(test_files) else 1

if __name__ == "__main__":
    sys.exit(main())
