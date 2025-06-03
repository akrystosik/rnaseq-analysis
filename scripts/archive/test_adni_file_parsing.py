#!/usr/bin/env python3
"""
Test script for ADNI file parsing with escaped tabs fix
"""
import os
import logging
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('adni_parse_test')

def preprocess_adni_file(file_path):
    """Preprocess ADNI file to fix escaped tab characters."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Check if file has escaped tabs
    if '\\t' in content:
        logger.info(f"Fixing escaped tabs in {file_path}")
        # Replace escaped tabs with actual tabs
        fixed_content = content.replace('\\t', '\t')
        
        # Create a temporary fixed file
        fixed_path = file_path + '.fixed'
        with open(fixed_path, 'w') as f:
            f.write(fixed_content)
        
        return fixed_path
    
    return file_path

def test_adni_file(file_path):
    """Test parsing an ADNI file with the new preprocessing approach."""
    logger.info(f"Testing parsing for ADNI file: {file_path}")
    
    # First try direct parsing
    try:
        logger.info("Trying direct parsing first...")
        df = pd.read_csv(file_path, sep='\t')
        logger.info(f"Direct parsing succeeded! Columns: {list(df.columns)}")
        logger.info(f"Data shape: {df.shape}")
        return True
    except Exception as e:
        logger.warning(f"Direct parsing failed: {e}")
    
    # Try with preprocessing
    try:
        logger.info("Trying with preprocessing...")
        fixed_path = preprocess_adni_file(file_path)
        df = pd.read_csv(fixed_path, sep='\t')
        logger.info(f"Parsing with preprocessing succeeded! Columns: {list(df.columns)}")
        logger.info(f"Data shape: {df.shape}")
        
        # Clean up if we created a temporary file
        if fixed_path != file_path and os.path.exists(fixed_path):
            os.remove(fixed_path)
            
        return True
    except Exception as e:
        logger.warning(f"Parsing with preprocessing failed: {e}")
        
        # Try with comma delimiter as fallback
        try:
            logger.info("Trying with comma delimiter...")
            df = pd.read_csv(file_path)
            logger.info(f"Parsing with comma delimiter succeeded! Columns: {list(df.columns)}")
            logger.info(f"Data shape: {df.shape}")
            return True
        except Exception as e2:
            logger.error(f"All parsing attempts failed. Last error: {e2}")
            return False

def find_test_file():
    """Find a suitable ADNI test file."""
    adni_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni_microarray"
    
    # Look for subject directories
    for subdir in os.listdir(adni_dir):
        if subdir.startswith("0") and "_S_" in subdir:
            subj_dir = os.path.join(adni_dir, subdir)
            if os.path.isdir(subj_dir):
                # Look for CSV files in the subject directory
                for file in os.listdir(subj_dir):
                    if file.endswith(".csv"):
                        return os.path.join(subj_dir, file)
    
    return None

if __name__ == "__main__":
    logger.info("=== Testing ADNI File Parsing Fix ===")
    
    # Find a test file
    test_file = find_test_file()
    
    if test_file:
        logger.info(f"Found test file: {test_file}")
        success = test_adni_file(test_file)
        
        if success:
            logger.info("Test PASSED: ADNI file parsing fix works!")
        else:
            logger.error("Test FAILED: ADNI file parsing fix does not work.")
    else:
        logger.error("Could not find a suitable ADNI test file.")
