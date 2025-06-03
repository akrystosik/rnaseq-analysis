#!/usr/bin/env python3
"""
Review script for the current RNA-seq standardization pipeline
"""
import os
import re
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('pipeline_review')

# Define the base directory
BASE_DIR = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")

def review_script(script_name):
    """Review a specific script and identify key functions and issues."""
    file_path = os.path.join(SCRIPTS_DIR, script_name)
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return
    
    logger.info(f"Reviewing {script_name}...")
    
    # Read the file content
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Log file size and basic info
    logger.info(f"File size: {len(content):,} bytes")
    
    # Find all function definitions
    functions = re.findall(r'def ([a-zA-Z0-9_]+)\(', content)
    logger.info(f"Found {len(functions)} functions: {', '.join(functions[:10])}...")
    
    # Find potential issue points
    
    # 1. ADNI parsing - look for CSV reading
    adni_pattern = r'def process_adni_(.+?)def '
    adni_matches = re.findall(adni_pattern, content, re.DOTALL)
    if adni_matches:
        logger.info("Found ADNI processing code:")
        csv_read_pattern = r'pd\.read_csv\((.+?)\)'
        csv_reads = re.findall(csv_read_pattern, adni_matches[0])
        logger.info(f"  CSV reading operations: {len(csv_reads)}")
        for i, read in enumerate(csv_reads[:3]):
            logger.info(f"  Read {i+1}: {read.strip()}")
    
    # 2. Categorical comparisons
    categorical_pattern = r'is_categorical_dtype'
    categorical_matches = re.findall(categorical_pattern, content)
    logger.info(f"Found {len(categorical_matches)} categorical dtype checks")
    
    # 3. AnnData serialization
    save_pattern = r'def save_anndata(.+?)def '
    save_matches = re.findall(save_pattern, content, re.DOTALL)
    if save_matches:
        logger.info("Found save_anndata function:")
        # Look for serialization handling
        serialization_patterns = ['str(', 'astype(str)', 'to_dict(', 'tolist(']
        for pattern in serialization_patterns:
            pattern_count = save_matches[0].count(pattern)
            logger.info(f"  {pattern} operations: {pattern_count}")
    
    # 4. Placeholder IDs
    placeholder_pattern = r'PLACEHOLDER_'
    placeholder_matches = re.findall(placeholder_pattern, content)
    logger.info(f"Found {len(placeholder_matches)} references to placeholder IDs")
    
    logger.info(f"Review of {script_name} complete")

def review_main_scripts():
    """Review the main scripts in the pipeline."""
    main_scripts = [
        "standardize_datasets.py",
        "preprocess_dataset_gene_ids.py",
        "run_rnaseq_pipeline.sh"
    ]
    
    for script in main_scripts:
        review_script(script)

if __name__ == "__main__":
    logger.info("=== RNA-seq Pipeline Review ===")
    review_main_scripts()
