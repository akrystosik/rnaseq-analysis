#!/usr/bin/env python3
"""
Fix Categorical Columns in AnnData Objects

This script fixes categorical column comparisons by converting categorical columns
to string types to ensure proper comparisons.
"""

import os
import sys
import logging
import argparse
import scanpy as sc
import pandas as pd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('categorical_fix')

def fix_categorical_columns(adata, columns=None):
    """
    Convert categorical columns to string type to avoid comparison issues.
    
    Args:
        adata: AnnData object
        columns: List of columns to convert, or None to check all columns
    
    Returns:
        Modified AnnData object
    """
    # If no columns specified, check all var columns
    if columns is None:
        columns = adata.var.columns
    
    # Track modified columns
    modified_columns = []
    
    for col in columns:
        if col in adata.var.columns and isinstance(adata.var[col].dtype, pd.CategoricalDtype):
            logger.info(f"Converting categorical column {col} to string")
            adata.var[col] = adata.var[col].astype(str)
            modified_columns.append(col)
    
    if modified_columns:
        logger.info(f"Converted {len(modified_columns)} columns to string: {', '.join(modified_columns)}")
    else:
        logger.info("No categorical columns needed conversion")
    
    return adata

def process_file(input_file, output_file=None):
    """
    Process a single AnnData file to fix categorical columns.
    
    Args:
        input_file: Path to input AnnData file
        output_file: Path to output AnnData file (optional)
    
    Returns:
        Path to the output file
    """
    if output_file is None:
        output_file = input_file.replace('.h5ad', '_fixed_categorical.h5ad')
    
    logger.info(f"Processing {input_file}")
    
    # Load AnnData
    adata = sc.read_h5ad(input_file)
    logger.info(f"Loaded AnnData with {adata.n_obs} observations and {adata.n_vars} variables")
    
    # Fix categorical columns
    adata = fix_categorical_columns(adata)
    
    # Save modified AnnData
    logger.info(f"Saving to {output_file}")
    adata.write_h5ad(output_file)
    logger.info(f"Saved fixed AnnData to {output_file}")
    
    return output_file

def process_directory(input_dir, output_dir=None, file_pattern='*.h5ad'):
    """
    Process all AnnData files in a directory.
    
    Args:
        input_dir: Path to input directory
        output_dir: Path to output directory (optional)
        file_pattern: Pattern to match AnnData files
    
    Returns:
        Number of files processed
    """
    if output_dir is None:
        output_dir = input_dir
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all matching files
    import glob
    h5ad_files = glob.glob(os.path.join(input_dir, file_pattern))
    
    if not h5ad_files:
        logger.warning(f"No files found matching {file_pattern} in {input_dir}")
        return 0
    
    logger.info(f"Found {len(h5ad_files)} files to process")
    
    # Process each file
    processed_files = 0
    for input_file in h5ad_files:
        output_file = os.path.join(output_dir, os.path.basename(input_file))
        try:
            process_file(input_file, output_file)
            processed_files += 1
        except Exception as e:
            logger.error(f"Error processing {input_file}: {e}")
            import traceback
            logger.error(traceback.format_exc())
    
    logger.info(f"Processed {processed_files} files")
    return processed_files

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Fix categorical columns in AnnData objects')
    parser.add_argument('--input', required=True, help='Input AnnData file or directory')
    parser.add_argument('--output', help='Output AnnData file or directory')
    parser.add_argument('--pattern', default='*.h5ad', help='File pattern when processing directories')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Check if input is a file or directory
    if os.path.isfile(args.input):
        # Process single file
        process_file(args.input, args.output)
    elif os.path.isdir(args.input):
        # Process directory
        process_directory(args.input, args.output, args.pattern)
    else:
        logger.error(f"Input not found: {args.input}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
