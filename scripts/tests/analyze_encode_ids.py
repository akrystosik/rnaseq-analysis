#!/usr/bin/env python3
"""
ENCODE Original Gene ID Collector

This script collects all original gene IDs from ENCODE TSV files
and saves them for use in the preprocessing pipeline.
"""

import os
import glob
import pandas as pd
import numpy as np
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('encode_gene_id_collector')

def collect_gene_ids_from_file(file_path):
    """Extract gene IDs from a single TSV file."""
    try:
        df = pd.read_csv(file_path, sep='\t')
        
        if 'gene_id' not in df.columns:
            logger.warning(f"No gene_id column in {file_path}")
            return []
        
        # Get unique gene IDs
        gene_ids = df['gene_id'].unique().tolist()
        logger.info(f"Collected {len(gene_ids)} unique gene IDs from {os.path.basename(file_path)}")
        return gene_ids
    
    except Exception as e:
        logger.error(f"Error processing {file_path}: {e}")
        return []

def collect_all_gene_ids():
    """Collect gene IDs from all ENCODE TSV files."""
    encode_dir = Path('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data')
    
    # Find all TSV files
    tsv_files = []
    for root, dirs, files in os.walk(encode_dir):
        for file in files:
            if file.endswith('.tsv'):
                tsv_files.append(os.path.join(root, file))
    
    logger.info(f"Found {len(tsv_files)} TSV files in ENCODE directory")
    
    # Process each file and collect unique gene IDs
    all_gene_ids = set()
    for file in tsv_files:
        gene_ids = collect_gene_ids_from_file(file)
        all_gene_ids.update(gene_ids)
    
    # Convert to list
    gene_id_list = list(all_gene_ids)
    logger.info(f"Collected {len(gene_id_list)} unique gene IDs from all files")
    
    # Analyze gene ID types
    ensembl_ids = sum(1 for id in gene_id_list if str(id).startswith('ENSG'))
    numeric_ids = sum(1 for id in gene_id_list if str(id).isdigit())
    other_ids = len(gene_id_list) - ensembl_ids - numeric_ids
    
    logger.info(f"Gene ID analysis:")
    logger.info(f"  - Ensembl IDs: {ensembl_ids} ({ensembl_ids/len(gene_id_list)*100:.2f}%)")
    logger.info(f"  - Numeric IDs: {numeric_ids} ({numeric_ids/len(gene_id_list)*100:.2f}%)")
    logger.info(f"  - Other IDs: {other_ids} ({other_ids/len(gene_id_list)*100:.2f}%)")
    
    return gene_id_list

def save_gene_ids(gene_ids, output_file):
    """Save the collected gene IDs to a file."""
    # Create DataFrame with gene IDs
    df = pd.DataFrame({
        'gene_id': gene_ids,
        'is_ensembl': [str(id).startswith('ENSG') for id in gene_ids],
        'is_numeric': [str(id).isdigit() for id in gene_ids]
    })
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(gene_ids)} gene IDs to {output_file}")

def main():
    # Define output file path
    output_dir = Path('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping')
    output_file = output_dir / 'encode_original_gene_ids.csv'
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True)
    
    # Collect gene IDs
    gene_ids = collect_all_gene_ids()
    
    # Save gene IDs
    save_gene_ids(gene_ids, output_file)

if __name__ == "__main__":
    main()