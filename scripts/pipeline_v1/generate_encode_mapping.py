#!/usr/bin/env python3
"""
ENCODE Gene ID Mapping Generator

This script processes ENCODE gene IDs from raw TSV files and creates a mapping
between original IDs and standard Ensembl IDs without version numbers.
"""

import os
import glob
import pandas as pd
import numpy as np
import logging
import argparse
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('encode_id_mapper')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate ENCODE gene ID mapping')
    parser.add_argument('--encode-dir', type=str, required=True,
                        help='Directory containing raw ENCODE data')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Directory to save mapping files')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration of mapping even if exists')
    
    return parser.parse_args()

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

def collect_all_gene_ids(encode_dir):
    """Collect gene IDs from all ENCODE TSV files."""
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
    
    return gene_id_list

def extract_ensembl_id(id_str):
    """Extract Ensembl gene ID from various formats."""
    id_str = str(id_str)
    
    # Direct Ensembl ID
    if id_str.startswith('ENSG'):
        return id_str.split('.')[0]  # Remove version number
    
    # Pipe-delimited composite entry
    if '|' in id_str:
        fields = id_str.split('|')
        if len(fields) > 1:
            # Look for ENSG in the second field
            ensembl_field = fields[1]
            if ensembl_field.startswith('ENSG'):
                return ensembl_field.split('.')[0]  # Remove version number
    
    # Numeric ID (potential Entrez)
    if id_str.isdigit():
        return f"ENTREZ:{id_str}"
    
    # Special case: spike-in controls
    if id_str.startswith('gSpikein'):
        return id_str
    
    # Other format - can't extract Ensembl ID
    return None

def process_gene_ids(gene_ids):
    """Process collected gene IDs and extract Ensembl IDs."""
    # Create DataFrame with gene IDs
    gene_ids_df = pd.DataFrame({'gene_id': gene_ids})
    
    # Apply the extraction function
    gene_ids_df['extracted_ensembl_id'] = gene_ids_df['gene_id'].apply(extract_ensembl_id)
    
    # Calculate statistics
    total_ids = len(gene_ids_df)
    extracted_count = gene_ids_df['extracted_ensembl_id'].notna().sum()
    extraction_rate = extracted_count / total_ids * 100
    
    direct_ensembl = sum(1 for id in gene_ids_df['gene_id'] if str(id).startswith('ENSG'))
    composite_with_ensembl = sum(1 for i, row in gene_ids_df.iterrows() 
                                if not str(row['gene_id']).startswith('ENSG') and 
                                row['extracted_ensembl_id'] is not None and
                                str(row['extracted_ensembl_id']).startswith('ENSG'))
    entrez_ids = sum(1 for id in gene_ids_df['extracted_ensembl_id'] if id is not None and str(id).startswith('ENTREZ:'))
    spikein_ids = sum(1 for id in gene_ids_df['extracted_ensembl_id'] if id is not None and str(id).startswith('gSpikein'))
    unmapped_ids = total_ids - extracted_count
    
    # Report statistics
    logger.info(f"Extraction statistics:")
    logger.info(f"  - Total IDs processed: {total_ids}")
    logger.info(f"  - Successfully extracted: {extracted_count} ({extraction_rate:.2f}%)")
    logger.info(f"  - Direct Ensembl IDs: {direct_ensembl} ({direct_ensembl/total_ids*100:.2f}%)")
    logger.info(f"  - Composite entries with Ensembl: {composite_with_ensembl} ({composite_with_ensembl/total_ids*100:.2f}%)")
    logger.info(f"  - Entrez IDs: {entrez_ids} ({entrez_ids/total_ids*100:.2f}%)")
    logger.info(f"  - Spike-in controls: {spikein_ids} ({spikein_ids/total_ids*100:.2f}%)")
    logger.info(f"  - Unmapped IDs: {unmapped_ids} ({unmapped_ids/total_ids*100:.2f}%)")
    
    return gene_ids_df

def create_mapping_file(processed_df, output_dir):
    """Create mapping files from processed gene IDs."""
    # Save processed IDs with details
    processed_file = os.path.join(output_dir, 'encode_processed_gene_ids.csv')
    processed_df.to_csv(processed_file, index=False)
    logger.info(f"Saved processed gene IDs to {processed_file}")
    
    # Create a mapping file for integration with pipeline
    mapping_df = processed_df[processed_df['extracted_ensembl_id'].notna()].copy()
    mapping_df = mapping_df[['gene_id', 'extracted_ensembl_id']].rename(
        columns={'gene_id': 'original_id', 'extracted_ensembl_id': 'ensembl_id'}
    )
    
    # Save mapping file
    mapping_file = os.path.join(output_dir, 'encode_id_to_ensembl_mapping.csv')
    mapping_df.to_csv(mapping_file, index=False)
    logger.info(f"Saved ID mapping to {mapping_file}")
    
    return mapping_file

def main():
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    mapping_file = os.path.join(args.output_dir, 'encode_id_to_ensembl_mapping.csv')
    if os.path.exists(mapping_file) and not args.force:
        logger.info(f"Mapping file {mapping_file} already exists. Use --force to regenerate.")
        return mapping_file
    
    # Collect gene IDs from ENCODE files
    gene_ids = collect_all_gene_ids(args.encode_dir)
    
    # Process gene IDs
    processed_df = process_gene_ids(gene_ids)
    
    # Create mapping file
    mapping_file = create_mapping_file(processed_df, args.output_dir)
    
    return mapping_file

if __name__ == "__main__":
    main()