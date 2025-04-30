#!/usr/bin/env python3
"""
ENCODE Gene ID Processor

This script processes the original ENCODE gene IDs, extracting Ensembl IDs
from composite entries and preparing them for mapping.
"""

import pandas as pd
import numpy as np
import re
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('encode_id_processor')

# Load the collected gene IDs
gene_ids_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/encode_original_gene_ids.csv'
gene_ids_df = pd.read_csv(gene_ids_file)

logger.info(f"Processing {len(gene_ids_df)} ENCODE gene IDs")

# Function to extract Ensembl gene ID from composite entries
def extract_ensembl_id(id_str):
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

# Save processed IDs
output_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/encode_processed_gene_ids.csv'
gene_ids_df.to_csv(output_file, index=False)
logger.info(f"Saved processed gene IDs to {output_file}")

# Create a mapping file for integration with our pipeline
mapping_df = gene_ids_df[gene_ids_df['extracted_ensembl_id'].notna()].copy()
mapping_df = mapping_df[['gene_id', 'extracted_ensembl_id']].rename(
    columns={'gene_id': 'original_id', 'extracted_ensembl_id': 'ensembl_id'}
)

# Save mapping file
mapping_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/encode_id_to_ensembl_mapping.csv'
mapping_df.to_csv(mapping_file, index=False)
logger.info(f"Saved ID mapping to {mapping_file}")