#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis/scripts/entrez-to-ensembl-mapping.py
#!/usr/bin/env python3
"""
Entrez to Ensembl Mapping Creator

This script creates a mapping between Entrez Gene IDs and Ensembl IDs
using NCBI's gene2ensembl file, which is the definitive source for this mapping.

Usage:
  python entrez-to-ensembl-mapping.py --output /path/to/output.csv
"""

import os
import sys
import argparse
import pandas as pd
import requests
import gzip
import io
import logging
import time
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('entrez_to_ensembl')

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis")
DEFAULT_OUTPUT = BASE_DIR / "metadata/entrez_to_ensembl_mapping.csv"

def create_gene_mapping(output_file=DEFAULT_OUTPUT):
    """
    Create a mapping file between Entrez Gene IDs and Ensembl IDs.
    
    Parameters:
    -----------
    output_file : str
        Path to save the mapping file
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Download NCBI gene2ensembl file
        # This is the official mapping between Entrez Gene IDs and Ensembl
        logger.info("Downloading gene2ensembl from NCBI...")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
        
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        # Prepare to store the mapping data
        mapping_data = []
        total_lines = 0
        human_count = 0
        
        # Parse the file
        logger.info("Parsing gene2ensembl file...")
        with gzip.open(io.BytesIO(response.content)) as f:
            # Skip header
            next(f)
            
            # Process each line
            for line in f:
                total_lines += 1
                line = line.decode('utf-8').strip()
                fields = line.split('\t')
                
                # Check if this is a human gene (tax_id = 9606)
                if fields[0] == '9606':
                    entrez_id = fields[1]
                    ensembl_gene = fields[2]
                    
                    # Store the mapping
                    mapping_data.append({
                        'entrez_id': entrez_id,
                        'ensembl_id': ensembl_gene
                    })
                    
                    human_count += 1
                
                # Show progress
                if total_lines % 100000 == 0:
                    logger.info(f"Processed {total_lines:,} lines, found {human_count:,} human gene mappings")
        
        # Create a DataFrame and save to CSV
        mapping_df = pd.DataFrame(mapping_data)
        mapping_df.to_csv(output_file, index=False)
        
        logger.info(f"Created mapping file with {len(mapping_df):,} entries")
        logger.info(f"Saved mapping to {output_file}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error creating gene mapping: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Entrez to Ensembl Mapping Creator')
    parser.add_argument('--output', default=DEFAULT_OUTPUT, help='Path to save the mapping file')
    
    args = parser.parse_args()
    
    start_time = time.time()
    success = create_gene_mapping(args.output)
    
    if success:
        logger.info(f"Successfully created mapping in {time.time() - start_time:.2f} seconds")
    else:
        logger.error("Failed to create mapping")
        sys.exit(1)

if __name__ == '__main__':
    main()