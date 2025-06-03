#!/usr/bin/env python3
"""
Entrez to Ensembl Mapping Creator

This script creates a mapping between Entrez Gene IDs and Ensembl IDs
using NCBI's gene2ensembl file, which is the definitive source for this mapping.

Usage:
  python entrez-to-ensembl-mapping.py --output /path/to/output.csv [--species human]
"""

import os
import sys
import argparse
import pandas as pd
import urllib.request
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
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
DEFAULT_OUTPUT = BASE_DIR / "metadata/json/entrez_to_ensembl_mapping.csv"

def get_taxon_id(species):
    """Get NCBI taxonomy ID for a species."""
    species_tax_map = {
        'human': '9606',
        'mouse': '10090',
        'rat': '10116',
        'fly': '7227',
        'worm': '6239'
    }
    return species_tax_map.get(species.lower())

def create_gene_mapping(output_file=DEFAULT_OUTPUT, species="human"):
    """
    Create a mapping file between Entrez Gene IDs and Ensembl IDs.
    
    Parameters:
    -----------
    output_file : str
        Path to save the mapping file
    species : str
        Species to generate mapping for ('human', 'mouse', etc.)
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Get taxon ID for species
        taxon_id = get_taxon_id(species)
        if not taxon_id:
            logger.error(f"Invalid species: {species}")
            return False
            
        logger.info(f"Creating mapping for {species} (taxon ID {taxon_id})")
        
        # Download NCBI gene2ensembl file directly to memory
        logger.info("Downloading gene2ensembl from NCBI...")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
        
        response = urllib.request.urlopen(url)
        compressed_file = io.BytesIO(response.read())
        
        # Prepare to store the mapping data
        mapping_data = []
        total_lines = 0
        species_count = 0
        
        # Parse the file
        logger.info("Parsing gene2ensembl file...")
        with gzip.open(compressed_file, 'rt') as f:
            # Skip header
            next(f)
            
            # Process each line
            for line in f:
                total_lines += 1
                fields = line.strip().split('\t')
                
                # Check if this is the species we want
                if fields[0] == taxon_id:
                    entrez_id = fields[1]
                    ensembl_gene = fields[2]
                    
                    # Store the mapping
                    mapping_data.append({
                        'entrez_id': entrez_id,
                        'ensembl_id': ensembl_gene
                    })
                    
                    species_count += 1
                
                # Show progress
                if total_lines % 1000000 == 0:
                    logger.info(f"Processed {total_lines:,} lines, found {species_count:,} {species} gene mappings")
        
        # Create a DataFrame directly from the mapping data
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
    parser.add_argument('--species', default='human', choices=['human', 'mouse', 'rat', 'fly', 'worm'],
                       help='Species to generate mapping for')
    parser.add_argument('--force', action='store_true', help='Force regeneration even if output file exists')
    
    args = parser.parse_args()
    
    # Check if output already exists and we're not forcing regeneration
    if os.path.exists(args.output) and not args.force:
        logger.info(f"Output file {args.output} already exists. Use --force to regenerate.")
        return
    
    start_time = time.time()
    success = create_gene_mapping(args.output, args.species)
    
    if success:
        logger.info(f"Successfully created mapping in {time.time() - start_time:.2f} seconds")
    else:
        logger.error("Failed to create mapping")
        sys.exit(1)

if __name__ == '__main__':
    main()