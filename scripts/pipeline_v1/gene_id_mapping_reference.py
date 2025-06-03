#!/usr/bin/env python3
"""
Gene ID Reference Mapping Generator

This script creates a comprehensive gene ID reference mapping database
that serves as the single source of truth for all gene identifiers.
It combines GENCODE v24 annotations with ENCODE/ENTEx numeric ID mappings
to create a unified mapping between different gene ID formats.

Usage:
    python gene_id_mapping_reference.py \
        --encode-dir /path/to/encode/data \
        --entex-dir /path/to/entex/data \
        --gencode-gtf /path/to/gencode.v24.gtf \
        --entrez-mapping /path/to/entrez_to_ensembl_mapping.csv \
        --output /path/to/gene_id_reference_mapping.csv
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import gzip
import re
from pathlib import Path
import glob
import urllib.request
import shutil
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gene_id_mapping_generator')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate comprehensive gene ID reference mapping')
    parser.add_argument('--encode-dir', type=str, required=True,
                        help='Directory containing ENCODE data')
    parser.add_argument('--entex-dir', type=str, required=True,
                        help='Directory containing ENTEx data')
    parser.add_argument('--gencode-gtf', type=str,
                        help='Path to GENCODE v24 GTF file (will download if not provided)')
    parser.add_argument('--entrez-mapping', type=str, required=True,
                        help='Path to Entrez to Ensembl mapping CSV')
    parser.add_argument('--output', type=str, required=True,
                        help='Output file for gene ID reference mapping')
    parser.add_argument('--temp-dir', type=str, default='/tmp',
                        help='Temporary directory for downloaded files')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration of mapping even if output exists')
    
    return parser.parse_args()

def download_gencode_gtf(output_path, temp_dir='/tmp'):
    """Download GENCODE v24 GTF file if not already present."""
    logger.info("Downloading GENCODE v24 GTF file...")
    
    # GENCODE v24 URL (hg38)
    gencode_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz"
    
    # Create temp directory if it doesn't exist
    os.makedirs(temp_dir, exist_ok=True)
    
    # Download gzipped file
    temp_file = os.path.join(temp_dir, "gencode.v24.annotation.gtf.gz")
    try:
        urllib.request.urlretrieve(gencode_url, temp_file)
        
        # Decompress and save to output path
        with gzip.open(temp_file, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        logger.info(f"GENCODE v24 GTF file downloaded and saved to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Error downloading GENCODE v24 GTF file: {e}")
        return False

def parse_gencode_gtf(gtf_file):
    """Parse GENCODE GTF file to extract gene information."""
    logger.info(f"Parsing GENCODE GTF file: {gtf_file}")
    
    genes = {}
    gene_count = 0
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                
                # Extract gene information
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                
                # Parse attributes
                attr_dict = {}
                attr_string = fields[8]
                for attr in attr_string.split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    
                    key_value = attr.split(' ', 1)
                    if len(key_value) != 2:
                        continue
                    
                    key, value = key_value
                    value = value.strip('"')
                    attr_dict[key] = value
                
                # Required attributes for a gene
                if 'gene_id' not in attr_dict:
                    continue
                
                original_gene_id = attr_dict['gene_id']  # Keep full ID with version
                gene_id = original_gene_id.split('.')[0]  # Base ID for mapping                
                gene_name = attr_dict.get('gene_name', '')
                gene_type = attr_dict.get('gene_type', '')

                genes[gene_id] = {
                    'gene_id': gene_id,
                    'original_gene_id': original_gene_id,  # Add this field
                    'gene_name': gene_name,                
                    'gene_type': gene_type,
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                }
                
                gene_count += 1
                if gene_count % 10000 == 0:
                    logger.info(f"Processed {gene_count} genes...")
        
        logger.info(f"Parsed {gene_count} genes from GENCODE GTF file")
        return genes
    except Exception as e:
        logger.error(f"Error parsing GENCODE GTF file: {e}")
        return {}

def load_encode_numeric_ids(encode_dir):
    """Load numeric IDs from ENCODE dataset."""
    logger.info(f"Loading numeric IDs from ENCODE dataset: {encode_dir}")
    
    numeric_ids = set()
    
    # Look for TSV files in ENCODE directory
    tsv_files = glob.glob(os.path.join(encode_dir, '**/*.tsv'), recursive=True)
    
    if not tsv_files:
        logger.warning(f"No TSV files found in {encode_dir}")
        return numeric_ids
    
    # Process a sample of files to identify gene IDs
    for file_path in tsv_files[:10]:  # Process first 10 files
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#')
            if 'gene_id' in df.columns:
                # Check if these are numeric IDs
                sample_ids = df['gene_id'].iloc[:5].astype(str).tolist()
                if any(id.isdigit() for id in sample_ids):
                    numeric_ids.update(df['gene_id'].astype(str))
                    logger.info(f"Found numeric IDs in {file_path}")
        except Exception as e:
            logger.warning(f"Error reading {file_path}: {e}")
    
    logger.info(f"Found {len(numeric_ids)} numeric IDs in ENCODE dataset")
    return numeric_ids

def load_entex_gene_ids(entex_dir):
    """Load gene IDs from ENTEx dataset."""
    logger.info(f"Loading gene IDs from ENTEx dataset: {entex_dir}")
    
    gene_ids = defaultdict(set)  # Dictionary to store different ID types
    
    # Look for TSV files in ENTEx directory
    tsv_files = glob.glob(os.path.join(entex_dir, '**/*.tsv'), recursive=True)
    
    if not tsv_files:
        logger.warning(f"No TSV files found in {entex_dir}")
        return gene_ids
    
    # Process a sample of files to identify gene IDs
    for file_path in tsv_files[:10]:  # Process first 10 files
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#')
            if 'gene_id' in df.columns:
                # Categorize IDs by type
                for id_str in df['gene_id'].astype(str):
                    if id_str.startswith('ENSG'):
                        gene_ids['ensembl'].add(id_str.split('.')[0])  # Remove version
                    elif id_str.isdigit():
                        gene_ids['numeric'].add(id_str)
                    elif id_str.startswith('gSpikein'):
                        gene_ids['spike_in'].add(id_str)
                    else:
                        gene_ids['other'].add(id_str)
        except Exception as e:
            logger.warning(f"Error reading {file_path}: {e}")
    
    # Log counts by ID type
    for id_type, ids in gene_ids.items():
        logger.info(f"Found {len(ids)} {id_type} IDs in ENTEx dataset")
    
    return gene_ids

def load_entrez_to_ensembl_mapping(mapping_file):
    """Load Entrez to Ensembl ID mapping."""
    logger.info(f"Loading Entrez to Ensembl mapping from {mapping_file}")
    
    try:
        # Load mapping file
        mapping_df = pd.read_csv(mapping_file)
        logger.info(f"Loaded mapping with {len(mapping_df)} entries")
        
        # Create mapping dictionary
        entrez_to_ensembl = {}
        for _, row in mapping_df.iterrows():
            entrez_id = str(row['entrez_id'])
            ensembl_id = row['ensembl_id']
            
            entrez_to_ensembl[entrez_id] = ensembl_id
        
        logger.info(f"Created Entrez to Ensembl mapping with {len(entrez_to_ensembl)} entries")
        return entrez_to_ensembl
    
    except Exception as e:
        logger.error(f"Error loading Entrez to Ensembl mapping: {e}")
        return {}

def load_standardized_datasets(data_dir):
    """Load standardized datasets to extract gene ID information."""
    logger.info(f"Loading standardized datasets from {data_dir}")
    
    datasets = {}
    
    # Find standardized h5ad files
    h5ad_files = glob.glob(os.path.join(data_dir, '*_standardized_v*.h5ad'))
    
    for file_path in h5ad_files:
        dataset_name = os.path.basename(file_path).split('_')[0]
        logger.info(f"Loading {dataset_name} dataset from {file_path}")
        
        try:
            adata = sc.read_h5ad(file_path)
            datasets[dataset_name] = adata
            
            # Log dataset information
            logger.info(f"Loaded {dataset_name} with {adata.n_obs} samples and {adata.n_vars} genes")
            
            # Check if var contains gene_id column
            if 'gene_id' in adata.var.columns:
                logger.info(f"Found gene_id column in {dataset_name} var DataFrame")
            
            # Check var_names format
            var_names_sample = list(adata.var_names[:5])
            logger.info(f"Sample var_names in {dataset_name}: {var_names_sample}")
            
        except Exception as e:
            logger.error(f"Error loading {file_path}: {e}")
    
    return datasets

def create_gene_id_mapping(gencode_genes, encode_numeric_ids, entex_gene_ids, entrez_to_ensembl, standardized_datasets):
    """Create comprehensive gene ID mapping."""
    logger.info("Creating comprehensive gene ID mapping...")
    
    # Initialize mapping DataFrame with GENCODE genes
    mapping_data = []
    for gene_id, gene_info in gencode_genes.items():
        mapping_data.append({
            'gene_id': gene_id,
            'gene_name': gene_info['gene_name'],
            'gene_type': gene_info['gene_type'],
            'chromosome': gene_info['chromosome'],
            'start': gene_info['start'],
            'end': gene_info['end'],
            'strand': gene_info['strand'],
            'numeric_id': None,
            'source': 'GENCODE',
            'mapping_confidence': 'high'
        })
    
    # Add Entrez to Ensembl mappings
    entrez_added = set()
    
    for entrez_id, ensembl_id in entrez_to_ensembl.items():
        # Skip if Ensembl ID is not in GENCODE (which would be unusual)
        if ensembl_id not in gencode_genes:
            continue
        
        # Find the GENCODE entry for this Ensembl ID
        for item in mapping_data:
            if item['gene_id'] == ensembl_id:
                item['numeric_id'] = entrez_id
                entrez_added.add(entrez_id)
                break
    
    logger.info(f"Added {len(entrez_added)} Entrez to Ensembl mappings")
    
    # Add remaining Entrez IDs (ones not in GENCODE)
    missing_entrez = encode_numeric_ids - entrez_added
    logger.info(f"Found {len(missing_entrez)} Entrez IDs not in GENCODE mapping")
    
    for entrez_id in missing_entrez:
        # Add as placeholder
        mapping_data.append({
            'gene_id': f"PLACEHOLDER_{entrez_id}",
            'gene_name': f"Unknown_{entrez_id}",
            'gene_type': 'unknown',
            'chromosome': '',
            'start': 0,
            'end': 0,
            'strand': '',
            'numeric_id': entrez_id,
            'source': 'ENCODE',
            'mapping_confidence': 'low'
        })
    
    # Handle ENTEx gene IDs
    # 1. Add Ensembl IDs that aren't in GENCODE
    for ensembl_id in entex_gene_ids['ensembl']:
        if ensembl_id not in gencode_genes:
            mapping_data.append({
                'gene_id': ensembl_id,
                'gene_name': '',
                'gene_type': 'unknown',
                'chromosome': '',
                'start': 0,
                'end': 0,
                'strand': '',
                'numeric_id': None,
                'source': 'ENTEx',
                'mapping_confidence': 'medium'
            })
    
    # 2. Add spike-in controls with special handling
    for spike_in_id in entex_gene_ids['spike_in']:
        mapping_data.append({
            'gene_id': spike_in_id,
            'gene_name': spike_in_id,
            'gene_type': 'spike_in_control',
            'chromosome': 'spike_in',
            'start': 0,
            'end': 0,
            'strand': '',
            'numeric_id': None,
            'source': 'ENTEx',
            'mapping_confidence': 'high'
        })
    
    # Create DataFrame from mapping data
    mapping_df = pd.DataFrame(mapping_data)
    logger.info(f"Created mapping with {len(mapping_df)} entries")
    
    # Log some statistics
    logger.info(f"Mapped numeric IDs: {sum(mapping_df['numeric_id'].notna())}")
    logger.info(f"Mapping confidence counts: {mapping_df['mapping_confidence'].value_counts().to_dict()}")
    
    return mapping_df

def main():
    """Main function to create gene ID reference mapping."""
    args = parse_arguments()
    
    # Check if output already exists and we're not forcing regeneration
    if os.path.exists(args.output) and not args.force:
        logger.info(f"Output file {args.output} already exists. Use --force to regenerate.")
        return
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Step 1: Get GENCODE GTF file
    gencode_gtf_file = args.gencode_gtf
    if not gencode_gtf_file or not os.path.exists(gencode_gtf_file):
        temp_gtf_file = os.path.join(args.temp_dir, "gencode.v24.annotation.gtf")
        if not download_gencode_gtf(temp_gtf_file, args.temp_dir):
            logger.error("Failed to download GENCODE GTF file. Please provide it manually.")
            return
        gencode_gtf_file = temp_gtf_file
    
    # Step 2: Parse GENCODE GTF to extract gene information
    gencode_genes = parse_gencode_gtf(gencode_gtf_file)
    if not gencode_genes:
        logger.error("Failed to parse GENCODE GTF file.")
        return
    
    # Step 3: Load NCBI Entrez to Ensembl mapping
    entrez_to_ensembl = load_entrez_to_ensembl_mapping(args.entrez_mapping)
    if not entrez_to_ensembl:
        logger.error("Failed to load Entrez to Ensembl mapping.")
        return
    
    # Step 4: Load numeric IDs from ENCODE dataset
    encode_numeric_ids = load_encode_numeric_ids(args.encode_dir)
    
    # Step 5: Load gene IDs from ENTEx dataset
    entex_gene_ids = load_entex_gene_ids(args.entex_dir)
    
    # Step 6: Load standardized datasets for additional information
    standardized_dir = os.path.dirname(args.output)
    standardized_datasets = load_standardized_datasets(standardized_dir)
    
    # Step 7: Create comprehensive gene ID mapping
    mapping_df = create_gene_id_mapping(
        gencode_genes, encode_numeric_ids, entex_gene_ids, entrez_to_ensembl, standardized_datasets
    )
    
    # Step 8: Save mapping to output file
    mapping_df.to_csv(args.output, index=False)
    logger.info(f"Gene ID reference mapping saved to {args.output}")
    
    # Step 9: Create a JSON version for easier loading
    json_output = args.output.replace('.csv', '.json')
    mapping_df.to_json(json_output, orient='records')
    logger.info(f"Gene ID reference mapping also saved as JSON to {json_output}")

if __name__ == "__main__":
    main()