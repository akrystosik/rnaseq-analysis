#!/usr/bin/env python3
"""
ENCODE Processing with Local Gene Mapping

This script converts ENCODE cell line TPM data to AnnData format using 
locally created gene ID mappings instead of relying on web APIs.

Usage:
  python encode-with-local-mapping.py --input /path/to/encode/data --output /path/to/output.h5ad
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import anndata as ad
import requests
import glob
import re
import logging
import time
import json
import gzip
import io
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('encode_to_anndata')

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis")
DEFAULT_INPUT = BASE_DIR / "raw_data/gene_quantification"
DEFAULT_OUTPUT = BASE_DIR / "processed_data/anndata/encode_cell_lines_gencode.h5ad"
GENCODE_MAPPING_FILE = BASE_DIR / "gtex_analysis/metadata/gencode_v24_complete_mapping.csv"
ENTREZ_MAPPING_FILE = BASE_DIR / "metadata/entrez_to_ensembl_mapping.csv"

# Cell line metadata - Add more as needed
CELL_LINE_INFO = {
    'A549': {
        'tissue_of_origin': 'lung',
        'disease': 'adenocarcinoma',
        'cell_type': 'epithelial',
        'sex': 'male',
        'organism': 'human',
        'age': "58", 
        'ethnicity': 'European',
        'geo-id':'SAMN05733878'
    },    
    'K562': {
        'tissue_of_origin': 'bone marrow',
        'disease': 'chronic myelogenous leukemia',
        'cell_type': 'lymphoblast',
        'sex': 'female',
        'age': "53", 
        'organism': 'human',
        'ethnicity': 'unknown',
        'geo-id':'SAMN04284550'
    },
    'HepG2': {
        'tissue_of_origin': 'liver',
        'disease': 'hepatocellular carcinoma',
        'cell_type': 'epithelial',
        'sex': 'male',
        'organism': 'human',
        'age': "15", 
        'ethnicity' :'European',
        'geo-id':'SAMN04284581'
    },
    'GM23248': {
        'tissue_of_origin': '"skin of body", "limb", "connective tissue"',
        'disease': 'normal',
        'cell_type': 'Fibroblast',
        'sex': 'male',
        'organism': 'human',
        'age': "53", 
        'ethnicity' :'European',
        'geo-id':'SAMN04284514'
    },
    'Caki2': {
        'tissue_of_origin': 'kidney',
        'disease': 'clear cell renal cell carcinoma',
        'cell_type': 'epithelial',
        'sex': 'male',
        'organism': 'human',
        'age': "69", 
        'ethnicity' :'European',
        'geo-id':'SAMN04284635'
    },
    'NCI-H460': {
        'tissue_of_origin': 'lung',
        'disease': 'large cell lung carcinoma',
        'cell_type': 'epithelial',
        'sex': 'male',
        'organism': 'human',
        'age': "",  # Changed from None to empty string
        'ethnicity': 'unknown',
        'geo-id':'SAMN07791353'
        
    },
    'Panc1': {
        'tissue_of_origin': 'pancreas',
        'disease': 'pancreatic carcinoma',
        'cell_type': 'ductal cell',
        'sex': 'male',
        'organism': 'human',        
        'age': "56", 
        'ethnicity' :'European',
        'geo-id':'SAMN05733879'
    }
}

def sanitize_value(value):
    """Convert any value to a string-compatible type for h5ad storage."""
    if value is None:
        return ""  # Convert None to empty string
    elif isinstance(value, (str, int, float, bool)):
        return str(value)  # Convert ALL values to strings for consistency
    elif isinstance(value, (list, tuple)):
        # Convert lists to comma-separated strings
        return ', '.join(str(sanitize_value(v)) for v in value if v is not None)
    elif isinstance(value, dict):
        # Convert dicts to JSON strings
        try:
            return json.dumps(value)
        except:
            return str(value)
    else:
        # For other types, convert to string
        return str(value)

def sanitize_metadata(metadata):
    """Sanitize all values in a metadata dictionary to ensure h5ad compatibility."""
    sanitized = {}
    for key, value in metadata.items():
        sanitized[key] = sanitize_value(value)
    return sanitized

def prepare_for_anndata(obs_df, var_df, data_df):
    """Prepare dataframes for AnnData creation to avoid h5ad compatibility issues."""
    
    # Convert any None values to empty strings in observation metadata
    for col in obs_df.columns:
        obs_df[col] = obs_df[col].apply(lambda x: "" if x is None else str(x))
    
    # Convert any None values to empty strings in variable metadata
    for col in var_df.columns:
        var_df[col] = var_df[col].apply(lambda x: "" if x is None else str(x))
    
    # Ensure data matrix has no NaN values
    data_df = data_df.fillna(0)
    
    return obs_df, var_df, data_df

def extract_metadata(file_path, encode_metadata=None):
    """Extract metadata from file path and ENCODE metadata."""
    # Initialize metadata dictionary
    metadata = {}
    
    # Get file basename and extract file ID
    file_name = os.path.basename(file_path)
    file_id = file_name.split('.')[0]
    metadata['file_id'] = file_id
    
    # Extract cell line from path
    path_parts = file_path.split(os.sep)
    cell_line_dir = [part for part in path_parts if part in CELL_LINE_INFO or any(cl in part for cl in CELL_LINE_INFO.keys())]
    
    if cell_line_dir:
        # Handle directories like 'A549_polyA_plus'
        cell_line_full = cell_line_dir[0]
        
        # Extract cell line name
        cell_line = None
        for cl in CELL_LINE_INFO.keys():
            if cell_line_full.startswith(cl):
                cell_line = cl
                break
        
        if cell_line:
            metadata['cell_line'] = cell_line
            
            # Add basic cell line information
            if cell_line in CELL_LINE_INFO:
                for key, value in CELL_LINE_INFO[cell_line].items():
                    metadata[key] = value
            
            # Extract RNA extraction method if in directory name
            if '_polyA_plus' in cell_line_full:
                metadata['extraction_method'] = 'polyA_plus'
            elif '_polyA_minus' in cell_line_full:
                metadata['extraction_method'] = 'polyA_minus'
            elif '_total' in cell_line_full:
                metadata['extraction_method'] = 'total'
    
    # Sanitize all metadata values for h5ad compatibility
    return sanitize_metadata(metadata)

def read_tpm_file(file_path):
    """Read an ENCODE TPM file."""
    try:
        # Read the TSV file
        df = pd.read_csv(file_path, sep='\t')
        
        # Check expected columns
        if 'gene_id' in df.columns and 'TPM' in df.columns:
            # Set gene_id as index and handle numeric IDs
            df['gene_id'] = df['gene_id'].astype(str)
            df = df.set_index('gene_id')
            return df[['TPM']]
        else:
            # Try to detect columns
            gene_id_col = [col for col in df.columns if 'gene' in col.lower() and 'id' in col.lower()]
            tpm_col = [col for col in df.columns if 'tpm' in col.lower()]
            
            if gene_id_col and tpm_col:
                # Convert gene IDs to strings to ensure consistent types
                df[gene_id_col[0]] = df[gene_id_col[0]].astype(str)
                df = df.set_index(gene_id_col[0])
                return df[[tpm_col[0]]]
            else:
                logger.warning(f"Could not identify gene_id and TPM columns in {file_path}")
                return None
    
    except Exception as e:
        logger.warning(f"Error reading TPM file {file_path}: {e}")
        return None

def create_entrez_mapping():
    """
    Create a comprehensive mapping between Entrez Gene IDs and Ensembl IDs
    using multiple NCBI sources for better coverage.
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with mapping between Entrez Gene IDs and Ensembl IDs
    """
    try:
        os.makedirs(os.path.dirname(ENTREZ_MAPPING_FILE), exist_ok=True)
        
        logger.info("Creating comprehensive Entrez to Ensembl mapping file...")
        
        # Primary approach: gene2ensembl (the canonical mapping file)
        mapping_data = []
        
        # Download NCBI gene2ensembl file
        logger.info("Downloading gene2ensembl from NCBI...")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
        
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        # Parse the file
        logger.info("Parsing gene2ensembl file...")
        with gzip.open(io.BytesIO(response.content)) as f:
            # Skip header
            next(f)
            
            # Process each line
            count = 0
            for line in f:
                line = line.decode('utf-8').strip()
                fields = line.split('\t')
                
                # Check if this is a human gene (tax_id = 9606)
                if fields[0] == '9606':
                    entrez_id = fields[1]
                    ensembl_gene = fields[2]
                    
                    # Store the mapping
                    mapping_data.append({
                        'entrez_id': entrez_id,
                        'ensembl_id': ensembl_gene,
                        'source': 'gene2ensembl'
                    })
                    
                    count += 1
                    if count % 10000 == 0:
                        logger.info(f"Processed {count} human gene mappings from gene2ensembl")
        
        logger.info(f"Found {len(mapping_data)} mappings from gene2ensembl")
        
        # Secondary approach: gene_info for additional mappings
        # This file contains dbXrefs which include Ensembl IDs
        logger.info("Downloading gene_info from NCBI for additional mappings...")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
        
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        logger.info("Parsing gene_info file...")
        count = 0
        with gzip.open(io.BytesIO(response.content)) as f:
            # Skip header
            next(f)
            
            for line in f:
                line = line.decode('utf-8').strip()
                fields = line.split('\t')
                
                entrez_id = fields[1]
                dbxrefs = fields[5]
                
                # Extract Ensembl IDs from dbXrefs
                for ref in dbxrefs.split('|'):
                    if ref.startswith('Ensembl:'):
                        ensembl_id = ref.replace('Ensembl:', '')
                        mapping_data.append({
                            'entrez_id': entrez_id,
                            'ensembl_id': ensembl_id,
                            'source': 'gene_info'
                        })
                        count += 1
                
                if count % 10000 == 0 and count > 0:
                    logger.info(f"Processed {count} additional mappings from gene_info")
        
        logger.info(f"Found {count} additional mappings from gene_info")
        
        # Tertiary approach: history file for deprecated/replaced IDs
        logger.info("Downloading gene_history from NCBI for historical mappings...")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz"
        
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        # First, build a map of old to new Entrez IDs
        logger.info("Building historical Entrez ID mapping...")
        entrez_history = {}
        with gzip.open(io.BytesIO(response.content)) as f:
            # Skip header
            next(f)
            
            for line in f:
                line = line.decode('utf-8').strip()
                fields = line.split('\t')
                
                # Check if this is a human gene (tax_id = 9606)
                if fields[0] == '9606':
                    # Fields: tax_id, GeneID, Discontinued_GeneID, Discontinued_Symbol, Discontinue_Date
                    if len(fields) >= 3 and fields[2] != '-':
                        old_id = fields[2]
                        new_id = fields[1]
                        entrez_history[old_id] = new_id
        
        logger.info(f"Found {len(entrez_history)} historical Entrez ID mappings")
        
        # Create a set of all entrez IDs we've seen
        # (to check for missing IDs that might need history lookup)
        mapped_entrez_ids = set(item['entrez_id'] for item in mapping_data)
        
        # Create DataFrame and save
        mapping_df = pd.DataFrame(mapping_data)
        
        # Remove duplicates but keep track of sources
        mapping_df = mapping_df.drop_duplicates(subset=['entrez_id', 'ensembl_id'])
        
        # Save the comprehensive mapping
        mapping_df.to_csv(ENTREZ_MAPPING_FILE, index=False)
        
        logger.info(f"Created mapping with {len(mapping_df)} entries between Entrez and Ensembl IDs")
        logger.info(f"Saved mapping to {ENTREZ_MAPPING_FILE}")
        
        return mapping_df, entrez_history
    
    except Exception as e:
        logger.error(f"Error creating Entrez to Ensembl mapping: {e}")
        return pd.DataFrame(columns=['entrez_id', 'ensembl_id', 'source']), {}

def load_entrez_mapping():
    """
    Load or create mapping between Entrez Gene IDs and Ensembl IDs.
    
    Returns:
    --------
    dict
        Dictionary mapping Entrez Gene IDs to Ensembl IDs
    dict
        Dictionary mapping old Entrez IDs to new Entrez IDs (for history lookup)
    """
    # Check if mapping file exists
    if os.path.exists(ENTREZ_MAPPING_FILE):
        try:
            logger.info(f"Loading Entrez to Ensembl mapping from {ENTREZ_MAPPING_FILE}")
            mapping_df = pd.read_csv(ENTREZ_MAPPING_FILE)
            
            # Create mapping dictionary
            mapping_dict = {}
            for _, row in mapping_df.iterrows():
                entrez_id = str(row['entrez_id'])
                ensembl_id = str(row['ensembl_id'])
                
                if entrez_id not in mapping_dict:
                    mapping_dict[entrez_id] = []
                
                if ensembl_id not in mapping_dict[entrez_id]:
                    mapping_dict[entrez_id].append(ensembl_id)
            
            logger.info(f"Loaded mapping for {len(mapping_dict)} Entrez IDs")
            
            # Create empty history dict since we're loading from file
            entrez_history = {}
            
            return mapping_dict, entrez_history
        
        except Exception as e:
            logger.error(f"Error loading Entrez mapping: {e}")
            logger.info("Creating new mapping file...")
    
    # Create new mapping
    return create_entrez_mapping()

def load_gencode_mapping():
    """
    Load GENCODE mapping from CSV file.
    
    Returns:
    --------
    dict
        Dictionary mapping Ensembl IDs to gene information
    """
    try:
        if not os.path.exists(GENCODE_MAPPING_FILE):
            logger.warning(f"GENCODE mapping file not found: {GENCODE_MAPPING_FILE}")
            return {}
        
        logger.info(f"Loading GENCODE mapping from {GENCODE_MAPPING_FILE}")
        gencode_df = pd.read_csv(GENCODE_MAPPING_FILE)
        
        # Create dictionary for faster lookups
        gencode_dict = {}
        
        for _, row in gencode_df.iterrows():
            if 'gene_id' in row:
                # Get the gene ID (with or without version)
                gene_id = row['gene_id']
                
                # Also extract the base ID without version
                base_id = gene_id.split('.')[0] if '.' in gene_id else gene_id
                
                info = {
                    'gene_name': row['gene_name'] if 'gene_name' in row else '',
                    'gene_type': row['gene_type'] if 'gene_type' in row else '',
                    'chromosome': row['chromosome'] if 'chromosome' in row else ''
                }
                
                # Map both the full ID and base ID
                gencode_dict[gene_id] = info
                gencode_dict[base_id] = info
        
        logger.info(f"Loaded GENCODE mapping for {len(gencode_dict)} gene IDs")
        return gencode_dict
    
    except Exception as e:
        logger.warning(f"Error loading GENCODE mapping: {e}")
        return {}

def process_encode_files(input_dir, output_file, cell_lines=None):
    """
    Process ENCODE TPM files that already use Ensembl IDs and create an AnnData object.
    
    Parameters:
    -----------
    input_dir : str
        Path to the directory containing ENCODE data
    output_file : str
        Path to save the AnnData file
    cell_lines : list, optional
        List of cell lines to process (if None, process all)
    
    Returns:
    --------
    anndata.AnnData
        AnnData object with TPM data and metadata
    """
    logger.info(f"Processing ENCODE files from {input_dir}")
    
    # Find all TSV files
    if cell_lines:
        # Only specified cell lines
        tsv_files = []
        for cl in cell_lines:
            tsv_files.extend(glob.glob(os.path.join(input_dir, f"{cl}*", "*.tsv")))
            tsv_files.extend(glob.glob(os.path.join(input_dir, cl, "*.tsv")))
    else:
        # All cell lines
        tsv_files = glob.glob(os.path.join(input_dir, "**", "*.tsv"), recursive=True)
    
    # Remove duplicates
    tsv_files = list(set(tsv_files))
    
    if not tsv_files:
        logger.error(f"No TSV files found in {input_dir}")
        return None
    
    logger.info(f"Found {len(tsv_files)} TSV files to process")
    
    # Load GENCODE mapping for gene annotation
    gencode_mapping = load_gencode_mapping()
    
    # Process each file
    sample_dfs = {}  # Use dictionary to track sample data frames
    sample_metadata = {}  # Use dictionary to track sample metadata
    
    # Track all gene IDs we've seen
    all_gene_ids = set()
    
    for file_path in tsv_files:
        # Get file ID
        file_id = os.path.basename(file_path).split('.')[0]
        
        logger.info(f"Processing {file_path}")
        
        # Read TPM data
        tpm_df = read_tpm_file(file_path)
        if tpm_df is None:
            continue
        
        # Rename TPM column to file ID
        col_name = tpm_df.columns[0]
        tpm_df = tpm_df.rename(columns={col_name: file_id})
        
        # Add to collections
        sample_dfs[file_id] = tpm_df
        
        # Extract metadata
        metadata = extract_metadata(file_path)
        sample_metadata[file_id] = metadata
        
        # Track all gene IDs
        all_gene_ids.update(tpm_df.index)
    
    if not sample_dfs:
        logger.error("No valid TPM data found")
        return None
    
    # Process gene IDs
    logger.info(f"Processing {len(all_gene_ids)} gene IDs")
    
    # Standardize Ensembl IDs - handle version numbers
    gene_id_mapping = {}
    base_ids = {}
    version_stats = {'with_version': 0, 'without_version': 0}
    
    for gene_id in all_gene_ids:
        # Check if this is an Ensembl ID with version
        if '.' in gene_id and gene_id.startswith('ENSG'):
            base_id = gene_id.split('.')[0]
            version_stats['with_version'] += 1
        else:
            base_id = gene_id
            version_stats['without_version'] += 1
        
        # Store both forms of the ID
        gene_id_mapping[gene_id] = base_id
        base_ids[base_id] = base_id
    
    logger.info(f"Found {version_stats['with_version']} Ensembl IDs with version numbers")
    logger.info(f"Found {version_stats['without_version']} Ensembl IDs without version numbers")
    
    # Create a single unified expression matrix
    logger.info("Creating unified expression matrix")
    
    # Use base Ensembl IDs (without version) for the unified matrix
    unique_base_ids = sorted(base_ids.keys())
    
    # OPTIMIZATION: Create DataFrame efficiently
    # Create a dictionary to collect data for each sample
    sample_data = {}
    
    # Add expression data for each sample using base IDs
    for file_id, tpm_df in sample_dfs.items():
        # Map expression values to base Ensembl IDs
        ensembl_expr = {}
        
        for gene_id, tpm in zip(tpm_df.index, tpm_df[file_id]):
            base_id = gene_id_mapping.get(gene_id, gene_id)
            
            # In case of multiple versions mapping to same base ID, use the maximum TPM value
            if base_id in ensembl_expr:
                ensembl_expr[base_id] = max(ensembl_expr[base_id], tpm)
            else:
                ensembl_expr[base_id] = tpm
        
        # Store in the dictionary
        sample_data[file_id] = ensembl_expr
    
    # Create the DataFrame efficiently
    data_array = []
    for base_id in unique_base_ids:
        row = [sample_data[file_id].get(base_id, 0.0) for file_id in sample_dfs.keys()]
        data_array.append(row)
    
    ensembl_df = pd.DataFrame(data_array, index=unique_base_ids, columns=sample_dfs.keys())
    
    # Create observation metadata
    obs_df = pd.DataFrame([sample_metadata[file_id] for file_id in ensembl_df.columns])
    obs_df.index = ensembl_df.columns
    
    # Create gene metadata (var)
    var_df = pd.DataFrame(index=ensembl_df.index)
    var_df['ensembl_id'] = var_df.index
    
    # Add GENCODE information if available
    if gencode_mapping:
        logger.info("Adding GENCODE gene information")
        
        gene_names = []
        gene_types = []
        chromosomes = []
        
        for ensembl_id in var_df.index:
            # Check if this ID is in the GENCODE mapping
            if ensembl_id in gencode_mapping:
                info = gencode_mapping[ensembl_id]
                gene_names.append(info['gene_name'])
                gene_types.append(info['gene_type'])
                chromosomes.append(info['chromosome'])
            else:
                gene_names.append("")
                gene_types.append("")
                chromosomes.append("")
        
        var_df['gene_name'] = gene_names
        var_df['gene_type'] = gene_types
        var_df['chromosome'] = chromosomes
    
    # Create AnnData object
    logger.info("Creating AnnData object")
    
    # Final processing of the DataFrames to avoid potential issues
    obs_df, var_df, ensembl_df = prepare_for_anndata(obs_df, var_df, ensembl_df)
    
    # Make sure all column names are unique strings
    obs_df.columns = [str(col) for col in obs_df.columns]
    var_df.columns = [str(col) for col in var_df.columns]
    
    # Transpose the data matrix (samples as rows, genes as columns)
    X = ensembl_df.T.values
    
    # Create AnnData object
    adata = ad.AnnData(
        X=X,
        obs=obs_df,
        var=var_df
    )
    
    # Add additional information to uns
    adata.uns['dataset'] = 'ENCODE_cell_lines'
    adata.uns['expression_type'] = 'TPM'
    adata.uns['conversion_date'] = pd.Timestamp.now().strftime('%Y-%m-%d')
    adata.uns['gene_id_type'] = 'Ensembl'
    
    # Save the AnnData object
    logger.info(f"Saving AnnData object to {output_file}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    adata.write_h5ad(output_file)
    
    logger.info(f"Created AnnData object with {adata.n_obs} samples and {adata.n_vars} genes")
    
    # Print summary of cell lines
    if 'cell_line' in adata.obs.columns:
        cell_line_counts = adata.obs['cell_line'].value_counts()
        logger.info("Cell line distribution:")
        for cell_line, count in cell_line_counts.items():
            logger.info(f"  {cell_line}: {count} samples")
    
    return adata
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='ENCODE Cell Line Data to AnnData Conversion with GENCODE Mapping')
    parser.add_argument('--input', default=DEFAULT_INPUT, help='Path to the ENCODE data directory')
    parser.add_argument('--output', default=DEFAULT_OUTPUT, help='Path to save the AnnData file')
    parser.add_argument('--cell-lines', nargs='+', help='Cell lines to process (if not specified, process all)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    start_time = time.time()
    adata = process_encode_files(args.input, args.output, args.cell_lines)
    
    if adata is not None:
        logger.info(f"Successfully created AnnData object in {time.time() - start_time:.2f} seconds")

if __name__ == '__main__':
    main()