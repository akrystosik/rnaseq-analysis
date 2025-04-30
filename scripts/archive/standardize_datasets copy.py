#!/usr/bin/env python3
"""
Multi-Dataset Standardization Pipeline

This script processes RNA-seq and expression data from multiple sources 
(ENCODE, GTEx, MAGE, ADNI) into standardized AnnData objects.

Usage:
  python standardize_datasets.py --encode-dir /path/to/encode/data \
                                 --gtex-file /path/to/gtex.gct.gz \
                                 --mage-dir /path/to/mage/data \
                                 --adni-dir /path/to/adni/data \
                                 --output-dir /path/to/output
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import anndata as ad
import glob
import re
import logging
import time
import gzip
import json
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('data_standardizer')

# Define paths and constants
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
GENCODE_MAPPING_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gencode_v24_complete_mapping.csv")
DEFAULT_OUTPUT_DIR = BASE_DIR / "standardized_data"
ENTEX_METADATA_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/metadata/entex_manual_metadata.json")

# Common metadata fields across all datasets
CORE_METADATA_FIELDS = [
    'sample_id',       # Unique identifier for each sample
    'subject_id',      # Identifier for the subject/donor
    'sex',             # Biological sex (standardized as 'male', 'female', 'unknown')
    'age',             # Age at sample collection
    'tissue',          # Tissue or brain region
    'dataset',         # Source dataset identifier
    'data_type',       # Type of data (RNA-seq, microarray, etc.)
    'expression_unit'  # Unit of measurement (TPM, FPKM, counts, etc.)
]

# ENCODE cell line metadata
ENCODE_CELL_INFO = {
    'A549': {
        'tissue': 'lung',
        'disease': 'adenocarcinoma',
        'cell_type': 'epithelial',
        'sex': 'male',
        'organism': 'human',
        'age': "58", 
        'ethnicity': 'European',
        'subject_id': 'ENCDO451RUA',
        'geo_id': 'SAMN05733878'
    },    
    'K562': {
        'tissue': 'bone marrow',
        'disease': 'chronic myelogenous leukemia',
        'cell_type': 'lymphoblast',
        'sex': 'female',
        'age': "53", 
        'organism': 'human',
        'ethnicity': 'unknown',
        'subject_id': 'ENCDO000AAL',
        'geo_id': 'SAMN04284550'
    },
    'HepG2': {
        'tissue': 'liver',
        'disease': 'hepatocellular carcinoma',
        'cell_type': 'epithelial',
        'sex': 'male',
        'organism': 'human',
        'age': "15", 
        'ethnicity': 'European',
        'subject_id': 'ENCDO886MPB',
        'geo_id': 'SAMN04284581'
    },
    'GM23248': {
        'tissue': 'skin',  # Simplified from multiple tissues
        'disease': 'normal',
        'cell_type': 'Fibroblast',
        'sex': 'male',
        'organism': 'human',
        'age': "53", 
        'ethnicity': 'European',
        'subject_id': 'ENCDO467QPX',
        'geo_id': 'SAMN04284514'
    },
    'Caki2': {
        'tissue': 'kidney',
        'disease': 'clear cell renal cell carcinoma',
        'cell_type': 'epithelial',
        'sex': 'male',
        'organism': 'human',
        'age': "69", 
        'ethnicity': 'European',
        'subject_id': 'ENCDO869AAI',
        'geo_id': 'SAMN04284635'
    },
    'NCI-H460': {
        'tissue': 'lung',
        'disease': 'large cell lung carcinoma',
        'cell_type': 'epithelial',
        'sex': 'male',
        'organism': 'human',
        'age': "", 
        'ethnicity': 'unknown',
        'subject_id': 'ENCDO105AAA',
        'geo_id': 'SAMN07791353'
    },
    'Panc1': {
        'tissue': 'pancreas',
        'disease': 'pancreatic carcinoma',
        'cell_type': 'ductal cell',
        'sex': 'male',
        'organism': 'human',        
        'age': "56", 
        'ethnicity': 'European',
        'subject_id': 'ENCDO845WKR',
        'geo_id': 'SAMN05733879'
    }
}

# GTEx metadata files
GTEx_SAMPLE_ATTRS = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
GTEx_SUBJECT_ATTRS = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")


# ====================================================================================
# Shared Utility Functions
# ====================================================================================

def standardize_ensembl_id(gene_id):
    """Standardize Ensembl gene ID by removing version number."""
    if pd.isna(gene_id):
        return None
    
    # Convert to string if needed
    gene_id = str(gene_id)
    
    # Remove version if present
    if '.' in gene_id and gene_id.startswith('ENSG'):
        return gene_id.split('.')[0]
    return gene_id

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
                base_id = standardize_ensembl_id(gene_id)
                
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

def add_gencode_annotations(var_df, gencode_mapping):
    """
    Add GENCODE annotations to variable DataFrame - optimized version.
    """
    # Pre-allocate arrays for efficiency
    gene_count = len(var_df)
    gene_names = np.empty(gene_count, dtype=object)
    gene_types = np.empty(gene_count, dtype=object)
    chromosomes = np.empty(gene_count, dtype=object)
    mapping_sources = np.empty(gene_count, dtype=object)
    
    for i, (index, _) in enumerate(var_df.iterrows()):
        ensembl_id = index
        
        # Try exact match first
        if ensembl_id in gencode_mapping:
            info = gencode_mapping[ensembl_id]
            gene_names[i] = info['gene_name']
            gene_types[i] = info['gene_type']
            chromosomes[i] = info['chromosome']
            mapping_sources[i] = 'exact_match'
        else:
            # Try base ID match (without version)
            base_id = standardize_ensembl_id(ensembl_id)
            if base_id in gencode_mapping:
                info = gencode_mapping[base_id]
                gene_names[i] = info['gene_name']
                gene_types[i] = info['gene_type']
                chromosomes[i] = info['chromosome']
                mapping_sources[i] = 'base_id_match'
            else:
                # No mapping found
                gene_names[i] = ""
                gene_types[i] = ""
                chromosomes[i] = ""
                mapping_sources[i] = 'unmapped'
    
    # Add annotations to DataFrame efficiently
    var_df['gene_name'] = gene_names
    var_df['gene_type'] = gene_types
    var_df['chromosome'] = chromosomes
    var_df['mapping_source'] = mapping_sources
    
    return var_df

def standardize_metadata(metadata_df, dataset_name):
    """
    Standardize metadata across datasets to ensure consistent fields and values.
    
    Parameters:
    -----------
    metadata_df : pandas.DataFrame
        Metadata DataFrame
    dataset_name : str
        Source dataset name
    
    Returns:
    --------
    pandas.DataFrame
        Standardized metadata DataFrame
    """
    # Make a copy to avoid modifying the original
    df = metadata_df.copy()
    
    # Add dataset identifier
    df['dataset'] = dataset_name
    
    # Ensure all core fields exist
    for field in CORE_METADATA_FIELDS:
        if field not in df.columns:
            df[field] = ""
    
    # Standardize sex values
    if 'sex' in df.columns:
        # Map various encodings to standard values
        sex_mapping = {
            1: 'male', '1': 'male', 'male': 'male', 'm': 'male', 'Male': 'male',
            2: 'female', '2': 'female', 'female': 'female', 'f': 'female', 'Female': 'female',
            0: 'unknown', '0': 'unknown', 'unknown': 'unknown', 'u': 'unknown', np.nan: 'unknown'
        }
        df['sex'] = df['sex'].map(lambda x: sex_mapping.get(x, 'unknown'))
    
    # Convert age to string to ensure consistent format
    if 'age' in df.columns:
        df['age'] = df['age'].astype(str)
        df['age'] = df['age'].replace('nan', '')
    
    # Convert categorical fields to categorical data type
    categorical_fields = ['sex', 'tissue', 'dataset', 'cell_type', 'disease', 'ethnicity', 'data_type']
    for field in categorical_fields:
        if field in df.columns:
            df[field] = pd.Categorical(df[field].astype(str))
    
    return df

def prepare_for_anndata(obs_df, var_df, data_df):
    """
    Prepare dataframes for AnnData creation with consistent data types.
    Updated to fix categorical deprecation warnings.
    
    Parameters:
    -----------
    obs_df : pandas.DataFrame
        Observation metadata (samples)
    var_df : pandas.DataFrame
        Variable metadata (genes)
    data_df : pandas.DataFrame
        Expression data
    
    Returns:
    --------
    tuple
        Processed observation, variable, and data DataFrames
    """
    # Ensure all observation metadata columns have appropriate types
    for col in obs_df.columns:
        if pd.api.types.is_numeric_dtype(obs_df[col]):
            # Keep numeric columns as is, but fill NaNs
            obs_df[col] = obs_df[col].fillna(0)
        elif isinstance(obs_df[col].dtype, pd.CategoricalDtype):  # Updated check
            # For categorical columns, we need to convert to string first,
            # fill NaNs, then convert back to categorical
            categories = list(obs_df[col].cat.categories)
            # Check if empty string is already in categories
            if "" not in categories:
                categories = categories + [""]
            
            # Convert to string, fill NaNs with empty string, then back to categorical with updated categories
            obs_df[col] = pd.Categorical(
                obs_df[col].astype(str).fillna(""),
                categories=categories
            )
        else:
            # For object/string columns, simply fill NaNs with empty string
            obs_df[col] = obs_df[col].fillna("")
    
    # Ensure all variable metadata columns have appropriate types
    for col in var_df.columns:
        if col != 'gene_id':  # Don't convert index/gene_id
            if isinstance(var_df[col].dtype, pd.CategoricalDtype):  # Updated check
                # For categorical columns, handle the same way as above
                categories = list(var_df[col].cat.categories)
                if "" not in categories:
                    categories = categories + [""]
                
                var_df[col] = pd.Categorical(
                    var_df[col].astype(str).fillna(""),
                    categories=categories
                )
            elif pd.api.types.is_numeric_dtype(var_df[col]):
                # Keep numeric columns as is, but fill NaNs
                var_df[col] = var_df[col].fillna(0)
            else:
                # For object/string columns, fill NaNs with empty string
                var_df[col] = var_df[col].fillna("")
    
    # Ensure data matrix has no NaN values
    data_df = data_df.fillna(0)
    
    return obs_df, var_df, data_df

def create_standard_anndata(data_df, obs_df, var_df, dataset_info):
    """
    Create standardized AnnData object with consistent structure.
    
    Parameters:
    -----------
    data_df : pandas.DataFrame
        Expression data (samples x genes)
    obs_df : pandas.DataFrame
        Observation metadata (samples)
    var_df : pandas.DataFrame
        Variable metadata (genes)
    dataset_info : dict
        Dataset information for uns slot
    
    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    # Prepare data for AnnData
    obs_df, var_df, data_df = prepare_for_anndata(obs_df, var_df, data_df)
    
    # Create AnnData object (assuming data_df is samples x genes)
    adata = ad.AnnData(
        X=data_df.values,
        obs=obs_df,
        var=var_df
    )
    
    # Add dataset information to uns
    adata.uns['dataset_info'] = dataset_info
    
    # Add processing timestamp
    adata.uns['processing_date'] = pd.Timestamp.now().strftime('%Y-%m-%d')
    
    return adata

def save_anndata(adata, file_path):
    """Save AnnData object to file."""
    try:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        
        # Save the AnnData object
        adata.write_h5ad(file_path)
        
        logger.info(f"Saved AnnData to {file_path}")
        return True
    
    except Exception as e:
        logger.error(f"Error saving AnnData: {e}")
        return False

# ====================================================================================
# ENCODE-specific Processing
# ====================================================================================
def load_entex_metadata():
    """
    Load ENTEx metadata from JSON file.
    
    Returns:
    --------
    dict
        Dictionary containing ENTEx sample metadata and donor information
    """
    try:
        if not ENTEX_METADATA_FILE.exists():
            logger.warning(f"ENTEx metadata file not found: {ENTEX_METADATA_FILE}")
            return {"donor_map": {}, "entex_metadata": []}
        
        logger.info(f"Loading ENTEx metadata from {ENTEX_METADATA_FILE}")
        with open(ENTEX_METADATA_FILE, 'r') as f:
            metadata = json.load(f)
        
        # Create a lookup dictionary for quick access
        sample_lookup = {}
        for sample in metadata.get('entex_metadata', []):
            sample_id = sample.get('sample_id')
            if sample_id:
                sample_lookup[sample_id] = sample
        
        metadata['sample_lookup'] = sample_lookup
        
        logger.info(f"Loaded metadata for {len(sample_lookup)} ENTEx samples")
        return metadata
    
    except Exception as e:
        logger.error(f"Error loading ENTEx metadata: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {"donor_map": {}, "entex_metadata": [], "sample_lookup": {}}

def read_encode_tpm_file(file_path):
    """Read an ENCODE TPM file - optimized version."""
    try:
        # Use the low_memory=False option to avoid DtypeWarning
        df = pd.read_csv(file_path, sep='\t', low_memory=False)
        
        # Check expected columns
        if 'gene_id' in df.columns and 'TPM' in df.columns:
            # Set gene_id as index and handle numeric IDs
            df['gene_id'] = df['gene_id'].astype(str)
            # Extract only the columns we need before setting index (more efficient)
            df = df[['gene_id', 'TPM']].set_index('gene_id')
            return df
        else:
            # Try to detect columns
            gene_id_col = [col for col in df.columns if 'gene' in col.lower() and 'id' in col.lower()]
            tpm_col = [col for col in df.columns if 'tpm' in col.lower()]
            
            if gene_id_col and tpm_col:
                # Extract only the columns we need before conversion (more efficient)
                df = df[[gene_id_col[0], tpm_col[0]]]
                # Convert gene IDs to strings to ensure consistent types
                df[gene_id_col[0]] = df[gene_id_col[0]].astype(str)
                df = df.set_index(gene_id_col[0])
                return df
            else:
                logger.warning(f"Could not identify gene_id and TPM columns in {file_path}")
                return None
    
    except Exception as e:
        logger.warning(f"Error reading TPM file {file_path}: {e}")
        return None

def extract_encode_metadata(file_path, entex_metadata=None):
    """
    Extract metadata from ENCODE file path, now handling ENTEx samples.
    
    Parameters:
    -----------
    file_path : str
        Path to the ENCODE file
    entex_metadata : dict, optional
        ENTEx metadata with sample_lookup dictionary
    
    Returns:
    --------
    dict
        Dictionary of metadata
    """
    metadata = {}
    
    # Get file basename and extract file ID
    file_name = os.path.basename(file_path)
    file_id = file_name.split('.')[0]  # e.g., ENCFF008RKC
    
    # Store original sample ID
    metadata['original_sample_id'] = file_id
    
    # Check if this is an ENTEx file with available metadata
    if entex_metadata and 'sample_lookup' in entex_metadata and 'entex' in file_path.lower():
        if file_id in entex_metadata['sample_lookup']:
            sample_info = entex_metadata['sample_lookup'][file_id]
            donor_id = sample_info.get('donor_id')
            
            # Add basic metadata from lookup
            metadata['tissue'] = sample_info.get('tissue', '')
            metadata['subject_id'] = donor_id
            metadata['assay_type'] = sample_info.get('assay_type', 'total')
            metadata['genome_annotation'] = sample_info.get('genome_annotation', '')
            metadata['genome_assembly'] = sample_info.get('genome_assembly', '')
            metadata['data_source'] = 'ENTEx'
            
            # Add donor information if available
            if donor_id and donor_id in entex_metadata.get('donor_map', {}):
                donor_info = entex_metadata['donor_map'][donor_id]
                metadata['sex'] = donor_info.get('sex', '')
                metadata['age'] = donor_info.get('age', '')
                
            logger.debug(f"Found ENTEx metadata for sample {file_id}")
            
        else:
            # Extract minimal metadata from the directory structure
            logger.warning(f"No metadata found for ENTEx sample {file_id}, extracting from directory")
            
            # Extract tissue from directory name (e.g., "ENTEx_adrenal gland_unknown")
            path_parts = file_path.split(os.sep)
            entex_parts = [part for part in path_parts if part.startswith('ENTEx_')]
            
            if entex_parts:
                # Parse tissue from directory name
                tissue_parts = entex_parts[0].split('_')
                if len(tissue_parts) > 1:
                    metadata['tissue'] = tissue_parts[1]
                
                # Add ENTEx flag
                metadata['data_source'] = 'ENTEx'
            
            # Try to extract donor ID from path
            donor_parts = [part for part in path_parts if part.startswith('entex_tc_')]
            if donor_parts:
                metadata['subject_id'] = donor_parts[0].replace('entex_tc_', '')
    
    else:
        # Regular ENCODE cell line processing
        # Extract cell line from path
        path_parts = file_path.split(os.sep)
        cell_line_dir = [part for part in path_parts if part in ENCODE_CELL_INFO or 
                         any(cl in part for cl in ENCODE_CELL_INFO.keys())]
        
        if cell_line_dir:
            # Handle directories like 'A549_polyA_plus'
            cell_line_full = cell_line_dir[0]
            
            # Extract cell line name
            cell_line = None
            for cl in ENCODE_CELL_INFO.keys():
                if cell_line_full.startswith(cl):
                    cell_line = cl
                    break
            
            if cell_line:
                metadata['cell_line'] = cell_line
                
                # Add basic cell line information
                if cell_line in ENCODE_CELL_INFO:
                    for key, value in ENCODE_CELL_INFO[cell_line].items():
                        metadata[key] = value
                
                # Extract RNA extraction method if in directory name
                if '_polyA_plus' in cell_line_full:
                    metadata['extraction_method'] = 'polyA_plus'
                elif '_polyA_minus' in cell_line_full:
                    metadata['extraction_method'] = 'polyA_minus'
                elif '_total' in cell_line_full:
                    metadata['extraction_method'] = 'total'
    
    # Add standard fields
    metadata['data_type'] = 'RNA-seq'
    metadata['expression_unit'] = 'TPM'
    
    return metadata

def process_encode_data(input_dir, entex_dir=None, output_file=None):
    """
    Process ENCODE RNA-seq data into standardized AnnData with optimizations.
    Now includes ENTEx samples.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing ENCODE cell line TPM files
    entex_dir : str, optional
        Directory containing ENCODE ENTEx TPM files
    output_file : str
        Path to save the standardized AnnData
    
    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    logger.info(f"Processing ENCODE data from {input_dir}")
    
    # Load ENTEx metadata
    entex_metadata = load_entex_metadata() if entex_dir else None
    
    # Find all TSV files in cell lines
    tsv_files = glob.glob(os.path.join(input_dir, "**", "*.tsv"), recursive=True)
    
    # Also look for ENTEx files if directory is provided
    if entex_dir:
        logger.info(f"Checking ENTEx directory: {entex_dir}")
        entex_files = glob.glob(os.path.join(entex_dir, "**", "*.tsv"), recursive=True)
        logger.info(f"Found {len(entex_files)} TSV files in ENTEx directory")
        
        # Check metadata file
        if ENTEX_METADATA_FILE.exists():
            logger.info(f"ENTEx metadata file exists: {ENTEX_METADATA_FILE}")
        else:
            logger.error(f"ENTEx metadata file not found: {ENTEX_METADATA_FILE}")
        
        # Filter for GENCODE v24 files to match our standardization
        if entex_metadata and 'sample_lookup' in entex_metadata:
            v24_samples = []
            for file_path in entex_files:
                file_name = os.path.basename(file_path)
                file_id = file_name.split('.')[0]
                
                if file_id in entex_metadata['sample_lookup']:
                    sample_info = entex_metadata['sample_lookup'][file_id]
                    if sample_info.get('genome_annotation') == 'V24':
                        v24_samples.append(file_path)
                        
            logger.info(f"Selected {len(v24_samples)} ENTEx files with GENCODE V24 annotation")
            tsv_files.extend(v24_samples)
        else:
            # If no metadata filtering possible, add all ENTEx files
            tsv_files.extend(entex_files)
    
    if not tsv_files:
        logger.error(f"No TSV files found")
        return None
    
    logger.info(f"Found {len(tsv_files)} TPM files to process")
    
    # Process each file
    sample_dfs = {}  # Data frames for each sample
    sample_metadata = {}  # Dictionary for metadata, keyed by sample_id
    
    # Track all gene IDs
    all_gene_ids = set()
    
    for file_path in tsv_files:
        # Read TPM data
        tpm_df = read_encode_tpm_file(file_path)
        if tpm_df is None:
            continue
        
        # Extract sample ID
        file_name = os.path.basename(file_path)
        sample_id = file_name.split('.')[0]
        
        # Rename TPM column to sample ID
        col_name = tpm_df.columns[0]
        tpm_df = tpm_df.rename(columns={col_name: sample_id})
        
        # Add to collections
        sample_dfs[sample_id] = tpm_df
        
        # Extract and store metadata with ENTEx support
        metadata = extract_encode_metadata(file_path, entex_metadata)
        sample_metadata[sample_id] = metadata
        
        # Track gene IDs
        all_gene_ids.update(tpm_df.index)
    
    if not sample_dfs:
        logger.error("No valid TPM data found")
        return None
    
    
    # Standardize gene IDs and create expression matrix
    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs")
    
    # Map original IDs to standardized IDs
    gene_id_mapping = {gene_id: standardize_ensembl_id(gene_id) for gene_id in all_gene_ids}
    
    # Get unique standardized IDs
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")
    
    # OPTIMIZED: Create a numpy array to hold the data
    # Shape: (unique_std_ids x sample_ids)
    num_genes = len(unique_std_ids)
    sample_ids = list(sample_dfs.keys())
    num_samples = len(sample_ids)
    
    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    sample_id_to_idx = {id: idx for idx, id in enumerate(sample_ids)}
    
    # Pre-allocate the expression matrix with zeros (as floats)
    expr_matrix = np.zeros((num_genes, num_samples), dtype=np.float32)
    
    # Fill the matrix efficiently
    logger.info("Creating unified expression matrix")
    for sample_id, tpm_df in sample_dfs.items():
        sample_idx = sample_id_to_idx[sample_id]
        
        # Create temporary dictionary to collect values for this sample
        sample_data = {}
        
        # Process each gene in this sample
        for gene_id, tpm_val in zip(tpm_df.index, tpm_df[sample_id]):
            if gene_id in gene_id_mapping:
                std_id = gene_id_mapping[gene_id]
                
                # If multiple original IDs map to same standardized ID, use maximum value
                if std_id in sample_data:
                    sample_data[std_id] = max(sample_data[std_id], tpm_val)
                else:
                    sample_data[std_id] = tpm_val
        
        # Update the matrix all at once
        for std_id, value in sample_data.items():
            if std_id in std_id_to_idx:  # Handle case if std_id isn't in our target set
                gene_idx = std_id_to_idx[std_id]
                expr_matrix[gene_idx, sample_idx] = value
    
    # Create the DataFrame from the matrix
    unified_df = pd.DataFrame(
        expr_matrix, 
        index=unique_std_ids, 
        columns=sample_ids
    )
    
    # Create observation metadata DataFrame - ensuring it matches the sample IDs in unified_df
    # This is critical - the rows in obs_df must exactly match the columns in unified_df
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            # Make a copy of the metadata to avoid modifying the original
            metadata_copy = sample_metadata[sample_id].copy()
            
            # If there's a sample_id column, rename it to avoid conflict with index
            if 'sample_id' in metadata_copy:
                metadata_copy['original_sample_id'] = metadata_copy.pop('sample_id')
                
            # Add the sample_id to each record
            metadata_copy['_sample_id'] = sample_id
            obs_data.append(metadata_copy)
        else:
            # If we somehow have a sample in the matrix but not in metadata,
            # create minimal metadata for it
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({'_sample_id': sample_id, 'dataset': 'ENCODE'})
    
    # Create obs_df with explicit index to ensure alignment
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index('_sample_id', inplace=True)
    
    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        logger.error(f"Metadata has {len(obs_df)} samples, expression matrix has {len(unified_df.columns)} samples")
        
        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]
        logger.info(f"Fixed alignment - now using {len(obs_df)} samples")
    
    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df['gene_id'] = var_df.index
    
    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)
    
    # Add original IDs to variable metadata
    var_df['original_ids'] = var_df.index.map(lambda x: ';'.join(reverse_mapping.get(x, [])))
    
    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)
    
    # Standardize metadata
    obs_df = standardize_metadata(obs_df, 'ENCODE')
    
    # Dataset info for uns
    dataset_info = {
        'source': 'ENCODE',
        'gencode_version': 24,
        'data_type': 'RNA-seq',
        'expression_unit': 'TPM',
        'samples': len(obs_df),
        'genes': len(var_df)
    }
    
    # Final dimension check before creating AnnData
    logger.info(f"Final dimensions check: {len(obs_df)} samples x {len(var_df)} genes")
    logger.info(f"Expression matrix shape: {unified_df.T.shape}")
    
    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)
    
    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed ENCODE data with {adata.n_obs} samples and {adata.n_vars} genes")
    
    return adata

# ====================================================================================
# GTEx-specific Processing - Enhanced Version
# ====================================================================================

def load_gtex_metadata():
    """Load and standardize GTEx metadata with improved tissue mapping."""
    logger.info("Loading GTEx metadata")
    
    try:
        # Check if files exist
        if not GTEx_SAMPLE_ATTRS.exists():
            logger.error(f"GTEx sample attributes file not found: {GTEx_SAMPLE_ATTRS}")
            return pd.DataFrame()
            
        if not GTEx_SUBJECT_ATTRS.exists():
            logger.error(f"GTEx subject attributes file not found: {GTEx_SUBJECT_ATTRS}")
            return pd.DataFrame()
        
        # Load sample attributes
        logger.info(f"Loading sample attributes from {GTEx_SAMPLE_ATTRS}")
        sample_attrs = pd.read_csv(GTEx_SAMPLE_ATTRS, sep='\t')
        
        # Load subject attributes
        logger.info(f"Loading subject attributes from {GTEx_SUBJECT_ATTRS}")
        subject_attrs = pd.read_csv(GTEx_SUBJECT_ATTRS, sep='\t')
        
        # Extract SUBJID from SAMPID
        sample_attrs['SUBJID'] = sample_attrs['SAMPID'].str.extract(r'(GTEX-[A-Z0-9]+)')
        
        # Merge sample and subject attributes
        merged_attrs = pd.merge(sample_attrs, subject_attrs, on='SUBJID', how='left')
        
        # Standardize column names
        column_mapping = {
            'SAMPID': 'sample_id',
            'SUBJID': 'subject_id',
            'SMTSD': 'tissue',  # Detailed tissue site
            'SMTS': 'broad_tissue',  # Broad tissue type
            'SEX': 'sex',
            'AGE': 'age',
            'SMGEBTCH': 'batch',
            'SMRIN': 'rna_integrity_number',
            'SMTSISCH': 'ischemic_time',
            'SMNABTCH': 'array_batch'
        }
        
        # Rename columns that exist in the DataFrame
        existing_cols = [col for col in column_mapping.keys() if col in merged_attrs.columns]
        merged_attrs = merged_attrs.rename(columns={col: column_mapping[col] for col in existing_cols})
        
        # Ensure we have a tissue column
        if 'tissue' not in merged_attrs.columns and 'broad_tissue' in merged_attrs.columns:
            merged_attrs['tissue'] = merged_attrs['broad_tissue']
            logger.info("Using broad_tissue as tissue column")
        
        # Set sample ID as index
        merged_attrs.set_index('sample_id', inplace=True)
        
        logger.info(f"Loaded metadata for {len(merged_attrs)} GTEx samples")
        
        # Extract tissue statistics
        if 'tissue' in merged_attrs.columns:
            tissue_counts = merged_attrs['tissue'].value_counts()
            logger.info(f"Found {len(tissue_counts)} unique tissue types")
            logger.info(f"Top 5 tissues by sample count:")
            for tissue, count in tissue_counts.head(5).items():
                logger.info(f"  - {tissue}: {count} samples")
        
        return merged_attrs
    
    except Exception as e:
        logger.error(f"Error loading GTEx metadata: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return pd.DataFrame()


def load_gtex_expression(file_path):
    """
    Load GTEx expression data from GCT file - memory optimized.
    
    Parameters:
    -----------
    file_path : str
        Path to the GTEx GCT file
    
    Returns:
    --------
    pandas.DataFrame
        Expression data with genes as rows and samples as columns
    """
    logger.info(f"Loading GTEx expression data from {file_path}")
    
    try:
        # First, read header to get dimensions
        with gzip.open(file_path, 'rt') as f:
            # Skip version line
            next(f)
            # Read dimensions line
            dims = next(f).strip().split()
            n_genes, n_samples = int(dims[0]), int(dims[1])
            logger.info(f"GTEx file dimensions: {n_genes} genes, {n_samples} samples")
            
            # Read column headers
            header_line = next(f).strip().split('\t')
            sample_ids = header_line[2:]  # Skip Name and Description columns
        
        # Now read the data in chunks
        chunk_size = 1000  # Adjust based on memory constraints
        
        logger.info(f"Reading GTEx data in chunks of {chunk_size} genes")
        
        # Initialize the gene list and expression list
        all_genes = []
        all_data = []
        
        with gzip.open(file_path, 'rt') as f:
            # Skip header lines
            next(f)  # Version
            next(f)  # Dimensions
            next(f)  # Column headers
            
            # Read data line by line
            line_count = 0
            for line in f:
                line_count += 1
                
                fields = line.strip().split('\t')
                gene_id = fields[0]  # First column is the gene ID
                # Skip Description column (position 1)
                expression = [float(x) for x in fields[2:]]  # Expression values start at position 2
                
                all_genes.append(gene_id)
                all_data.append(expression)
                
                # Log progress
                if line_count % 10000 == 0:
                    logger.info(f"Processed {line_count} genes")
        
        # Create the DataFrame
        data_df = pd.DataFrame(all_data, index=all_genes, columns=sample_ids)
        logger.info(f"Loaded expression data with {data_df.shape[0]} genes and {data_df.shape[1]} samples")
        
        return data_df
    
    except Exception as e:
        logger.error(f"Error loading GTEx expression data: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return pd.DataFrame()


def process_gtex_data(input_file, output_file):
    """
    Process GTEx RNA-seq data into standardized AnnData - enhanced version with improved gene mapping.
    
    Parameters:
    -----------
    input_file : str
        Path to GTEx expression file (GCT format)
    output_file : str
        Path to save the standardized AnnData
    
    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    logger.info(f"Processing GTEx data from {input_file}")
    
    # Load expression data
    expr_df = load_gtex_expression(input_file)
    if expr_df.empty:
        logger.error("Failed to load GTEx expression data")
        return None
    
    # Load metadata
    metadata_df = load_gtex_metadata()
    if metadata_df.empty:
        logger.warning("No GTEx metadata found, creating minimal metadata")
        metadata_df = pd.DataFrame(index=expr_df.columns)
        metadata_df['sample_id'] = metadata_df.index
    
    # Find common samples between expression data and metadata
    common_samples = sorted(set(expr_df.columns).intersection(set(metadata_df.index)))
    
    if not common_samples:
        logger.warning("No samples in common between expression data and metadata")
        # Create minimal metadata for all samples
        metadata_df = pd.DataFrame(index=expr_df.columns)
        common_samples = list(expr_df.columns)
    else:
        logger.info(f"Found {len(common_samples)} samples with both expression data and metadata")
        # Filter expression data to only include samples with metadata
        expr_df = expr_df[common_samples]
        # Filter metadata to only include samples with expression data
        metadata_df = metadata_df.loc[common_samples]
    
    # Standardize gene IDs
    logger.info(f"Standardizing {expr_df.shape[0]} gene IDs")
    
    # Map original IDs to standardized IDs
    gene_id_mapping = {gene_id: standardize_ensembl_id(gene_id) for gene_id in expr_df.index}
    
    # Get unique standardized IDs
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")
    
    # Create a mapping from std_ids to row indices
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    
    # Create the expression matrix
    sample_count = len(common_samples)
    gene_count = len(unique_std_ids)
    expr_matrix = np.zeros((gene_count, sample_count), dtype=np.float32)
    
    # Fill the matrix
    for sample_idx, sample_id in enumerate(common_samples):
        if sample_idx % 1000 == 0 or sample_idx == len(common_samples) - 1:
            logger.info(f"Processing sample {sample_idx + 1}/{len(common_samples)}")
        
        # Process genes for this sample
        for gene_id, expr_val in zip(expr_df.index, expr_df[sample_id]):
            std_id = gene_id_mapping[gene_id]
            std_idx = std_id_to_idx[std_id]
            
            # If multiple original IDs map to same standardized ID, use maximum value
            expr_matrix[std_idx, sample_idx] = max(expr_matrix[std_idx, sample_idx], expr_val)
    
    # Create DataFrame with the results
    std_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=common_samples)
    
    # Free memory
    del expr_df
    del expr_matrix
    
    # Create variable (gene) metadata DataFrame
    var_df = pd.DataFrame(index=unique_std_ids)
    var_df['gene_id'] = var_df.index
    
    # Add original ID mapping information
    var_df['original_ids'] = ""
    
    # Create a reverse mapping from std_id to original ids
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)
    
    # Add the original IDs as a concatenated string
    for std_id in var_df.index:
        if std_id in reverse_mapping:
            var_df.loc[std_id, 'original_ids'] = ";".join(reverse_mapping[std_id])
    
    # Add GENCODE annotations
    logger.info("Adding GENCODE annotations")
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)
    
    # Standardize observation metadata
    logger.info("Standardizing metadata")
    obs_df = standardize_metadata(metadata_df, 'GTEx')
    
    # Add missing standard fields
    if 'data_type' not in obs_df.columns:
        obs_df['data_type'] = 'RNA-seq'
    if 'expression_unit' not in obs_df.columns:
        obs_df['expression_unit'] = 'TPM'
    
    # Dataset info for uns
    dataset_info = {
        'source': 'GTEx',
        'version': 'v10',  # Updated to v10
        'gencode_version': 24,
        'data_type': 'RNA-seq',
        'expression_unit': 'TPM',
        'samples': len(obs_df),
        'genes': len(var_df),
        'tissue_count': len(obs_df['tissue'].unique()) if 'tissue' in obs_df.columns else 0,
        'subject_count': len(obs_df['subject_id'].unique()) if 'subject_id' in obs_df.columns else 0
    }
    
    # Create AnnData object
    logger.info("Creating AnnData object")
    adata = create_standard_anndata(std_expr_df.T, obs_df, var_df, dataset_info)
    
    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed GTEx data with {adata.n_obs} samples and {adata.n_vars} genes")
    
    return adata

# ====================================================================================
# MAGE-specific Processing
# ====================================================================================

def process_mage_dir(file_path, donor_id, tissue, sample_dfs, sample_metadata, all_gene_ids):
    """
    Process a single MAGE expression file with improved handling for file format.
    Also handles column names properly per standardized format.
    """
    try:
        file_name = os.path.basename(file_path)
        sample_id = f"{donor_id}_{os.path.splitext(file_name)[0]}"
        
        logger.debug(f"Processing MAGE file: {file_path}")
        
        # Check if file exists
        if not os.path.exists(file_path):
            logger.warning(f"File does not exist: {file_path}")
            return False
            
        # Check file size
        file_size = os.path.getsize(file_path)
        if file_size == 0:
            logger.warning(f"Empty file: {file_path}")
            return False
            
        # Try to read the file - MAGE files are tab-separated even though they have .csv extension
        try:
            # First check the delimiter by reading a few lines
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                
            # Determine delimiter based on content
            if '\t' in first_line:
                delimiter = '\t'
            elif ',' in first_line:
                delimiter = ','
            else:
                # Default to tab if can't determine
                delimiter = '\t'
                
            logger.debug(f"Using delimiter '{delimiter}' for file {file_path}")
            
            # Now read with the appropriate delimiter
            df = pd.read_csv(file_path, delimiter=delimiter)
            
            # Check if empty after reading
            if df.empty:
                logger.warning(f"Empty data after reading: {file_path}")
                return False
                
            logger.debug(f"Read data with columns: {df.columns.tolist()}")
            
            # Check for expected columns
            if 'gene_id' in df.columns and 'TPM' in df.columns:
                # File has the expected structure, but TPM is actually normalized microarray intensity
                # Rename TPM to normalized_intensity if processing ADNI data
                if 'ADNI' in file_path:
                    df = df.rename(columns={'TPM': 'normalized_intensity'})
                    expr_col = 'normalized_intensity'
                    expr_type = 'normalized_intensity'
                else:  # For other datasets like MAGE, keep TPM
                    expr_col = 'TPM'
                    expr_type = 'TPM'
                
                # Remove FPKM column if it exists (as per your requirement)
                if 'FPKM' in df.columns:
                    df = df.drop('FPKM', axis=1)
                
                # Create simplified DataFrame with gene_id and expression values
                simple_df = pd.DataFrame()
                simple_df[expr_col] = df[expr_col]
                simple_df.index = df['gene_id']
                
                # Store the data
                sample_dfs[sample_id] = simple_df
                
                # Collect gene IDs
                all_gene_ids.update(df['gene_id'])
                
                # Check if gene IDs are Ensembl IDs
                is_ensembl = all(str(gene_id).startswith("ENSG") for gene_id in df['gene_id'].iloc[:5]) if len(df) >= 5 else False
                
                # Prepare metadata
                metadata = {
                    'sample_id': sample_id,
                    'donor_id': donor_id,
                    'subject_id': donor_id,
                    'tissue': tissue,
                    'dataset': 'MAGE',
                    'data_type': 'RNA-seq',
                    'expression_unit': expr_type,  # This will be TPM for RNA-seq or normalized_intensity for microarray
                    'is_gencode': 'gencode' in file_name.lower(),
                    'is_ensembl': is_ensembl
                }
                
                # Store metadata
                sample_metadata[sample_id] = metadata
                logger.debug(f"Successfully processed {file_path} with {len(simple_df)} genes")
                return True
            else:
                # Unexpected columns - try to identify gene ID and expression columns
                logger.warning(f"Unexpected columns in {file_path}: {df.columns.tolist()}")
                
                # Look for gene ID column - try common names
                gene_id_col = None
                for col in df.columns:
                    if 'gene' in col.lower() and 'id' in col.lower():
                        gene_id_col = col
                        break
                
                # Look for expression column - try common names
                expr_col = None
                for col in df.columns:
                    if col.upper() in ['TPM', 'FPKM', 'COUNTS', 'EXPRESSION', 'NORMALIZED_INTENSITY']:
                        expr_col = col
                        break
                
                # If not found, try to use first and second columns as fallback
                if gene_id_col is None and len(df.columns) >= 1:
                    gene_id_col = df.columns[0]
                    logger.debug(f"Using first column '{gene_id_col}' as gene ID column")
                
                if expr_col is None and len(df.columns) >= 2:
                    expr_col = df.columns[1]
                    logger.debug(f"Using second column '{expr_col}' as expression column")
                
                # If we have identified columns, process the file
                if gene_id_col is not None and expr_col is not None:
                    # Create simplified DataFrame
                    simple_df = pd.DataFrame()
                    simple_df[expr_col] = df[expr_col]
                    simple_df.index = df[gene_id_col]
                    
                    # Store the data
                    sample_dfs[sample_id] = simple_df
                    
                    # Collect gene IDs
                    all_gene_ids.update(df[gene_id_col])
                    
                    # Check if gene IDs are Ensembl IDs
                    is_ensembl = all(str(gene_id).startswith("ENSG") for gene_id in df[gene_id_col].iloc[:5]) if len(df) >= 5 else False
                    
                    # Determine expression type based on column name
                    if 'TPM' in expr_col.upper():
                        expr_type = "TPM"
                    elif 'FPKM' in expr_col.upper():
                        expr_type = "FPKM"
                    elif 'COUNT' in expr_col.upper():
                        expr_type = "counts"
                    elif 'NORMALIZED' in expr_col.upper() or 'INTENSITY' in expr_col.upper():
                        expr_type = "normalized_intensity"
                    else:
                        expr_type = "TPM"  # Default
                    
                    # Prepare metadata
                    metadata = {
                        'sample_id': sample_id,
                        'donor_id': donor_id,
                        'subject_id': donor_id,
                        'tissue': tissue,
                        'dataset': 'MAGE',
                        'data_type': 'RNA-seq',
                        'expression_unit': expr_type,
                        'is_gencode': 'gencode' in file_name.lower(),
                        'is_ensembl': is_ensembl
                    }
                    
                    # Store metadata
                    sample_metadata[sample_id] = metadata
                    logger.debug(f"Successfully processed {file_path} with {len(simple_df)} genes")
                    return True
                else:
                    logger.warning(f"Could not identify required columns in {file_path}")
                    return False
        
        except Exception as e:
            logger.warning(f"Error reading file {file_path}: {e}")
            return False
    
    except Exception as e:
        logger.error(f"Error processing MAGE file {file_path}: {e}")
        return False

def process_mage_data(input_dir, output_file):
    """
    Process MAGE RNA-seq data into standardized AnnData - updated to prioritize original files.
    """
    logger.info(f"Processing MAGE data from {input_dir}")
    
    # Track processing statistics
    processed_donors = 0
    processed_samples = 0
    
    # Create collections to store data
    sample_dfs = {}  # Data frames for each sample
    sample_metadata = {}  # Metadata for each sample
    all_gene_ids = set()  # Track all gene IDs
    
    # Find all donor directories - includes both NA* and HG* patterns
    donor_dirs = glob.glob(os.path.join(input_dir, "NA*")) + glob.glob(os.path.join(input_dir, "HG*"))
    
    if not donor_dirs:
        logger.error(f"No donor directories found in {input_dir}")
        return None
    
    logger.info(f"Found {len(donor_dirs)} donor directories")
    
    # Process each donor directory
    for donor_dir in donor_dirs:
        donor_id = os.path.basename(donor_dir)
        processed_donors += 1
        
        logger.info(f"Processing donor {donor_id}")
        
        # Check for lymphoblast directory which contains the actual CSV files
        lymphoblast_dir = os.path.join(donor_dir, "lymphoblast")
        if os.path.isdir(lymphoblast_dir):
            tissue_name = "lymphoblast"
            
            # Find all CSV files in the directory
            csv_files = [f for f in os.listdir(lymphoblast_dir) if f.endswith('.csv')]
            
            # Separate original and gencode files
            original_files = [f for f in csv_files if "_gencode_" not in f]
            gencode_files = [f for f in csv_files if "_gencode_" in f]
            
            # Process original files first (prioritize them)
            for file_name in original_files:
                file_path = os.path.join(lymphoblast_dir, file_name)
                if process_mage_dir(file_path, donor_id, tissue_name, sample_dfs, 
                                    sample_metadata, all_gene_ids):
                    processed_samples += 1
                    
            # Only process gencode files if no original files were successfully processed
            if not any(donor_id in sample_id for sample_id in sample_dfs.keys()):
                for file_name in gencode_files:
                    file_path = os.path.join(lymphoblast_dir, file_name)
                    if process_mage_dir(file_path, donor_id, tissue_name, sample_dfs, 
                                        sample_metadata, all_gene_ids):
                        processed_samples += 1
        else:
            # Try to find any tissue directories
            tissue_dirs = [d for d in glob.glob(os.path.join(donor_dir, "*")) 
                          if os.path.isdir(d)]
            
            for tissue_dir in tissue_dirs:
                tissue_name = os.path.basename(tissue_dir)
                
                # Find all CSV files in the directory
                csv_files = [f for f in os.listdir(tissue_dir) if f.endswith('.csv')]
                
                # Separate original and gencode files
                original_files = [f for f in csv_files if "_gencode_" not in f]
                gencode_files = [f for f in csv_files if "_gencode_" in f]
                
                # Process original files first
                for file_name in original_files:
                    file_path = os.path.join(tissue_dir, file_name)
                    if process_mage_dir(file_path, donor_id, tissue_name, sample_dfs, 
                                        sample_metadata, all_gene_ids):
                        processed_samples += 1
                
                # Only process gencode files if no original files were successfully processed for this tissue
                tissue_sample_ids = [sample_id for sample_id in sample_dfs.keys() 
                                   if donor_id in sample_id and tissue_name in sample_id]
                if not tissue_sample_ids:
                    for file_name in gencode_files:
                        file_path = os.path.join(tissue_dir, file_name)
                        if process_mage_dir(file_path, donor_id, tissue_name, sample_dfs, 
                                            sample_metadata, all_gene_ids):
                            processed_samples += 1
    
    if not sample_dfs:
        logger.error("No valid expression data files found")
        return None
    
    logger.info(f"Processed {processed_donors} donors with {processed_samples} samples")
        
    # Standardize gene IDs and create expression matrix
    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs")
    
    # Map original IDs to standardized IDs
    gene_id_mapping = {gene_id: standardize_ensembl_id(gene_id) for gene_id in all_gene_ids}
    
    # Get unique standardized IDs
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")
    
    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    sample_ids = list(sample_dfs.keys())
    sample_id_to_idx = {id: idx for idx, id in enumerate(sample_ids)}
    
    # Pre-allocate the expression matrix
    num_genes = len(unique_std_ids)
    num_samples = len(sample_ids)
    expr_matrix = np.zeros((num_genes, num_samples), dtype=np.float32)
    
    # Fill the matrix efficiently
    logger.info("Creating unified expression matrix")
    for sample_id, expr_df in sample_dfs.items():
        sample_idx = sample_id_to_idx[sample_id]
        
        # Check if the dataframe is valid
        if expr_df.empty or len(expr_df.columns) == 0:
            logger.warning(f"Skipping empty dataframe for sample {sample_id}")
            continue
            
        # Create temporary dictionary to collect values for this sample
        sample_data = {}
        
        try:
            # Get the first column name
            first_col = expr_df.columns[0]
            
            # Process each gene in this sample
            for gene_id, expr_val in zip(expr_df.index, expr_df[first_col]):
                if gene_id in gene_id_mapping:
                    std_id = gene_id_mapping[gene_id]
                    
                    # If multiple original IDs map to same standardized ID, use maximum value
                    if std_id in sample_data:
                        sample_data[std_id] = max(sample_data[std_id], expr_val)
                    else:
                        sample_data[std_id] = expr_val
            
            # Update the matrix
            for std_id, value in sample_data.items():
                if std_id in std_id_to_idx:
                    gene_idx = std_id_to_idx[std_id]
                    expr_matrix[gene_idx, sample_idx] = value
        except Exception as e:
            logger.error(f"Error processing sample {sample_id}: {e}")
            
    
    # Create the DataFrame from the matrix
    unified_df = pd.DataFrame(
        expr_matrix, 
        index=unique_std_ids, 
        columns=sample_ids
    )
    
    # Create observation metadata DataFrame
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            metadata_copy = sample_metadata[sample_id].copy()
            metadata_copy['_sample_id'] = sample_id
            obs_data.append(metadata_copy)
        else:
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({'_sample_id': sample_id, 'dataset': 'MAGE'})
    
    # Create obs_df with explicit index
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index('_sample_id', inplace=True)
    
    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]
    
    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df['gene_id'] = var_df.index
    
    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)
    
    # Add original IDs to variable metadata
    var_df['original_ids'] = var_df.index.map(lambda x: ';'.join(reverse_mapping.get(x, [])))
    
    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)
    
    # Standardize metadata
    obs_df = standardize_metadata(obs_df, 'MAGE')
    
    # Dataset info for uns
    dataset_info = {
        'source': 'MAGE',
        'gencode_version': 24,  # Assuming same version as other datasets
        'data_type': 'RNA-seq',
        'expression_unit': 'TPM',  # Update if different
        'samples': len(obs_df),
        'genes': len(var_df),
        'donor_count': processed_donors,
        'tissue_count': len(set(obs_df['tissue'] if 'tissue' in obs_df.columns else []))
    }
    
    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)
    
    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed MAGE data with {adata.n_obs} samples and {adata.n_vars} genes")
    
    return adata

# ====================================================================================
# ADNI-specific Processing
# ====================================================================================

def process_adni_data(input_dir, output_file, dict_file=None):
    """
    Process ADNI gene expression data into standardized AnnData.
    Updated to handle individual sample files with Ensembl IDs.
    """
    logger.info(f"Processing ADNI data from {input_dir}")
    
    # Track processing statistics
    processed_subjects = 0
    processed_samples = 0
    
    # Create collections to store data
    sample_dfs = {}  # Data frames for each sample
    sample_metadata = {}  # Metadata for each sample
    all_gene_ids = set()  # Track all gene IDs
    
    # Find all subject directories
    subject_dirs = glob.glob(os.path.join(input_dir, "*_S_*"))
    
    if not subject_dirs:
        logger.error(f"No subject directories found in {input_dir}")
        return None
    
    logger.info(f"Found {len(subject_dirs)} subject directories")
    
    # Process each subject directory
    for subject_dir in subject_dirs:
        subject_id = os.path.basename(subject_dir)
        processed_subjects += 1
        
        logger.info(f"Processing subject {subject_id}")
        
        # Find all CSV files in the directory
        csv_files = glob.glob(os.path.join(subject_dir, "*.csv"))
        
        if not csv_files:
            logger.warning(f"No CSV files found for subject {subject_id}")
            continue
        
        # Process each CSV file
        for file_path in csv_files:
            file_name = os.path.basename(file_path)
            sample_id = f"{subject_id}_{os.path.splitext(file_name)[0]}"
            
            logger.debug(f"Processing file: {file_path}")
            
            try:
                # Read the CSV file
                df = pd.read_csv(file_path, sep='\t')
                
                # Check for required columns
                if 'gene_id' in df.columns and 'TPM' in df.columns:
                    # Create simplified DataFrame with gene_id and TPM values
                    expr_df = pd.DataFrame()
                    expr_df['TPM'] = df['TPM']
                    expr_df.index = df['gene_id']
                    
                    # Store the data
                    sample_dfs[sample_id] = expr_df
                    
                    # Collect gene IDs
                    all_gene_ids.update(df['gene_id'])
                    
                    # Extract metadata from filename and directory structure
                    # For ADNI, we'll extract subject ID from the directory name
                    metadata = {
                        'sample_id': sample_id,
                        'subject_id': subject_id,
                        'dataset': 'ADNI',
                        'data_type': 'RNA-seq',  # These appear to be TPM values from RNA-seq
                        'expression_unit': 'TPM',
                        'tissue': 'blood',  # Default, update if tissue information is available
                        'processing': 'gencode_v24' if 'gencode_v24' in file_name else 'unknown'
                    }
                    
                    # Store metadata
                    sample_metadata[sample_id] = metadata
                    processed_samples += 1
                    
                    logger.debug(f"Successfully processed {file_path} with {len(expr_df)} genes")
                else:
                    logger.warning(f"Missing required columns in {file_path}")
            
            except Exception as e:
                logger.error(f"Error processing file {file_path}: {e}")
                continue
    
    if not sample_dfs:
        logger.error("No valid expression data files found")
        return None
    
    logger.info(f"Processed {processed_subjects} subjects with {processed_samples} samples")
    
    # Standardize gene IDs and create expression matrix
    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs")
    
    # Map original IDs to standardized IDs (removing version numbers)
    gene_id_mapping = {gene_id: standardize_ensembl_id(gene_id) for gene_id in all_gene_ids}
    
    # Get unique standardized IDs
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")
    
    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    sample_ids = list(sample_dfs.keys())
    sample_id_to_idx = {id: idx for idx, id in enumerate(sample_ids)}
    
    # Pre-allocate the expression matrix
    num_genes = len(unique_std_ids)
    num_samples = len(sample_ids)
    expr_matrix = np.zeros((num_genes, num_samples), dtype=np.float32)
    
    # Fill the matrix efficiently
    logger.info("Creating unified expression matrix")
    for sample_id, expr_df in sample_dfs.items():
        sample_idx = sample_id_to_idx[sample_id]
        
        # Check if the dataframe is valid
        if expr_df.empty or len(expr_df.columns) == 0:
            logger.warning(f"Skipping empty dataframe for sample {sample_id}")
            continue
            
        # Create temporary dictionary to collect values for this sample
        sample_data = {}
        
        try:
            # Get the first column name (TPM)
            first_col = expr_df.columns[0]
            
            # Process each gene in this sample
            for gene_id, expr_val in zip(expr_df.index, expr_df[first_col]):
                if gene_id in gene_id_mapping:
                    std_id = gene_id_mapping[gene_id]
                    
                    # If multiple original IDs map to same standardized ID, use maximum value
                    if std_id in sample_data:
                        sample_data[std_id] = max(sample_data[std_id], expr_val)
                    else:
                        sample_data[std_id] = expr_val
            
            # Update the matrix
            for std_id, value in sample_data.items():
                if std_id in std_id_to_idx:
                    gene_idx = std_id_to_idx[std_id]
                    expr_matrix[gene_idx, sample_idx] = value
        except Exception as e:
            logger.error(f"Error processing sample {sample_id}: {e}")
            
    
    # Create the DataFrame from the matrix
    unified_df = pd.DataFrame(
        expr_matrix, 
        index=unique_std_ids, 
        columns=sample_ids
    )
    
    # Create observation metadata DataFrame
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            metadata_copy = sample_metadata[sample_id].copy()
            metadata_copy['_sample_id'] = sample_id
            obs_data.append(metadata_copy)
        else:
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({'_sample_id': sample_id, 'dataset': 'ADNI'})
    
    # Create obs_df with explicit index
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index('_sample_id', inplace=True)
    
    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]
    
    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df['gene_id'] = var_df.index
    
    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)
    
    # Add original IDs to variable metadata
    var_df['original_ids'] = var_df.index.map(lambda x: ';'.join(reverse_mapping.get(x, [])))
    
    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)
    
    # Standardize metadata
    obs_df = standardize_metadata(obs_df, 'ADNI')
    
    # Dataset info for uns
    dataset_info = {
        'source': 'ADNI',
        'data_type': 'RNA-seq',  # Updated from microarray since these appear to be RNA-seq TPM values
        'gencode_version': 24,
        'expression_unit': 'TPM',
        'samples': len(obs_df),
        'genes': len(var_df),
        'subject_count': processed_subjects
    }
    
    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)
    
    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed ADNI data with {adata.n_obs} samples and {adata.n_vars} genes")
    
    return adata

# ====================================================================================
# Main Processing Pipeline
# ====================================================================================
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Multi-Dataset Standardization Pipeline')
    parser.add_argument('--encode-dir', help='Directory containing ENCODE cell line TPM files')
    parser.add_argument('--encode-entex-dir', help='Directory containing ENCODE ENTEx TPM files')
    parser.add_argument('--gtex-file', help='Path to GTEx expression file (GCT format)')
    parser.add_argument('--mage-dir', help='Directory containing MAGE expression files')
    parser.add_argument('--adni-dir', help='Directory containing ADNI sample directories with expression files')
    parser.add_argument('--output-dir', default=str(DEFAULT_OUTPUT_DIR), help='Output directory for standardized files')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    start_time = time.time()
    
    # Process each dataset if input is provided
    if args.encode_dir:
        encode_output = os.path.join(args.output_dir, 'encode_standardized.h5ad')
        logger.info(f"Processing ENCODE data from {args.encode_dir}")
        encode_data = process_encode_data(args.encode_dir, args.encode_entex_dir, encode_output)
        if encode_data is not None:
            logger.info(f"Successfully processed ENCODE data: {encode_data.n_obs} samples and {encode_data.n_vars} genes")
            # Check gene statistics
            logger.info("ENCODE gene statistics:")
            gencode_stats = encode_data.var['mapping_source'].value_counts()
            for source, count in gencode_stats.items():
                logger.info(f"  {source}: {count} genes ({count/encode_data.n_vars:.1%})")
    
    
    if args.gtex_file:
        gtex_output = os.path.join(args.output_dir, 'gtex_standardized.h5ad')
        logger.info(f"Processing GTEx data from {args.gtex_file}")
        gtex_data = process_gtex_data(args.gtex_file, gtex_output)
        if gtex_data is not None:
            logger.info(f"Successfully processed GTEx data: {gtex_data.n_obs} samples, {gtex_data.n_vars} genes")
            # Check gene statistics
            logger.info("GTEx gene statistics:")
            gencode_stats = gtex_data.var['mapping_source'].value_counts()
            for source, count in gencode_stats.items():
                logger.info(f"  {source}: {count} genes ({count/gtex_data.n_vars:.1%})")
        
    if args.mage_dir:
        mage_output = os.path.join(args.output_dir, 'mage_standardized.h5ad')
        logger.info(f"Processing MAGE data from {args.mage_dir}")
        mage_data = process_mage_data(args.mage_dir, mage_output)
        if mage_data is not None:
            logger.info(f"Successfully processed MAGE data: {mage_data.n_obs} samples, {mage_data.n_vars} genes")
    
    # Process ADNI data
    if args.adni_dir:
        adni_output = os.path.join(args.output_dir, 'adni_standardized.h5ad')
        logger.info(f"Processing ADNI data from directory: {args.adni_dir}")
        adni_data = process_adni_data(args.adni_dir, adni_output)
        if adni_data is not None:
            logger.info(f"Successfully processed ADNI data: {adni_data.n_obs} samples, {adni_data.n_vars} genes")
            
    # Check for common genes if multiple datasets were processed
    processed_datasets = []
    if args.encode_dir and 'encode_data' in locals() and encode_data is not None:
        processed_datasets.append(('ENCODE', encode_data))
    if args.gtex_file and 'gtex_data' in locals() and gtex_data is not None:
        processed_datasets.append(('GTEx', gtex_data))
    if args.mage_dir and 'mage_data' in locals() and mage_data is not None:
        processed_datasets.append(('MAGE', mage_data))
    if args.adni_dir and 'adni_data' in locals() and adni_data is not None:
        processed_datasets.append(('ADNI', adni_data))
    
    if len(processed_datasets) > 1:
        logger.info("Analyzing gene overlap between datasets:")
        
        # Compare each pair of datasets
        for i in range(len(processed_datasets)):
            for j in range(i+1, len(processed_datasets)):
                name1, data1 = processed_datasets[i]
                name2, data2 = processed_datasets[j]
                
                genes1 = set(data1.var_names)
                genes2 = set(data2.var_names)
                common_genes = genes1.intersection(genes2)
                
                logger.info(f"  {name1} ({len(genes1)} genes) vs {name2} ({len(genes2)} genes):")
                logger.info(f"    Common genes: {len(common_genes)}")
                logger.info(f"    Overlap: {len(common_genes)/len(genes1):.1%} of {name1}, {len(common_genes)/len(genes2):.1%} of {name2}")
                
                # Check overlap of highly expressed genes
                # Define highly expressed as genes in the top 10% by mean expression
                def get_top_genes(adata, percent=10):
                    mean_expr = adata.X.mean(axis=0)
                    cutoff = np.percentile(mean_expr, 100 - percent)
                    top_indices = np.where(mean_expr >= cutoff)[0]
                    return set(adata.var_names[top_indices])
                
                top_genes1 = get_top_genes(data1)
                top_genes2 = get_top_genes(data2)
                common_top = top_genes1.intersection(top_genes2)
                
                logger.info(f"    Top 10% expressed genes: {name1} ({len(top_genes1)}), {name2} ({len(top_genes2)})")
                logger.info(f"    Common top genes: {len(common_top)}")
                logger.info(f"    Top gene overlap: {len(common_top)/len(top_genes1):.1%} of top {name1}, {len(common_top)/len(top_genes2):.1%} of top {name2}")
        
        # If both ENCODE and GTEx were processed, create a gene compatibility report
        if 'encode_data' in locals() and 'gtex_data' in locals() and encode_data is not None and gtex_data is not None:
            report_file = os.path.join(args.output_dir, 'gene_compatibility_report.csv')
            logger.info(f"Creating gene compatibility report: {report_file}")
            
            # Create DataFrame with gene information from both datasets
            encode_genes = pd.DataFrame(index=encode_data.var_names)
            encode_genes['gene_name'] = encode_data.var['gene_name']
            encode_genes['gene_type'] = encode_data.var['gene_type']
            encode_genes['in_encode'] = True
            
            gtex_genes = pd.DataFrame(index=gtex_data.var_names)
            gtex_genes['gene_name'] = gtex_data.var['gene_name']
            gtex_genes['gene_type'] = gtex_data.var['gene_type']
            gtex_genes['in_gtex'] = True
            
            # Merge the DataFrames
            all_genes = encode_genes.join(gtex_genes, how='outer', lsuffix='_encode', rsuffix='_gtex')
            all_genes['in_encode'] = all_genes['in_encode'].fillna(False)
            all_genes['in_gtex'] = all_genes['in_gtex'].fillna(False)
            
            # Fill in missing gene names using the other dataset
            all_genes['gene_name'] = all_genes['gene_name_encode'].combine_first(all_genes['gene_name_gtex'])
            all_genes['gene_type'] = all_genes['gene_type_encode'].combine_first(all_genes['gene_type_gtex'])
            
            # Calculate average expression in each dataset
            all_genes['mean_expr_encode'] = 0.0
            all_genes['mean_expr_gtex'] = 0.0
            
            for gene in all_genes.index:
                if gene in encode_data.var_names:
                    gene_idx = encode_data.var_names.get_loc(gene)
                    all_genes.loc[gene, 'mean_expr_encode'] = encode_data.X[:, gene_idx].mean()
                if gene in gtex_data.var_names:
                    gene_idx = gtex_data.var_names.get_loc(gene)
                    all_genes.loc[gene, 'mean_expr_gtex'] = gtex_data.X[:, gene_idx].mean()
            
            # Save the report
            all_genes.to_csv(report_file)
            logger.info(f"Saved gene compatibility report with {len(all_genes)} genes")
    
    total_time = time.time() - start_time
    logger.info(f"Total processing time: {total_time:.2f} seconds")

if __name__ == '__main__':
    main()
                    