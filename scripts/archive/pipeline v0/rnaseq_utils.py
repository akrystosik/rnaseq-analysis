#!/usr/bin/env python3
"""
RNA-seq Data Processing Utilities

This module contains utility functions for RNA-seq data processing,
standardization, and validation used across the RNA-seq standardization pipeline.
"""

import os
import json
import logging
import pandas as pd
import numpy as np
from pathlib import Path

# Set up logging
logger = logging.getLogger('rnaseq_utils')

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
GENCODE_MAPPING_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gencode_v24_complete_mapping.csv")
METADATA_JSON_DIR = BASE_DIR / "metadata/json"

# Constants
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

# ===================================================================================
# File Loading Functions
# ===================================================================================

def load_json_mapping(filename, default=None):
    """
    Load a JSON mapping file from the metadata directory.
    
    Args:
        filename: Name of the JSON file
        default: Default value to return if file not found
        
    Returns:
        Dictionary with mapping data
    """
    filepath = METADATA_JSON_DIR / filename
    
    try:
        if filepath.exists():
            with open(filepath, 'r') as f:
                return json.load(f)
        else:
            logger.warning(f"Mapping file not found: {filepath}")
            return default if default is not None else {}
    except Exception as e:
        logger.error(f"Error loading mapping file {filepath}: {e}")
        return default if default is not None else {}

def load_mappings():
    """
    Load all standard mapping files into a single dictionary.
    
    Returns:
        Dictionary containing all mappings
    """
    mappings = {
        'TISSUE_TO_UBERON': load_json_mapping('tissue_to_uberon.json', {}),
        'ASSAY_TO_EFO': load_json_mapping('assay_to_efo.json', {}),
        'AGE_TO_HSAPDV': load_json_mapping('age_to_hsapdv.json', {}),
        'SEX_STANDARDIZATION': load_json_mapping('sex_standardization.json', {}),
        'SPECIES_TO_NCBI_TAXON': load_json_mapping('species_to_taxon.json', {})
    }
    
    # Load dataset-specific metadata
    for dataset in ['encode', 'gtex', 'mage', 'adni', 'entex']:
        key = f"{dataset.upper()}_METADATA"
        mappings[key] = load_json_mapping(f"{dataset}_metadata.json", {})
    
    return mappings

def load_gencode_mapping():
    """
    Load GENCODE mapping from CSV file.
    
    Returns:
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

# ===================================================================================
# Gene ID Processing Functions
# ===================================================================================

def standardize_ensembl_id(gene_id):
    """
    Standardize Ensembl gene ID by removing version number.
    
    Args:
        gene_id: Gene ID to standardize
        
    Returns:
        Standardized gene ID
    """
    if pd.isna(gene_id):
        return None
    
    # Convert to string if needed
    gene_id = str(gene_id)
    
    # Remove version if present
    if '.' in gene_id and gene_id.startswith('ENSG'):
        return gene_id.split('.')[0]
    return gene_id

def add_gencode_annotations(var_df, gencode_mapping):
    """
    Add GENCODE annotations to variable DataFrame.
    
    Args:
        var_df: DataFrame with gene information
        gencode_mapping: Dictionary mapping gene IDs to gene information
        
    Returns:
        DataFrame with added GENCODE annotations
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

# ===================================================================================
# Metadata Processing Functions
# ===================================================================================

def map_age_to_hsapdv(age_value, age_mapping=None):
    """
    Map an age value to the appropriate HsapDv ontology term.
    
    Args:
        age_value: Age value as a string
        age_mapping: Dictionary mapping age to HsapDv terms
        
    Returns:
        HsapDv ontology ID or empty string if mapping not found
    """
    if not age_mapping:
        age_mapping = load_json_mapping('age_to_hsapdv.json', {})
        
    if not age_value or age_value == '' or pd.isna(age_value):
        return ""
    
    # Check if age is already a developmental stage
    if age_value.lower() in age_mapping:
        return age_mapping[age_value.lower()]
    
    # Try to convert to integer age
    try:
        age_years = int(age_value)
        
        # Map to age ranges
        if age_years < 2:
            return age_mapping.get("0-1", "")  # infant
        elif age_years < 13:
            return age_mapping.get("2-12", "")  # child
        elif age_years < 20:
            return age_mapping.get("13-19", "")  # adolescent
        elif age_years < 40:
            return age_mapping.get("20-39", "")  # young adult
        elif age_years < 60:
            return age_mapping.get("40-59", "")  # middle aged adult
        else:
            return age_mapping.get("60+", "")  # elderly
    except ValueError:
        # Not a simple integer - check for ranges or other formats
        # This could be expanded with more pattern matching
        return ""


def map_tissue_to_ontology(tissue_name, mappings):
    """
    Map a tissue name to an ontology ID.
    
    Args:
        tissue_name: Name of the tissue
        mappings: Dictionary containing mapping information
        
    Returns:
        Tuple of (ontology_id, confidence)
    """
    if not tissue_name or pd.isna(tissue_name):
        return "", "none"
    
    # Normalize tissue name
    tissue_lower = str(tissue_name).lower().strip()
    
    # Get tissue mappings from the mappings dictionary
    tissue_to_uberon = mappings.get('TISSUE_TO_UBERON', {})
    
    # Check for direct match
    if tissue_lower in tissue_to_uberon:
        return tissue_to_uberon[tissue_lower], "high"
    
    # Check for case-insensitive match
    for known_tissue, ontology_id in tissue_to_uberon.items():
        if known_tissue.lower() == tissue_lower:
            return ontology_id, "high"
    
    # Try substring matching for approximate matches
    for known_tissue, ontology_id in tissue_to_uberon.items():
        if known_tissue.lower() in tissue_lower or tissue_lower in known_tissue.lower():
            return ontology_id, "medium"
    
    # No match found
    return "", "none"

def standardize_metadata(metadata_df, dataset_name, mappings=None):
    """
    Standardize metadata across datasets to ensure consistent fields and values.
    Maps fields to standard ontology terms where possible.
    
    Args:
        metadata_df: Metadata DataFrame
        dataset_name: Source dataset name
        mappings: Dictionary with mapping information
        
    Returns:
        Standardized metadata DataFrame with ontology mappings
    """
    # Load mappings if not provided
    if mappings is None:
        mappings = load_mappings()
    
    # Make a copy to avoid modifying the original
    df = metadata_df.copy()
    
    # Get specific mappings
    tissue_to_uberon = mappings.get('TISSUE_TO_UBERON', {})
    assay_to_efo = mappings.get('ASSAY_TO_EFO', {})
    sex_standardization = mappings.get('SEX_STANDARDIZATION', {})
    species_to_taxon = mappings.get('SPECIES_TO_NCBI_TAXON', {})
    
    # Add dataset identifier
    df['dataset'] = dataset_name
    
    # Ensure all core fields exist
    for field in CORE_METADATA_FIELDS:
        if field not in df.columns:
            df[field] = ""
    
    # Add species field if not present - default to human
    if 'species' not in df.columns:
        df['species'] = "human"
        df['species_ontology'] = species_to_taxon.get("human", "")
    else:
        # Map species to NCBI Taxon
        df['species_ontology'] = df['species'].map(
            lambda x: species_to_taxon.get(x, species_to_taxon.get("unknown", ""))
        )
    
    # Standardize sex values
    if 'sex' in df.columns:
        df['sex'] = df['sex'].map(lambda x: sex_standardization.get(x, "unknown"))
    
    # Standardize tissue to Uberon ontology
    if 'tissue' in df.columns:
        # Keep original tissue names
        df['tissue_original'] = df['tissue']
        # Normalize tissue names (lowercase, strip whitespace)
        normalized_tissue = df['tissue'].astype(str).str.lower().str.strip()
        
        # Map to standardized tissue names
        df['tissue'] = normalized_tissue.map(
            lambda x: x if pd.isna(x) else x
        )
        
        # Map to Uberon IDs
        df['tissue_ontology'] = normalized_tissue.map(
            lambda x: tissue_to_uberon.get(x, "")
        )
        
        # For missing mappings, try a more flexible approach
        mask = df['tissue_ontology'] == ""
        if mask.any():
            for tissue_name, uberon_id in tissue_to_uberon.items():
                # Find tissues that contain this tissue name as a substring
                substr_mask = normalized_tissue.str.contains(tissue_name.lower(), na=False)
                # Update only those that don't already have a mapping
                update_mask = mask & substr_mask
                if update_mask.any():
                    df.loc[update_mask, 'tissue_ontology'] = uberon_id
    
    # Standardize data_type to EFO ontology
    if 'data_type' in df.columns:
        # Keep original data type
        df['data_type_original'] = df['data_type']
        
        # Map to EFO IDs
        df['assay_ontology'] = df['data_type'].map(
            lambda x: assay_to_efo.get(x, "")
        )
        
        # Also try to map extraction_method if available
        if 'extraction_method' in df.columns:
            # For rows without assay_ontology, try mapping from extraction_method
            mask = df['assay_ontology'] == ""
            df.loc[mask, 'assay_ontology'] = df.loc[mask, 'extraction_method'].map(
                lambda x: assay_to_efo.get(x, "")
            )

    # Handle missing subject_id by using donor_id as fallback
    if 'subject_id' not in df.columns and 'donor_id' in df.columns:
        logger.info(f"Using donor_id as subject_id for {dataset_name}")
        df['subject_id'] = df['donor_id']
    elif 'subject_id' not in df.columns and 'donor' in df.columns:
        logger.info(f"Using donor as subject_id for {dataset_name}")
        df['subject_id'] = df['donor']    
    
    
    # For unmapped tissues, add them to a pending file for curator review
    if 'tissue' in df.columns:
        unmapped_tissues = []
        for tissue in df['tissue'].dropna().unique():
            tissue_ontology, confidence = map_tissue_to_ontology(tissue, mappings)
            if not tissue_ontology and tissue not in unmapped_tissues:
                unmapped_tissues.append(tissue)
        

        if unmapped_tissues:
            logger.warning(f"{dataset_name}: Unmapped tissues: {', '.join(unmapped_tissues)}")
            # Only log the unmapped tissues, don't try to save to file as mappings_dir might not be defined
            logger.info(f"Unmapped tissues for potential mapping: {', '.join(unmapped_tissues)}")        
            # Add to pending file for curator review
            # Just log the unmapped tissues instead of trying to write to a file
            logger.info(f"Unmapped tissues for {dataset_name}: {', '.join([str(t) for t in unmapped_tissues])}")
    

    
    # Standardize age - try to map age ranges to developmental stages
    if 'age' in df.columns:
        # Keep original age values
        df['age_original'] = df['age']
        
        # Convert age to string and clean
        df['age'] = df['age'].astype(str).replace('nan', '').replace('None', '')
        
        # Map age ranges to developmental stages
        df['developmental_stage_ontology'] = df['age'].apply(
            lambda x: map_age_to_hsapdv(x, mappings.get('AGE_TO_HSAPDV', {}))
        )
    
    # Convert categorical fields to categorical data type
    categorical_fields = ['sex', 'tissue', 'dataset', 'cell_type', 'disease', 'ethnicity', 'data_type']
    for field in categorical_fields:
        if field in df.columns:
            df[field] = pd.Categorical(df[field].astype(str))
    
    return df

# ===================================================================================
# AnnData Processing Functions
# ===================================================================================

def prepare_for_anndata(df, var_df, data_df):
    """
    Prepare dataframes for AnnData creation with consistent data types.
    
    Args:
        df: Observation metadata DataFrame (samples)
        var_df: Variable metadata DataFrame (genes)
        data_df: Expression data DataFrame
        
    Returns:
        Processed observation, variable, and data DataFrames
    """
    # Make copies to avoid modifying the original DataFrames
    df = df.copy()
    var_df = var_df.copy()
    
    # Handle NA values in observation metadata
    for col in df.columns:
        if pd.api.types.is_numeric_dtype(df[col]):
            # For numeric columns, fill NaNs with 0
            df[col] = df[col].fillna(0)
        elif isinstance(df[col].dtype, pd.CategoricalDtype):
            # For categorical columns, convert to string first
            # Then back to categorical with empty string as a category
            str_series = df[col].astype(str).fillna("")
            categories = list(df[col].cat.categories)
            if "" not in categories:
                categories.append("")
            df[col] = pd.Categorical(str_series, categories=categories)
        else:
            # For other columns (object/string), fill NaNs with empty string
            df[col] = df[col].fillna("")
    
    # Handle NA values in variable metadata
    for col in var_df.columns:
        if col != 'gene_id':  # Don't convert index/gene_id
            if isinstance(var_df[col].dtype, pd.CategoricalDtype):
                # Same approach as above for categorical columns
                str_series = var_df[col].astype(str).fillna("")
                categories = list(var_df[col].cat.categories)
                if "" not in categories:
                    categories.append("")
                var_df[col] = pd.Categorical(str_series, categories=categories)
            elif pd.api.types.is_numeric_dtype(var_df[col]):
                var_df[col] = var_df[col].fillna(0)
            else:
                var_df[col] = var_df[col].fillna("")
    
    # Ensure data matrix has no NaN values
    data_df = data_df.fillna(0)
    
    return df, var_df, data_df

# ===================================================================================
# Validation Functions
# ===================================================================================

def validate_metadata(metadata_df, dataset_name=None):
    """
    Validate standardized metadata and log warnings for missing or invalid values.
    
    Args:
        metadata_df: Standardized metadata DataFrame
        dataset_name: Dataset name for logging purposes
        
    Returns:
        (validated_df, validation_report) tuple
    """
    # Make a copy to avoid modifying the original
    df = metadata_df.copy()
    
    # Initialize validation report
    report = {
        'dataset': dataset_name or 'unknown',
        'sample_count': len(df),
        'missing_fields': [],
        'missing_values': {},
        'unmapped_tissues': [],
        'unmapped_assays': [],
        'unmapped_ages': [],
        'validation_status': 'passed'
    }
    
    # Check required fields
    required_fields = ['sample_id', 'subject_id', 'sex', 'tissue', 'dataset', 'data_type']
    missing_fields = [field for field in required_fields if field not in df.columns]
    
    if missing_fields:
        logger.warning(f"{dataset_name}: Missing required fields in metadata: {', '.join(missing_fields)}")
        report['missing_fields'] = missing_fields
        report['validation_status'] = 'warning'
    
    # Check for missing values in key fields
    for field in df.columns:
        missing_count = df[field].isna().sum()
        if missing_count > 0:
            logger.warning(f"{dataset_name}: Field '{field}' has {missing_count} missing values out of {len(df)} samples")
            report['missing_values'][field] = int(missing_count)
            report['validation_status'] = 'warning'
    
    # Check ontology mappings
    if 'tissue_ontology' in df.columns:
        missing_ontology = (df['tissue_ontology'] == "").sum()
        if missing_ontology > 0:
            logger.warning(f"{dataset_name}: Missing tissue ontology mappings for {missing_ontology} samples out of {len(df)}")
            # Log the unmapped tissues
            unmapped_tissues = df.loc[df['tissue_ontology'] == "", 'tissue'].unique()
            unmapped_tissues = [t for t in unmapped_tissues if t and not pd.isna(t)]
            logger.warning(f"{dataset_name}: Unmapped tissues: {', '.join(unmapped_tissues)}")
            report['unmapped_tissues'] = list(unmapped_tissues)
            report['validation_status'] = 'warning'
    
    if 'assay_ontology' in df.columns:
        missing_ontology = (df['assay_ontology'] == "").sum()
        if missing_ontology > 0:
            logger.warning(f"{dataset_name}: Missing assay ontology mappings for {missing_ontology} samples out of {len(df)}")
            # Log the unmapped assays
            unmapped_assays = df.loc[df['assay_ontology'] == "", 'data_type'].unique()
            unmapped_assays = [a for a in unmapped_assays if a and not pd.isna(a)]
            logger.warning(f"{dataset_name}: Unmapped assays: {', '.join(unmapped_assays)}")
            report['unmapped_assays'] = list(unmapped_assays)
            report['validation_status'] = 'warning'
    
    if 'developmental_stage_ontology' in df.columns:
        missing_ontology = (df['developmental_stage_ontology'] == "").sum()
        if missing_ontology > 0 and 'age' in df.columns:
            logger.warning(f"{dataset_name}: Missing developmental stage mappings for {missing_ontology} samples out of {len(df)}")
            # Log the unmapped ages
            unmapped_ages = df.loc[df['developmental_stage_ontology'] == "", 'age'].unique()
            unmapped_ages = [a for a in unmapped_ages if a and not pd.isna(a)]
            logger.warning(f"{dataset_name}: Unmapped ages: {', '.join(unmapped_ages)}")
            report['unmapped_ages'] = list(unmapped_ages)
            report['validation_status'] = 'warning'
    
    # Add validation report to the DataFrame as metadata
    df.attrs['validation_report'] = report
    
    return df, report

def validate_dataset_metadata(metadata, dataset_name):
    """
    Validate that dataset metadata meets requirements.
    
    Args:
        metadata: Dictionary with dataset metadata
        dataset_name: Name of the dataset
        
    Returns:
        True if valid, raises ValueError otherwise
    """
    # Check required fields
    required_fields = ['harmonized_gencode_version', 'harmonized_reference_genome']
    for field in required_fields:
        if field not in metadata:
            raise ValueError(f"Dataset {dataset_name} metadata is missing required field '{field}'")
    
    # Validate GENCODE version
    harmonized_gencode = metadata['harmonized_gencode_version']
    if harmonized_gencode not in ["24", "v24"]:
        raise ValueError(f"Dataset {dataset_name} must have harmonized_gencode_version set to 'v24', got '{harmonized_gencode}'")
    
    # Validate reference genome
    harmonized_genome = metadata['harmonized_reference_genome']
    if harmonized_genome not in ["hg38", "GRCh38"]:
        raise ValueError(f"Dataset {dataset_name} must have harmonized_reference_genome set to 'hg38' or 'GRCh38', got '{harmonized_genome}'")
    
    return True

def load_dataset_specific_metadata(metadata_dir, dataset_name):
    """
    Load dataset-specific metadata from JSON file.
    
    Args:
        metadata_dir: Directory containing dataset metadata JSON files
        dataset_name: Name of the dataset to load metadata for
        
    Returns:
        Dictionary with dataset metadata or None if not found
    """
    if not metadata_dir:
        return None
    
    # Try both lowercase and original case filenames
    potential_filenames = [
        f"{dataset_name.lower()}_metadata.json",
        f"{dataset_name}_metadata.json"
    ]
    
    for filename in potential_filenames:
        metadata_path = os.path.join(metadata_dir, filename)
        if os.path.exists(metadata_path):
            try:
                logger.info(f"Loading dataset-specific metadata from {metadata_path}")
                with open(metadata_path, 'r') as f:
                    metadata = json.load(f)
                
                # Validate metadata meets requirements
                validate_dataset_metadata(metadata, dataset_name)
                
                return metadata
            except ValueError as ve:
                # Re-raise validation errors
                raise ve
            except Exception as e:
                logger.warning(f"Error loading dataset metadata from {metadata_path}: {e}")
                return None
    
    logger.info(f"No dataset-specific metadata found for {dataset_name}")
    return None

def apply_dataset_specific_metadata(adata, metadata):
    """
    Apply dataset-specific metadata from JSON to AnnData object.
    
    Args:
        adata: AnnData object to update
        metadata: Dictionary with metadata to apply
        
    Returns:
        Updated AnnData object
    """
    if not metadata:
        return adata
    
    logger.info("Applying dataset-specific metadata from JSON file")
    
    # Apply top-level attributes to uns
    for key, value in metadata.items():
        if key not in ['obs_columns', 'metadata_source']:  # Handle these separately
            if isinstance(value, dict):
                # For dictionary values like rna_seq_protocol, genome_info, etc.
                if key not in adata.uns:
                    adata.uns[key] = {}
                adata.uns[key].update(value)
            else:
                # For simple values
                adata.uns[key] = value
    
    # Apply observation column updates if present
    if 'obs_columns' in metadata:
        for col_name, value in metadata['obs_columns'].items():
            logger.info(f"  Updating obs column {col_name}")
            adata.obs[col_name] = value
    
    # Add metadata source information
    if 'metadata_source' in metadata:
        if 'metadata_sources' not in adata.uns:
            adata.uns['metadata_sources'] = []

        # adata.uns['metadata_sources'].append(metadata['metadata_source'])

        if 'metadata_sources' in adata.uns:
            # If it's a numpy array, convert to a list first
            if isinstance(adata.uns['metadata_sources'], np.ndarray):
                current_sources = adata.uns['metadata_sources'].tolist()
                current_sources.append(metadata['metadata_source'])
                adata.uns['metadata_sources'] = current_sources
            elif isinstance(adata.uns['metadata_sources'], list):
                adata.uns['metadata_sources'].append(metadata['metadata_source'])
            else:
                # If it's something else, create a new list with both elements
                adata.uns['metadata_sources'] = [adata.uns['metadata_sources'], metadata['metadata_source']]
        else:
            # If it doesn't exist yet, create it as a list with a single element
            adata.uns['metadata_sources'] = [metadata['metadata_source']]        
        
    
    logger.info("Dataset-specific metadata applied successfully")
    return adata