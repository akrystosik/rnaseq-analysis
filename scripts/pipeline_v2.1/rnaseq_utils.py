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
import logging

# Set up logging
logger = logging.getLogger('rnaseq_utils')

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
GENCODE_MAPPING_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/gencode_v24_complete_mapping.csv")
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

NIH_RACE_CATEGORIES_LOWER = sorted(list(set([
    "american indian or alaska native", "asian", "black or african american",
    "native hawaiian or other pacific islander", "white", "more than one race",
    "multiethnic", "unknown or not reported", "other"
])))

NIH_ETHNICITY_CATEGORIES_LOWER = sorted(list(set([
    "hispanic or latino", "not hispanic or latino", "unknown or not reported"
])))

# Categories for the 'self_reported_ethnicity' string label column
# This can also be defined here, or constructed where needed if it combines the above.
# For consistency, let's define a base set.
SRE_BASE_CATEGORIES_LOWER = sorted(list(set(
    NIH_RACE_CATEGORIES_LOWER +
    ["hispanic or latino"] # Add if "hispanic or latino" can be a standalone primary category
    # Add other combined terms if you decide on them as standard string labels,
    # e.g., "white hispanic or latino"
)))

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
        'SPECIES_TO_NCBI_TAXON': load_json_mapping('species_to_taxon.json', {}),
        'ETHNICITY_TO_HANCESTRO': load_json_mapping('ethnicity_to_hancestro.json', {}) 
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

def ensure_serializable(obj):
    """
    Recursively convert an object to HDF5-compatible native Python types.
    Handles nested lists, tuples, dicts, numpy arrays, scalars, pandas NA, etc.
    Converts problematic structures (like list of dicts) to JSON strings.
    """
    # --- Handle None first ---
    if obj is None:
        return None

    # --- Handle basic types ---
    if isinstance(obj, (str, int, float, bool)):
        return obj

    # --- Handle numpy/pandas scalars early ---
    # Check specific numpy types before generic hasattr('item')
    if isinstance(obj, np.integer): return int(obj)
    if isinstance(obj, np.floating):
        return None if np.isnan(obj) else float(obj) # Handle NaN directly
    if isinstance(obj, np.bool_): return bool(obj)
    if isinstance(obj, np.void): return None # Or appropriate representation

    # --- Handle numpy arrays specifically ---
    if isinstance(obj, np.ndarray):
        if obj.size == 0: return [] # Empty list for empty array
        # Attempt conversion to list of basic types first
        try:
            if obj.dtype == 'object':
                 # Recursively serialize elements if object type
                 return [ensure_serializable(x) for x in obj.tolist()]
            else:
                 # For numeric/bool arrays, handle internal NaNs then convert
                 if np.issubdtype(obj.dtype, np.floating):
                     # Replace NaN with None representation
                     return [None if np.isnan(x) else float(x) for x in obj.tolist()]
                 elif np.issubdtype(obj.dtype, np.integer):
                     return [int(x) for x in obj.tolist()]
                 elif np.issubdtype(obj.dtype, np.bool_):
                      return [bool(x) for x in obj.tolist()]
                 else:
                     # Other basic types (like fixed-length strings) might be listable directly
                     return obj.tolist()
        except Exception as e:
             logger.warning(f"Could not convert numpy array (dtype: {obj.dtype}) to list: {e}. Trying JSON string.")
             # Fallback to JSON string for complex arrays
             try:
                  # Need a custom encoder for numpy types within JSON
                  class NumpyEncoder(json.JSONEncoder):
                      def default(self, obj):
                          if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                                              np.int16, np.int32, np.int64, np.uint8,
                                              np.uint16, np.uint32, np.uint64)):
                              return int(obj)
                          elif isinstance(obj, (np.float_, np.float16, np.float32,
                                                np.float64)):
                              if np.isnan(obj): return None
                              if np.isinf(obj): return str(obj)
                              return float(obj)
                          elif isinstance(obj, (np.bool_)):
                              return bool(obj)
                          elif isinstance(obj, (np.void)):
                              return None
                          elif isinstance(obj, np.ndarray):
                              # Crucial: Recursively call ensure_serializable here too
                              return ensure_serializable(obj.tolist())
                          return super(NumpyEncoder, self).default(obj)
                  return json.dumps(obj, cls=NumpyEncoder)
             except Exception as json_e:
                  logger.error(f"Could not serialize numpy array to JSON: {json_e}. Falling back to str().")
                  return str(obj)

    # --- Handle simple Iterables (list, tuple, set) ---
    if isinstance(obj, (list, tuple, set)):
        temp_list = list(obj)
        # Check if it *contains* any dictionaries - more robust
        if any(isinstance(item, dict) for item in temp_list):
            logger.debug("List contains dictionaries, converting elements and dumping to JSON string.")
            try:
                # Ensure inner dicts are fully serializable first
                serializable_inner_list = [ensure_serializable(item) for item in temp_list]
                return json.dumps(serializable_inner_list, default=str) # Add default=str for safety
            except Exception as e:
                logger.warning(f"Could not serialize list containing dicts to JSON: {e}. Falling back to str().")
                return str(obj) # Fallback
        else:
            # Otherwise, just serialize elements recursively
            return [ensure_serializable(x) for x in temp_list]

    # --- Handle Dictionaries ---
    if isinstance(obj, dict):
        # Recursively serialize keys (as strings) and values in dicts
        serializable_dict = {}
        for k, v in obj.items():
            serializable_dict[str(k)] = ensure_serializable(v)
        return serializable_dict

    # --- Handle other specific types AFTER basic/numpy/list/dict ---
    if isinstance(obj, pd.Timestamp): return obj.isoformat()
    if isinstance(obj, pd.Series):
        return ensure_serializable(obj.tolist()) # Convert to list first
    if isinstance(obj, pd.Categorical):
        # Convert categories to string list
        return [str(x) for x in obj.tolist()]
    if isinstance(obj, pd.DataFrame):
         logger.warning(f"Attempting to serialize DataFrame in uns. Converting to JSON string.")
         try:
             # Convert to list of records (dicts), ensuring values are serializable
             dict_list = [ensure_serializable(row.to_dict()) for _, row in obj.iterrows()]
             return json.dumps(dict_list, default=str) # Use default=str for safety
         except Exception as e:
              logger.error(f"Failed to convert DataFrame to JSON string: {e}. Using str().")
              return str(obj)

    # --- Robust pd.isna check (LAST check before final string conversion) ---
    # This should only get scalars or types pd.isna understands if logic above is correct
    try:
        # Avoid calling pd.isna on types it explicitly errors on (like multi-element arrays)
        # The checks above should have handled numpy arrays already.
        if not isinstance(obj, (np.ndarray)) and pd.isna(obj):
            return None
    except ValueError as ve:
        # Catch "The truth value of an array is ambiguous" if it still occurs
        logger.warning(f"pd.isna failed for type {type(obj)}: {ve}. Treating as non-NA.")
        pass # Let it fall through to string conversion

    # --- Fallback for any remaining types ---
    logger.warning(f"Converting unrecognized type {type(obj)} to string for serialization.")
    return str(obj)


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

# In rnaseq_utils.py

def standardize_metadata(metadata_df, dataset_name, mappings=None):
    """
    Standardize metadata across datasets to ensure consistent fields and values.
    Maps fields to standard ontology terms where possible.
    
    Args:
        metadata_df: Metadata DataFrame
        dataset_name: Source dataset name
        mappings: Dictionary with mapping information
        
    Returns:
        Standardized metadata DataFrame with ontology mappings, or None if a critical error occurs.
    """
    try:
        if metadata_df is None:
            logger.error(f"standardize_metadata called with None input DataFrame for {dataset_name}. Returning None.")
            return None

        # Load mappings if not provided
        if mappings is None:
            logger.info("standardize_metadata: Mappings not provided, loading defaults.")
            mappings = load_mappings() # Assuming load_mappings() is defined in this file or imported
    
        # Make a copy to avoid modifying the original
        df = metadata_df.copy()
        
        # Get specific mappings - ensure these keys exist in `mappings` or handle missing keys
        tissue_to_uberon = mappings.get('TISSUE_TO_UBERON', {})
        assay_to_efo = mappings.get('ASSAY_TO_EFO', {})
        sex_standardization = mappings.get('SEX_STANDARDIZATION', {})
        species_to_taxon = mappings.get('SPECIES_TO_NCBI_TAXON', {})
        hancestro_map_loaded = mappings.get('ETHNICITY_TO_HANCESTRO', {}) # For HANCESTRO IDs

        # Define CORE_METADATA_FIELDS if not globally available in this context
        # These are fields that are fundamental and should exist.
        CORE_METADATA_FIELDS = [
            'sample_id', 'subject_id', 'sex', 'age', 'tissue', 
            'dataset', 'data_type', 'expression_unit', 'species'
            # Add race/ethnicity fields once they are consistently populated by upstream processes
        ]

        # Add dataset identifier
        df['dataset'] = dataset_name
        
        # Ensure all core fields exist (initialize with appropriate defaults if missing)
        for field in CORE_METADATA_FIELDS:
            if field not in df.columns:
                logger.debug(f"Standardize_metadata for {dataset_name}: Core field '{field}' missing. Adding with default.")
                if field in ['age', 'expression_unit', 'tissue']: # Fields that might be empty string
                    df[field] = "" 
                else: # Fields that usually have an 'unknown' category
                    df[field] = "unknown" 
        
        # Standardize species
        if 'species' in df.columns:
            df['species_ontology'] = df['species'].astype(str).str.lower().map(
                lambda x: species_to_taxon.get(x, species_to_taxon.get("unknown", ""))
            ).fillna("")
        else: # Should have been added by CORE_METADATA_FIELDS check
            df['species'] = "human" # Default
            df['species_ontology'] = species_to_taxon.get("human", "")

        # Standardize sex values
        if 'sex' in df.columns:
            # Ensure 'sex' column is string type before mapping for robustness
            df['sex'] = df['sex'].astype(str).str.lower().map(
                lambda x: sex_standardization.get(x, "unknown")
            ).fillna("unknown")
        else: # Should have been added by CORE_METADATA_FIELDS check
            df['sex'] = "unknown"


        # Standardize tissue to Uberon ontology
        if 'tissue' in df.columns:
            tissue_source_data = df['tissue']
            if isinstance(tissue_source_data, pd.DataFrame): # Defensive check
                logger.warning(
                    f"standardize_metadata for {dataset_name}: df['tissue'] is a DataFrame. Using its first column for 'tissue_original'. "
                    f"Columns found: {tissue_source_data.columns.tolist()}"
                )
                if not tissue_source_data.empty and len(tissue_source_data.columns) > 0:
                    df['tissue_original'] = tissue_source_data.iloc[:, 0].copy()
                    df['tissue'] = tissue_source_data.iloc[:, 0].copy() 
                else:
                    logger.error(f"standardize_metadata for {dataset_name}: df['tissue'] is an empty DataFrame or has no columns.")
                    df['tissue_original'] = "" 
                    df['tissue'] = ""          
            else:
                df['tissue_original'] = df['tissue'].copy() 

            if isinstance(df['tissue'], pd.DataFrame): # Final check
                logger.error(f"standardize_metadata for {dataset_name}: df['tissue'] is still a DataFrame. Aborting tissue processing.")
                df['tissue_ontology'] = "" 
            else:
                normalized_tissue = df['tissue'].astype(str).str.lower().str.strip()
                df['tissue'] = normalized_tissue # Update main 'tissue' column
                
                # Ensure tissue_to_uberon keys are lowercase for mapping
                tissue_to_uberon_lower = {k.lower(): v for k,v in tissue_to_uberon.items()}
                df['tissue_ontology'] = normalized_tissue.map(tissue_to_uberon_lower).fillna("")
                
                # Fallback substring matching for unmapped tissues (optional, can be noisy)
                # unmapped_mask = (df['tissue_ontology'] == "") & (df['tissue'] != "")
                # if unmapped_mask.any():
                #     # ... (substring matching logic, careful with this as it can mis-map) ...
                #     pass
        else: # if 'tissue' column doesn't exist, ensure ontology column does
             df['tissue_ontology'] = ""
             if 'tissue_original' not in df.columns: df['tissue_original'] = ""


        # Standardize data_type to EFO ontology
        if 'data_type' in df.columns:
            df['data_type_original'] = df['data_type'].copy()
            df['assay_ontology'] = df['data_type'].astype(str).str.lower().map(assay_to_efo).fillna("")
            
            if 'extraction_method' in df.columns:
                mask = (df['assay_ontology'] == "") & (df['extraction_method'].notna())
                df.loc[mask, 'assay_ontology'] = df.loc[mask, 'extraction_method'].astype(str).str.lower().map(assay_to_efo).fillna("")
        else: # Should have been added by CORE_METADATA_FIELDS
            df['data_type_original'] = ""
            df['assay_ontology'] = ""


        # Handle subject_id fallback
        if 'subject_id' not in df.columns or df['subject_id'].fillna('').eq('').all():
            if 'donor_id' in df.columns and not df['donor_id'].fillna('').eq('').all():
                logger.info(f"Using donor_id as subject_id for {dataset_name}")
                df['subject_id'] = df['donor_id']
            elif 'donor' in df.columns and not df['donor'].fillna('').eq('').all():
                logger.info(f"Using donor as subject_id for {dataset_name}")
                df['subject_id'] = df['donor']    
            else: # Ensure subject_id column exists even if empty
                df['subject_id'] = "unknown_subject_" + df.index.astype(str) # Fallback if truly missing
        
        # Log unmapped tissues (using map_tissue_to_ontology helper if it exists and is preferred)
        # For now, simple check based on tissue_ontology column
        if 'tissue' in df.columns and 'tissue_ontology' in df.columns:
            unmapped_tissues_final = df.loc[(df['tissue_ontology'] == "") & (df['tissue'] != ""), 'tissue'].dropna().unique()
            if len(unmapped_tissues_final) > 0:
                logger.warning(f"{dataset_name}: Unmapped tissues after processing: {', '.join(unmapped_tissues_final[:20])}")
                    

        # Standardize age and map to HsapDv
        if 'age' in df.columns:
            # --- BEGIN MODIFICATION ---
            age_source_data = df['age']
            if isinstance(age_source_data, pd.DataFrame): # Defensive check for age
                logger.warning(
                    f"standardize_metadata for {dataset_name}: df['age'] is a DataFrame. Using its first sub-column."
                )
                if not age_source_data.empty and len(age_source_data.columns) > 0:
                    # Ensure age_original captures the original problematic structure's first column if needed, or just the scalar
                    df['age_original'] = age_source_data.iloc[:, 0].copy().astype(str) # Keep original as series of strings
                    df['age'] = age_source_data.iloc[:, 0].copy().astype(str).replace('nan', '').replace('None', '').str.strip() # df['age'] becomes a Series
                else: # DataFrame column exists but is empty or has no sub-columns
                    logger.warning(f"standardize_metadata for {dataset_name}: df['age'] DataFrame is empty/malformed. Setting to empty string.")
                    df['age_original'] = ""
                    df['age'] = "" # df['age'] becomes a Series of empty strings
            else: # age column is already a Series
                df['age_original'] = df['age'].copy() # df['age'] is a Series
                df['age'] = df['age'].astype(str).replace('nan', '').replace('None', '').str.strip()
            # --- END MODIFICATION ---

            age_to_hsapdv_map = mappings.get('AGE_TO_HSAPDV', {})
            # df['age'] is now guaranteed to be a Series
            df['developmental_stage_ontology'] = df['age'].apply(
                lambda x: map_age_to_hsapdv(x, age_to_hsapdv_map) # map_age_to_hsapdv needs to be defined/imported
            ).fillna("")
        else: # age column doesn't exist
             logger.debug(f"Standardize_metadata: Adding missing 'age' column with default for {dataset_name}")
             df['age_original'] = ""
             df['age'] = "" # Ensure age column exists as string Series
             df['developmental_stage_ontology'] = ""            
            
        

        # --- Standardize Self-Reported Ethnicity and Race ---
        logger.debug(f"Standardizing ethnicity/race for {dataset_name}")


        # Ensure 'race' and 'is_hispanic_or_latino' columns exist, populated by dataset-specific logic or defaults here
        if 'race' not in df.columns: df['race'] = 'unknown or not reported'
        if 'is_hispanic_or_latino' not in df.columns: df['is_hispanic_or_latino'] = 'unknown or not reported'

        df['race_original'] = df['race'].astype(str).fillna('')
        df['race'] = df['race'].astype(str).str.lower().str.strip().fillna('unknown or not reported')
        df['race'] = df['race'].apply(lambda x: x if x in NIH_RACE_CATEGORIES_LOWER else 'unknown or not reported')
        
        df['ethnicity_original_is_hispanic'] = df['is_hispanic_or_latino'].astype(str).fillna('')
        df['is_hispanic_or_latino'] = df['is_hispanic_or_latino'].astype(str).str.lower().str.strip().fillna('unknown or not reported')
        df['is_hispanic_or_latino'] = df['is_hispanic_or_latino'].apply(lambda x: x if x in NIH_ETHNICITY_CATEGORIES_LOWER else 'unknown or not reported')

        sre_values = []
        for index, row in df.iterrows():
            race_val = str(row.get('race', 'unknown or not reported')) # Use .get for safety
            is_hispanic_val = str(row.get('is_hispanic_or_latino', 'unknown or not reported'))

            if race_val == "more than one race":
                sre_values.append("multiethnic")
            elif race_val == "unknown or not reported" and is_hispanic_val == "unknown or not reported":
                sre_values.append("unknown or not reported")
            elif is_hispanic_val == "hispanic or latino":
                if race_val != "unknown or not reported" and race_val != "multiethnic":
                    # Example: if race is "white", self_reported_ethnicity could be "white hispanic or latino"
                    # Or, if your HANCESTRO map has a term for "White" and another for "Hispanic or Latino",
                    # this combination might be mapped to "multiethnic" for the ontology ID.
                    # For now, let's create a combined label or default to a broader term.
                    # This decision depends on your HANCESTRO mapping strategy.
                    # A simple approach for now:
                    sre_values.append(f"{race_val} (hispanic or latino)") # Or simply "hispanic or latino" or race_val or "multiethnic"
                else: # race is unknown or multiethnic, but is Hispanic
                    sre_values.append("hispanic or latino")
            else: # Not Hispanic or Latino, or ethnicity unknown. Default to race value.
                sre_values.append(race_val)
        
        # Define categories for self_reported_ethnicity based on potential values
        sre_categories = sorted(list(set(NIH_RACE_CATEGORIES_LOWER + \
                                         [f"{r} (hispanic or latino)" for r in NIH_RACE_CATEGORIES_LOWER if r not in ["unknown or not reported", "multiethnic"]] + \
                                         ["hispanic or latino"])))

        df['self_reported_ethnicity'] = pd.Categorical(
            sre_values,
            categories=sre_categories,
            ordered=False
        ).fillna("unknown or not reported")

        hancestro_map_from_json = mappings.get('ETHNICITY_TO_HANCESTRO', {}) # Renamed for clarity
        hancestro_map_lower = {k.lower(): v for k, v in hancestro_map_from_json.items()}
        
        df['self_reported_ethnicity_ontology_term_id'] = df['self_reported_ethnicity'].astype(str).str.lower().map(hancestro_map_lower).fillna('unknown')

        # Validate against schema-allowed values for the ontology ID
        # The values in hancestro_map_lower.values() are the HANCESTRO IDs
        temp_allowed_ids = list(hancestro_map_lower.values())
        temp_allowed_ids.append("multiethnic")
        temp_allowed_ids.append("unknown")
        schema_allowed_ontology_ids = sorted(list(set(temp_allowed_ids))) # <<< ENSURE UNIQUE AND SORTED
        
        invalid_ontology_mask = ~df['self_reported_ethnicity_ontology_term_id'].isin(schema_allowed_ontology_ids)
        if invalid_ontology_mask.any():
            # ... (warning log) ...
            df.loc[invalid_ontology_mask, 'self_reported_ethnicity_ontology_term_id'] = 'unknown'
        
        # Ensure the final column is categorical with these allowed_ontology_ids
        df['self_reported_ethnicity_ontology_term_id'] = pd.Categorical(
            df['self_reported_ethnicity_ontology_term_id'],
            categories=schema_allowed_ontology_ids, # Use the unique, sorted list
            ordered=False
        )
        
        # --- Finalize Categorical Dtypes for all relevant fields ---
        # Ensure all fields intended to be categorical are set as such, with consistent categories where possible.
        # This list should include all fields that are meant to be categorical in the final obs.

        final_categorical_fields = [
            'sex', 'tissue', 'dataset', 'species', 'data_type', 
            'race', 'is_hispanic_or_latino', 'self_reported_ethnicity', 
            'self_reported_ethnicity_ontology_term_id',
            'tissue_ontology', 'assay_ontology', 'species_ontology', 'developmental_stage_ontology'
        ]
        if 'cell_type' in df.columns: final_categorical_fields.append('cell_type')
        if 'disease' in df.columns: final_categorical_fields.append('disease')

        for field in final_categorical_fields:
            if field in df.columns:
                current_categories_for_field = None # Use specific lists if available
                default_val_for_cat = 'unknown'
                
                if field == 'race': 
                    current_categories_for_field = NIH_RACE_CATEGORIES_LOWER # Uses module global
                    default_val_for_cat = 'unknown or not reported'
                elif field == 'is_hispanic_or_latino': 
                    current_categories_for_field = NIH_ETHNICITY_CATEGORIES_LOWER # Uses module global
                    default_val_for_cat = 'unknown or not reported'
                elif field == 'self_reported_ethnicity':
                    current_categories_for_field = SRE_BASE_CATEGORIES_LOWER # Uses module global
                    default_val_for_cat = 'unknown or not reported'
                elif field == 'self_reported_ethnicity_ontology_term_id':
                    # This list is dynamically built based on hancestro_map earlier in this function
                    current_categories_for_field = schema_allowed_ontology_ids 
                
                try:
                    # Ensure all values are strings before converting to category, handle NaNs
                    str_series = df[field].astype(str).replace('nan', default_val_for_cat).fillna(default_val_for_cat)
                    if current_categories_for_field:
                        # Ensure all values in series are part of defined categories
                        # And that current_categories_for_field itself has unique values
                        unique_cats = sorted(list(set(current_categories_for_field)))
                        valid_series = str_series.apply(lambda x: x if x in unique_cats else default_val_for_cat)
                        df[field] = pd.Categorical(valid_series, categories=unique_cats, ordered=False)
                    else: # Let pandas infer categories for other fields like tissue_ontology, etc.
                        df[field] = str_series.astype('category')
                except Exception as e_cat_final:
                    logger.error(f"Error converting field '{field}' to category for {dataset_name}: {e_cat_final}. Leaving as string.")
                    df[field] = df[field].astype(str).fillna(default_val_for_cat)
            else:
                logger.debug(f"Standardize_metadata: Adding missing categorical field '{field}' as 'unknown' for {dataset_name}.")
                df[field] = pd.Categorical(['unknown'] * len(df), categories=['unknown'])
        return df
    except Exception as e_std_meta:
        logger.error(f"CRITICAL ERROR within standardize_metadata for {dataset_name}: {e_std_meta}")
        import traceback
        logger.error(traceback.format_exc())
        return None # Ensure None is returned on critical failure

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
            # Store directly in uns - ensure_serializable will handle it later
            adata.uns[key] = value
            # Note: Complex nested structures are handled by ensure_serializable before saving

    # Apply observation column updates if present
    if 'obs_columns' in metadata:
        for col_name, value in metadata['obs_columns'].items():
            logger.info(f"  Updating obs column {col_name}")
            # Check if value should be applied to all rows or is a map
            if isinstance(value, dict):
                 logger.warning(f"Applying dict values to obs column '{col_name}' is not directly supported here. Check metadata structure.")
                 # Potentially map using index if value is a dict keyed by sample ID
            elif isinstance(value, list) and len(value) == adata.n_obs:
                 adata.obs[col_name] = value # Assign list if length matches
            elif not isinstance(value, (dict, list)):
                 adata.obs[col_name] = value # Assign scalar value to all rows
            else:
                 logger.warning(f"Could not apply value of type {type(value)} to obs column '{col_name}'")


    # Add metadata source information
    if 'metadata_source' in metadata:
        if 'metadata_sources' not in adata.uns:
            adata.uns['metadata_sources'] = []
        elif not isinstance(adata.uns['metadata_sources'], list):
            adata.uns['metadata_sources'] = [str(adata.uns['metadata_sources'])]

        source_info = metadata['metadata_source'] # Get the actual value
        
        # Ensure all elements in source_info, if it's a list, are serializable
        # This is more for the 'extras' field which is a list of dicts.
        # For metadata_source, it's usually a single dict.
        
        if isinstance(source_info, dict):
            try:
                # Serialize the source_info dict itself
                serializable_source_info = ensure_serializable(source_info) # Use ensure_serializable from standardize_datasets.py
                source_str = json.dumps(serializable_source_info, default=str)
                adata.uns['metadata_sources'].append(source_str)
                logger.debug(f"Appended metadata_source as JSON string: {source_str[:100]}...")
            except TypeError as te: # Handles cases where ensure_serializable might not fully clean it for json.dumps
                logger.warning(f"Could not directly JSON dump complex metadata_source: {te}. Converting to string.")
                adata.uns['metadata_sources'].append(str(source_info))
            except Exception as e:
                logger.error(f"Could not serialize metadata_source dict to JSON: {e}. Appending as raw string.")
                adata.uns['metadata_sources'].append(str(source_info))
        elif isinstance(source_info, str): # If it's already a string
            adata.uns['metadata_sources'].append(source_info)
        else: # Other types, convert to string
            adata.uns['metadata_sources'].append(str(source_info))    
    

    logger.info("Dataset-specific metadata applied successfully")
    return adata
