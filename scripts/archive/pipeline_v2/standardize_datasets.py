#!/usr/bin/env python3
"""
Multi-Dataset Standardization Pipeline

This script processes RNA-seq and expression data from multiple sources 
(ENCODE, GTEx, MAGE, ADNI, ENTEx) into standardized AnnData objects.
It uses JSON configuration files for metadata and ontology mappings.

Usage:
  python standardize_datasets.py --encode-dir /path/to/encode/data \
                                 --gtex-file /path/to/gtex.gct.gz \
                                 --mage-dir /path/to/mage/data \
                                 --adni-dir /path/to/adni/data \
                                 --metadata-dir /path/to/metadata/json \
                                 --output-dir /path/to/output

  # For ENTEx processing:
  python standardize_datasets.py --metadata /path/to/entex_metadata.json \
                                 --dataset-type entex \
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
import scipy.sparse as sp 

# Import shared utilities
from rnaseq_utils import (
    standardize_ensembl_id,
    load_gencode_mapping,
    add_gencode_annotations,
    standardize_metadata,
    validate_metadata,
    prepare_for_anndata,
    load_mappings,
    load_json_mapping,
    load_dataset_specific_metadata,
    apply_dataset_specific_metadata,
)


# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("data_standardizer")

# Define paths and constants
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
DEFAULT_OUTPUT_DIR = BASE_DIR / "standardized_data"
DEFAULT_METADATA_DIR = BASE_DIR / "metadata/json"

# GTEx metadata files
GTEx_SAMPLE_ATTRS = Path(
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
)
GTEx_SUBJECT_ATTRS = Path(
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
)


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
    if isinstance(obj, (np.integer)): return int(obj)
    if isinstance(obj, (np.floating)):
        return None if np.isnan(obj) else float(obj) # Handle NaN directly
    if isinstance(obj, (np.bool_)): return bool(obj)

    # --- Handle numpy arrays specifically ---
    if isinstance(obj, np.ndarray):
        if obj.size == 0: return [] # Empty list for empty array
        if obj.dtype == 'object':
             # Recurse on elements if object type
             return [ensure_serializable(x) for x in obj.tolist()]
        else:
             # Handle NaNs within numeric arrays before tolist
             if np.issubdtype(obj.dtype, np.floating):
                 # Replace NaN with None representation
                 return [None if np.isnan(x) else float(x) for x in obj.tolist()]
             else:
                 # Other basic types can be listed
                 return obj.tolist()

    # --- Handle lists/tuples/sets ---
    if isinstance(obj, (list, tuple, set)):
         # Check if it's a list of dictionaries (problematic for h5py arrays)
        is_list_of_dicts = False
        temp_list = list(obj) # Convert set/tuple to list for check
        if temp_list and all(isinstance(item, dict) for item in temp_list):
            is_list_of_dicts = True

        if is_list_of_dicts:
            # Convert list of dicts to JSON string
            logger.debug("Converting list of dicts to JSON string.")
            try:
                # Ensure dicts within the list are also serializable before dumping
                serializable_inner_list = [ensure_serializable(item) for item in temp_list]
                return json.dumps(serializable_inner_list)
            except Exception as e:
                logger.warning(f"Could not serialize list of dicts to JSON: {e}. Falling back to str().")
                return str(obj)
        else:
            # Recursively serialize elements in other lists/tuples/sets
            return [ensure_serializable(x) for x in temp_list] # Use temp_list here


    # --- Handle Dictionaries ---
    if isinstance(obj, dict):
        # Recursively serialize keys (as strings) and values in dicts
        # Check if any value remains complex after serialization attempt
        serialized_dict = {}
        contains_complex = False
        for k, v in obj.items():
            serialized_val = ensure_serializable(v)
            serialized_dict[str(k)] = serialized_val
            if isinstance(serialized_val, (dict, list)): # If nesting remains
                 contains_complex = True

        # If dict still seems complex, consider converting the whole dict to JSON string
        # This is aggressive but might catch nested issues h5py dislikes
        # if contains_complex:
        #      logger.debug(f"Dictionary might still contain complex nested types. Converting to JSON string.")
        #      try:
        #           return json.dumps(serialized_dict)
        #      except Exception as e:
        #           logger.warning(f"Could not serialize complex dict to JSON: {e}. Falling back to str().")
        #           return str(obj)
        # else:
        return serialized_dict # Return the dict if values seem simple

    # --- Handle other specific types ---
    if isinstance(obj, pd.Series):
        return ensure_serializable(obj.tolist()) # Convert to list first
    if isinstance(obj, pd.Timestamp):
        return obj.isoformat()
    if isinstance(obj, pd.Categorical):
        return [str(x) for x in obj.tolist()]
    if isinstance(obj, pd.DataFrame):
         logger.warning(f"Attempting to serialize DataFrame in uns. Converting to JSON string.")
         try:
             dict_list = obj.astype(str).to_dict(orient='records')
             return json.dumps(dict_list) # Convert list of dicts to JSON
         except Exception as e:
              logger.error(f"Failed to convert DataFrame to JSON string: {e}. Using str().")
              return str(obj)

    # --- Robust pd.isna check (AFTER specific type checks) ---
    # This should now only get scalars or types pd.isna understands
    try:
        if pd.isna(obj):
            return None
    except ValueError: # Catch the "empty array" error if it somehow still gets here
        logger.warning(f"pd.isna failed for type {type(obj)}, likely empty array slipped through. Treating as None/empty.")
        return None


    # --- Fallback for any other type ---
    logger.warning(f"Converting unrecognized type {type(obj)} to string.")
    return str(obj)



def identify_columns(df):
    """Identify gene_id and TPM columns across different file formats"""

    # Initialize column mappings
    gene_id_col = None
    tpm_col = None

    # Check column names (case-insensitive)
    columns_lower = [col.lower() for col in df.columns]

    # Check for gene identifier column
    if "gene_id" in df.columns:
        gene_id_col = "gene_id"
    elif "target_id" in df.columns:
        gene_id_col = "target_id"

    # Check for TPM column (case-insensitive)
    if "tpm" in columns_lower:
        tpm_col = df.columns[columns_lower.index("tpm")]

    return gene_id_col, tpm_col


def get_encode_cell_info(metadata_dir=None):
    """Get ENCODE cell line info from JSON."""
    if metadata_dir:
        json_path = os.path.join(metadata_dir, "encode_metadata.json")
        if os.path.exists(json_path):
            try:
                with open(json_path, "r") as f:
                    metadata = json.load(f)
                    # Extract cell line info from the JSON
                    if "dataset_info" in metadata and "cell_lines" in metadata["dataset_info"]:
                        return metadata["dataset_info"]["cell_lines"]
            except Exception as e:
                logger.warning(f"Error loading ENCODE cell info: {e}")

    # Provide a minimal fallback dictionary instead of referring to the removed ENCODE_CELL_INFO
    logger.warning("Could not load cell line info from JSON, using minimal fallback values")
    return {
        "A549": {"tissue": "lung", "organism": "human"},
        "K562": {"tissue": "bone marrow", "organism": "human"},
        "HepG2": {"tissue": "liver", "organism": "human"},
        # Add minimal entries for other cell lines if needed
    }


def create_standard_anndata(data_df, obs_df, var_df, dataset_info):
    # Ensure all obs columns can be properly saved
    for col in obs_df.columns:
        if obs_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting obs column {col} to string")
            obs_df[col] = obs_df[col].astype(str)

    # Ensure all var columns can be properly saved
    for col in var_df.columns:
        if var_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting var column {col} to string")
            var_df[col] = var_df[col].astype(str)

    # Convert all dictionary values to strings
    safe_dataset_info = {}
    for key, value in dataset_info.items():
        if isinstance(value, dict):
            safe_dataset_info[key] = {k: str(v) if v is not None else "" for k, v in value.items()}
        elif value is None:
            safe_dataset_info[key] = ""
        else:
            safe_dataset_info[key] = str(value)
    """
    Create standardized AnnData object with consistent structure.
    Now with enhanced metadata including ontology mappings and validation.

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
    # Ensure all obs columns can be properly saved
    for col in obs_df.columns:
        if obs_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting obs column {col} to string")
            obs_df[col] = obs_df[col].astype(str)

    # Ensure all var columns can be properly saved
    for col in var_df.columns:
        if var_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting var column {col} to string")
            var_df[col] = var_df[col].astype(str)

    # Create AnnData object
    adata = ad.AnnData(X=data_df.values, obs=obs_df, var=var_df)

    # Ensure all uns values are serializable
    for key in dataset_info.keys():
        if isinstance(dataset_info[key], (pd.Series, pd.DataFrame, np.ndarray)):
            logger.debug(
                f"Converting uns[{key}] from {type(dataset_info[key])} to standard Python type"
            )
            dataset_info[key] = (
                dataset_info[key].tolist()
                if hasattr(dataset_info[key], "tolist")
                else str(dataset_info[key])
            )

    # Validate metadata before preparing for AnnData
    dataset_name = dataset_info.get("source", "unknown")
    validated_obs_df, validation_report = validate_metadata(obs_df, dataset_name)

    # Prepare data for AnnData
    validated_obs_df, var_df, data_df = prepare_for_anndata(validated_obs_df, var_df, data_df)

    # Create AnnData object (assuming data_df is samples x genes)
    logger.debug(
        f"Creating AnnData with X shape={data_df.values.shape}, obs shape={validated_obs_df.shape}, var shape={var_df.shape}"
    )
    adata = ad.AnnData(X=data_df.values, obs=validated_obs_df, var=var_df)

    # Verify metadata was transferred correctly
    if adata.var.shape[1] == 0:
        logger.error("ERROR: var DataFrame is empty after AnnData creation!")
        logger.error(
            f"Original var_df had shape {var_df.shape} with columns: {list(var_df.columns)}"
        )

        # Emergency fix - manually add the metadata if it was lost
        for col in var_df.columns:
            adata.var[col] = var_df[col].values

    if adata.obs.shape[1] == 0:
        logger.error("ERROR: obs DataFrame is empty after AnnData creation!")
        logger.error(
            f"Original obs_df had shape {validated_obs_df.shape} with columns: {list(validated_obs_df.columns)}"
        )

        # Emergency fix - manually add the metadata if it was lost
        for col in validated_obs_df.columns:
            adata.obs[col] = validated_obs_df[col].values

    # Add dataset information to uns
    # Ensure dataset_info is serializable
    dataset_info = ensure_serializable(dataset_info)
    adata.uns["dataset_info"] = dataset_info

    # Add reference genome info to uns
    adata.uns["reference_genome"] = "hg38"
    adata.uns["gencode_version"] = dataset_info.get("gencode_version", 24)

    # Add validation report to uns
    adata.uns["metadata_validation"] = validation_report

    # Add ontology mapping dictionaries to uns for reference
    mappings = load_mappings()
    adata.uns["ontology_mappings"] = {
        "tissue_to_uberon": mappings.get("TISSUE_TO_UBERON", {}),
        "assay_to_efo": mappings.get("ASSAY_TO_EFO", {}),
        "age_to_hsapdv": mappings.get("AGE_TO_HSAPDV", {}),
        "species_to_ncbi_taxon": mappings.get("SPECIES_TO_NCBI_TAXON", {}),
    }

    # Add processing timestamp
    adata.uns["processing_date"] = pd.Timestamp.now().strftime("%Y-%m-%d")

    return adata



def save_anndata(adata, file_path):
    """Save AnnData object to file with validation and improved serialization."""
    try:
        logger.info(f"Preparing to save AnnData to {file_path}")

        # --- Pre-Save Validation and Cleanup ---
        if adata.n_vars == 0:
            logger.error("Cannot save AnnData: var DataFrame is empty!")
            return False
        if 'gene_id' not in adata.var.columns and len(adata.var_names) > 0:
             logger.warning("Reconstructing potentially missing 'gene_id' column in var from index.")
             adata.var['gene_id'] = adata.var_names
        elif adata.var.shape[1] == 0 and len(adata.var_names) > 0:
             logger.error("Cannot save AnnData: var DataFrame has no columns! Attempting reconstruction.")
             adata.var['gene_id'] = adata.var_names


        # Ensure index names don't conflict with columns
        if adata.obs.index.name is not None and adata.obs.index.name in adata.obs.columns:
            new_name = f"original_obs_{adata.obs.index.name}"
            logger.warning(f"Fixing index/column name conflict in obs: Renaming column '{adata.obs.index.name}' to '{new_name}'")
            adata.obs = adata.obs.rename(columns={adata.obs.index.name: new_name})

        if adata.var.index.name is not None and adata.var.index.name in adata.var.columns:
            new_name = f"original_var_{adata.var.index.name}"
            logger.warning(f"Fixing index/column name conflict in var: Renaming column '{adata.var.index.name}' to '{new_name}'")
            adata.var = adata.var.rename(columns={adata.var.index.name: new_name})

        # Ensure obs/var indices are strings
        adata.obs.index = adata.obs.index.astype(str)
        adata.var.index = adata.var.index.astype(str)

        # --- Explicitly convert object columns in obs/var to string ---
        # This addresses errors like the one with 'cell_type_info' in obs
        logger.info("Converting object columns in obs/var to string before saving...")
        for df_name, df in [('obs', adata.obs), ('var', adata.var)]:
            for col in df.columns:
                if df[col].dtype == 'object':
                    try:
                        # Attempt conversion, filling potential NaNs introduced during conversion
                        # This handles cases where objects might be complex (like dicts or lists)
                        df[col] = df[col].astype(str).fillna('')
                        logger.debug(f"Converted {df_name} column '{col}' to string.")
                    except Exception as e:
                        logger.error(f"Failed to convert {df_name} object column '{col}' to string: {e}. Skipping column for safety.")
                        # Optionally drop the column if conversion fails catastrophically
                        # df = df.drop(columns=[col])

        # Reassign potentially modified dataframes back to adata
        adata.obs = df if df_name == 'obs' else adata.obs
        adata.var = df if df_name == 'var' else adata.var
        logger.info("Object column conversion complete.")
        # --- End object column conversion ---

        # --- Serialize .uns ---
        logger.info("Ensuring adata.uns is serializable...")
        original_uns = adata.uns.copy() # Keep a copy of the original uns
        try:
            serializable_uns = ensure_serializable(original_uns)
            adata.uns = serializable_uns
            logger.info("adata.uns converted to serializable format.")
        except Exception as e:
             logger.error(f"Failed during uns serialization: {e}", exc_info=True)
             logger.warning("Attempting to save without full uns serialization.")
             adata.uns = original_uns # Restore original on failure
        # --- End .uns serialization ---


        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        # --- Save the AnnData object ---
        logger.info(f"Writing AnnData to {file_path}...")
        try:
            adata.write_h5ad(file_path)
            logger.info(f"Successfully wrote AnnData to {file_path}")
            save_successful = True
        except Exception as write_e:
             logger.error(f"Error writing AnnData to {file_path}: {write_e}", exc_info=True)
             save_successful = False
        # --- End Save ---

        # Restore original uns in the in-memory object after saving attempt
        adata.uns = original_uns

        # --- Post-Save Verification ---
        if save_successful:
            logger.info(f"Verifying saved file: {file_path}")
            try:
                test_load = ad.read_h5ad(file_path)
                logger.info(
                    f"Verification successful: Loaded AnnData with shape={test_load.shape}, "
                    f"var columns={list(test_load.var.columns)}, obs columns={len(test_load.obs.columns)}"
                )
                if test_load.var.shape[1] == 0:
                     logger.error("Verification FAILED: Loaded AnnData has empty var DataFrame!")
                     return False
                if test_load.obs.shape[1] == 0:
                     logger.error("Verification FAILED: Loaded AnnData has empty obs DataFrame!")
                     return False
                return True # Return True on successful save and verification

            except Exception as load_e:
                logger.error(f"Verification FAILED: Error loading saved file {file_path}: {load_e}")
                return False # Saving failed if we can't reload it
        else:
            return False # Save failed

    except Exception as e:
        logger.error(f"Unhandled error in save_anndata for {file_path}: {e}", exc_info=True)
        # Attempt to restore original uns even on error
        if 'original_uns' in locals() and 'adata' in locals():
             adata.uns = original_uns
        return False # Return False on error

def export_metadata_to_json(metadata_df, output_file):
    """
    Export standardized metadata to a JSON file.

    Parameters:
    -----------
    metadata_df : pandas.DataFrame
        Standardized metadata DataFrame
    output_file : str
        Path to save the JSON file

    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Convert DataFrame to dictionary
        metadata_dict = {}

        # Create a record for each sample
        for idx, row in metadata_df.iterrows():
            sample_id = idx if not pd.isna(idx) else f"sample_{len(metadata_dict)}"
            metadata_dict[sample_id] = row.to_dict()

            # Convert numpy and pandas types to Python native types
            for key, value in metadata_dict[sample_id].items():
                if isinstance(value, (np.integer, np.floating)):
                    metadata_dict[sample_id][key] = value.item()
                elif pd.isna(value):
                    metadata_dict[sample_id][key] = None

        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Write to JSON file
        with open(output_file, "w") as f:
            json.dump(metadata_dict, f, indent=2)

        logger.info(f"Exported metadata to {output_file}")
        return True

    except Exception as e:
        logger.error(f"Error exporting metadata: {e}")
        import traceback

        logger.error(traceback.format_exc())
        return False


# ====================================================================================
# ENCODE-specific Processing
# ====================================================================================
def load_entex_metadata(metadata_file=None):
    """
    Load ENTEx metadata from JSON file.

    Parameters:
    -----------
    metadata_file : str or Path, optional
        Path to the ENTEx metadata JSON file

    Returns:
    --------
    dict
        Dictionary containing ENTEx sample metadata and donor information
    """
    try:
        if metadata_file and os.path.exists(metadata_file):
            metadata_path = metadata_file
        else:
            # Use default path if not specified
            metadata_path = Path(
                "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/metadata/entex_metadata.json"
            )

        if not os.path.exists(metadata_path):
            logger.warning(f"ENTEx metadata file not found: {metadata_path}")
            return {"donor_map": {}, "entex_metadata": [], "sample_lookup": {}}

        logger.info(f"Loading ENTEx metadata from {metadata_path}")
        with open(metadata_path, "r") as f:
            metadata = json.load(f)

        # Initialize sample_lookup dictionary
        sample_lookup = {}
        donor_map = metadata.get("donor_map", {})

        # First, load samples from the entex_metadata array if it exists
        if "entex_metadata" in metadata and isinstance(metadata["entex_metadata"], list):
            for sample in metadata["entex_metadata"]:
                sample_id = sample.get("sample_id")
                if sample_id:
                    sample_lookup[sample_id] = sample
            logger.info(f"Loaded {len(sample_lookup)} samples from entex_metadata array")

        # Next, check for individual sample entries at the root level (like ENCFF* keys)
        for key, value in metadata.items():
            if key.startswith("ENCFF") and isinstance(value, dict):
                # This looks like a sample entry directly in the root
                sample_id = key
                sample_lookup[sample_id] = value

                # Ensure it has sample_id field for consistency
                if "sample_id" not in value:
                    sample_lookup[sample_id]["sample_id"] = sample_id

        # Also incorporate any existing sample_lookup dictionary
        if "sample_lookup" in metadata and isinstance(metadata["sample_lookup"], dict):
            for sample_id, sample_info in metadata["sample_lookup"].items():
                sample_lookup[sample_id] = sample_info

        # Update the metadata with the complete sample_lookup
        metadata["sample_lookup"] = sample_lookup

        # If there's no donor_map, try to build it from samples
        if not donor_map:
            donor_map = {}
            for sample_id, info in sample_lookup.items():
                donor_id = info.get("donor", info.get("subject_id", info.get("donor_id")))
                if donor_id and donor_id not in donor_map:
                    donor_map[donor_id] = {"sex": info.get("sex", ""), "age": info.get("age", "")}
            metadata["donor_map"] = donor_map

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
        df = pd.read_csv(file_path, sep="\t", low_memory=False)

        # Use the identify_columns function to identify gene_id and TPM columns
        gene_id_col, tpm_col = identify_columns(df)

        if gene_id_col is not None and tpm_col is not None:
            # Set gene_id as index and handle numeric IDs
            df[gene_id_col] = df[gene_id_col].astype(str)
            # Extract only the columns we need and set index
            result_df = df[[gene_id_col, tpm_col]].set_index(gene_id_col)
            # Rename column to make it consistent
            result_df = result_df.rename(columns={tpm_col: "TPM"})
            return result_df
        else:
            logger.warning(f"Could not identify gene_id and TPM columns in {file_path}")
            logger.warning(f"Available columns: {list(df.columns)}")
            return None

    except Exception as e:
        logger.warning(f"Error reading TPM file {file_path}: {e}")
        return None


def extract_encode_metadata(file_path, entex_metadata=None, metadata_dir=None):
    """
    Extract metadata from ENCODE file path, now using JSON configurations.

    Parameters:
    -----------
    file_path : str
        Path to the ENCODE file
    entex_metadata : dict, optional
        ENTEx metadata with sample_lookup dictionary
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    dict
        Dictionary of metadata
    """
    metadata = {}

    # Get ENCODE cell line info
    encode_cell_info = get_encode_cell_info(metadata_dir)

    # Get file basename and extract file ID
    file_name = os.path.basename(file_path)
    file_id = file_name.split(".")[0]  # e.g., ENCFF008RKC

    # Store original sample ID
    metadata["original_sample_id"] = file_id

    # Check if this is an ENTEx file with available metadata
    if entex_metadata and "sample_lookup" in entex_metadata and "entex" in file_path.lower():
        if file_id in entex_metadata["sample_lookup"]:
            sample_info = entex_metadata["sample_lookup"][file_id]
            donor_id = sample_info.get(
                "donor_id", sample_info.get("donor", sample_info.get("subject_id"))
            )

            # Add basic metadata from lookup
            metadata["tissue"] = sample_info.get("tissue", "")
            metadata["subject_id"] = donor_id
            metadata["assay_type"] = sample_info.get("assay_type", "total")
            metadata["genome_annotation"] = sample_info.get("genome_annotation", "")
            metadata["genome_assembly"] = sample_info.get("genome_assembly", "")
            metadata["data_source"] = "ENTEx"

            # Add donor information if available
            if donor_id and donor_id in entex_metadata.get("donor_map", {}):
                donor_info = entex_metadata["donor_map"][donor_id]
                metadata["sex"] = donor_info.get("sex", "")
                metadata["age"] = donor_info.get("age", "")

            logger.debug(f"Found ENTEx metadata for sample {file_id}")

        else:
            # Extract minimal metadata from the directory structure
            logger.warning(
                f"No metadata found for ENTEx sample {file_id}, extracting from directory"
            )

            # Extract tissue from directory name (e.g., "ENTEx_adrenal gland_unknown")
            path_parts = file_path.split(os.sep)
            entex_parts = [part for part in path_parts if part.startswith("ENTEx_")]

            if entex_parts:
                # Parse tissue from directory name
                tissue_parts = entex_parts[0].split("_")
                if len(tissue_parts) > 1:
                    metadata["tissue"] = tissue_parts[1]

                # Add ENTEx flag
                metadata["data_source"] = "ENTEx"

            # Try to extract donor ID from path
            donor_parts = [part for part in path_parts if part.startswith("entex_tc_")]
            if donor_parts:
                metadata["subject_id"] = donor_parts[0].replace("entex_tc_", "")

    else:
        # Regular ENCODE cell line processing
        # Extract cell line from path - improved detection
        path_parts = file_path.split(os.sep)
        # Check each cell line more thoroughly
        cell_line_dir = []
        for part in path_parts:
            # Check exact match
            if part in encode_cell_info:
                cell_line_dir.append(part)
                continue
            
            # Check for cell line as a substring
            for cl in encode_cell_info.keys():
                if cl in part:
                    cell_line_dir.append(part)
                    break

        if cell_line_dir:
            # Handle directories like 'A549_polyA_plus'
            cell_line_full = cell_line_dir[0]

            # Extract cell line name
            cell_line = None
            for cl in encode_cell_info.keys():
                if cell_line_full.startswith(cl):
                    cell_line = cl
                    break

            if cell_line:
                metadata["cell_line"] = cell_line

                # Add basic cell line information
                if cell_line in encode_cell_info:
                    for key, value in encode_cell_info[cell_line].items():
                        metadata[key] = value

                # Extract RNA extraction method if in directory name
                if "_polyA_plus" in cell_line_full:
                    metadata["extraction_method"] = "polyA_plus"
                elif "_polyA_minus" in cell_line_full:
                    metadata["extraction_method"] = "polyA_minus"
                elif "_total" in cell_line_full:
                    metadata["extraction_method"] = "total"

    # Add standard fields
    metadata["data_type"] = "RNA-seq"
    metadata["expression_unit"] = "TPM"
    
    # Ensure tissue is never None/nan for cell lines
    if "cell_line" in metadata and not metadata.get("tissue"):
        cell_line = metadata["cell_line"]
        if cell_line in encode_cell_info and "tissue" in encode_cell_info[cell_line]:
            metadata["tissue"] = encode_cell_info[cell_line]["tissue"]
        else:
            # Fallback mapping for common cell lines
            fallback_mapping = {
                "A549": "lung",
                "HepG2": "liver", 
                "K562": "blood",
                "HeLa": "cervix",
                "IMR-90": "lung",
                "GM12878": "blood",
                "H1-hESC": "embryonic stem cell"
            }
            if cell_line in fallback_mapping:
                metadata["tissue"] = fallback_mapping[cell_line]

    return metadata


def process_encode_data(
    input_dir, entex_dir=None, output_file=None, entex_metadata_file=None, metadata_dir=None
):
    """
    Process ENCODE RNA-seq data into standardized AnnData, using JSON configurations.

    Parameters:
    -----------
    input_dir : str
        Directory containing ENCODE cell line TPM files
    entex_dir : str, optional
        Directory containing ENCODE ENTEx TPM files
    output_file : str
        Path to save the standardized AnnData
    entex_metadata_file : str, optional
        Path to ENTEx metadata JSON file
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    logger.info(f"Processing ENCODE data from {input_dir}")

    # Load ENTEx metadata
    entex_metadata = (
        load_entex_metadata(entex_metadata_file) if (entex_dir or entex_metadata_file) else None
    )

    # Find all TSV files in cell lines
    tsv_files = []
    if input_dir and os.path.exists(input_dir):
        tsv_files = glob.glob(os.path.join(input_dir, "**", "*.tsv"), recursive=True)

    # Also look for ENTEx files if directory is provided
    if entex_dir:
        logger.info(f"Checking ENTEx directory: {entex_dir}")
        entex_files = glob.glob(os.path.join(entex_dir, "**", "*.tsv"), recursive=True)
        logger.info(f"Found {len(entex_files)} TSV files in ENTEx directory")

        # Check metadata file
        if entex_metadata:
            logger.info(
                f"ENTEx metadata loaded with {len(entex_metadata.get('sample_lookup', {}))} samples"
            )
        else:
            logger.warning("No ENTEx metadata available")

        # Filter for gene quantification files for ENTEx
        if entex_metadata and "sample_lookup" in entex_metadata:
            # Group candidate files by donor and tissue.
            candidate_groups = {}
            for file_path in entex_files:
                file_name = os.path.basename(file_path)
                file_id = file_name.split(".")[0]

                if file_id in entex_metadata["sample_lookup"]:
                    sample_info = entex_metadata["sample_lookup"][file_id]
                    quant_method = sample_info.get("quantification_method", "").lower()
                    # Only consider gene quantification files (exclude transcript quantifications)
                    if "gene" not in quant_method or "transcript" in quant_method:
                        continue

                    # Create group key: use donor (or subject_id) and tissue (normalized to lower case)
                    donor = sample_info.get("donor") or sample_info.get("subject_id", "")
                    tissue = sample_info.get("tissue", "").lower()
                    group_key = (donor, tissue)
                    candidate_groups.setdefault(group_key, []).append(file_path)

            filtered_samples = []
            # Process each (donor, tissue) group
            for group_key, files in candidate_groups.items():
                if len(files) == 1:
                    # Only one candidate – keep it.
                    filtered_samples.append(files[0])
                    logger.debug(f"Keeping single file for group {group_key}: {files[0]}")
                else:
                    # More than one file: try to filter based on annotations.
                    preferred_files = []
                    hg38_files = []
                    for file_path in files:
                        file_name = os.path.basename(file_path)
                        file_id = file_name.split(".")[0]
                        sample_info = entex_metadata["sample_lookup"][file_id]
                        genome_annotation = sample_info.get("genome_annotation", "").upper()
                        genome_assembly = sample_info.get("genome_assembly", "").lower()
                        # Check if file meets both preferred criteria: V24 and hg38.
                        if genome_annotation == "V24" and genome_assembly == "hg38":
                            preferred_files.append(file_path)
                            logger.debug(f"Group {group_key}: file {file_id} meets V24 and hg38")
                        # Otherwise, check if it meets at least the hg38 criterion.
                        elif genome_assembly == "hg38":
                            hg38_files.append(file_path)
                            logger.debug(f"Group {group_key}: file {file_id} meets hg38 only")

                    if preferred_files:
                        # If at least one meets both criteria, use these.
                        filtered_samples.extend(preferred_files)
                    elif hg38_files:
                        # If none meet both, but some meet hg38, use these.
                        filtered_samples.extend(hg38_files)
                    else:
                        # Otherwise, fallback to keeping one candidate (choose the first).
                        filtered_samples.append(files[0])
                        logger.debug(
                            f"Group {group_key}: no file meets hg38; keeping first file: {files[0]}"
                        )

            logger.info(
                f"Selected {len(filtered_samples)} ENTEx gene quantification file(s) after grouping by donor and tissue"
            )
            tsv_files.extend(filtered_samples)

        else:
            # If no metadata filtering is possible, apply basic filtering to file paths
            gene_files = [
                f for f in entex_files if "gene" in f.lower() and "transcript" not in f.lower()
            ]
            logger.info(
                f"Selected {len(gene_files)} ENTEx files that appear to be gene quantifications"
            )
            tsv_files.extend(gene_files)

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
        sample_id = file_name.split(".")[0]

        # Rename TPM column to sample ID
        col_name = tpm_df.columns[0]
        tpm_df = tpm_df.rename(columns={col_name: sample_id})

        # Add to collections
        sample_dfs[sample_id] = tpm_df

        # Extract and store metadata with ENTEx support
        metadata = extract_encode_metadata(file_path, entex_metadata, metadata_dir)
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
    unified_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)

    # Create observation metadata DataFrame - ensuring it matches the sample IDs in unified_df
    # This is critical - the rows in obs_df must exactly match the columns in unified_df
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            # Make a copy of the metadata to avoid modifying the original
            metadata_copy = sample_metadata[sample_id].copy()

            # If there's a sample_id column, rename it to avoid conflict with index
            if "sample_id" in metadata_copy:
                metadata_copy["original_sample_id"] = metadata_copy.pop("sample_id")

            # Add the sample_id to each record
            metadata_copy["_sample_id"] = sample_id
            obs_data.append(metadata_copy)
        else:
            # If we somehow have a sample in the matrix but not in metadata,
            # create minimal metadata for it
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({"_sample_id": sample_id, "dataset": "ENCODE"})

    # Create obs_df with explicit index to ensure alignment
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index("_sample_id", inplace=True)

    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        logger.error(
            f"Metadata has {len(obs_df)} samples, expression matrix has {len(unified_df.columns)} samples"
        )

        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]
        logger.info(f"Fixed alignment - now using {len(obs_df)} samples")

    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df["gene_id"] = var_df.index

    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add original IDs to variable metadata
    var_df["original_ids"] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])))

    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize metadata
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "ENCODE", mappings)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        encode_metadata = load_dataset_specific_metadata(metadata_dir, "encode")
        if encode_metadata:
            # We'll apply this later to the AnnData object
            logger.info("Found ENCODE-specific metadata")

    # Dataset info for uns
    dataset_info = {
        "source": "ENCODE",
        "gencode_version": 24,
        "data_type": "RNA-seq",
        "expression_unit": "TPM",
        "samples": len(obs_df),
        "genes": len(var_df),
    }

    # Final dimension check before creating AnnData
    logger.info(f"Final dimensions check: {len(obs_df)} samples x {len(var_df)} genes")
    logger.info(f"Expression matrix shape: {unified_df.T.shape}")

    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        encode_metadata = load_dataset_specific_metadata(metadata_dir, "encode")
        if encode_metadata:
            adata = apply_dataset_specific_metadata(adata, encode_metadata)

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(
            f"Successfully processed ENCODE data with {adata.n_obs} samples and {adata.n_vars} genes"
        )

    return adata


def process_entex_data(metadata_file, output_file, metadata_dir=None):
    """
    Process ENTEx RNA-seq data into standardized AnnData using metadata.

    Parameters:
    -----------
    metadata_file : str
        Path to ENTEx metadata JSON file
    output_file : str
        Path to save the standardized AnnData
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    logger.info(f"Processing ENTEx data using metadata from {metadata_file}")

    # Load the metadata
    metadata = load_entex_metadata(metadata_file)
    if not metadata or not metadata.get("sample_lookup"):
        logger.error("No valid ENTEx metadata found")
        return None

    sample_lookup = metadata.get("sample_lookup", {})
    logger.info(f"Loaded metadata for {len(sample_lookup)} ENTEx samples")

    # Collect all file paths and process them
    sample_dfs = {}  # Data frames for each sample
    sample_metadata = {}  # Dictionary for metadata, keyed by sample_id
    all_gene_ids = set()  # Track all gene IDs
    processed_count = 0

    for file_id, sample_info in sample_lookup.items():
        # Get file path from metadata
        file_path = sample_info.get("file_path")

        if not file_path:
            logger.warning(f"No file path specified for {file_id}, attempting to locate file")
            # Try to find the file in standard locations
            possible_paths = [
                f"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/entex/{file_id}.tsv",
                f"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/entex/{file_id}.tsv",
            ]

            for path in possible_paths:
                if os.path.exists(path):
                    file_path = path
                    # Update the metadata for future runs
                    sample_info["file_path"] = file_path
                    logger.info(f"Found file for {file_id} at {file_path}")
                    break

            if not file_path:
                logger.warning(f"File not found for {file_id}")
                continue

        if not os.path.exists(file_path):
            logger.warning(f"File not found for {file_id}: {file_path}")
            continue

        # Read TPM data
        tpm_df = read_encode_tpm_file(file_path)
        if tpm_df is None:
            logger.warning(f"Could not read TPM data from {file_path}")
            continue

        # Rename TPM column to sample ID
        col_name = tpm_df.columns[0]
        tpm_df = tpm_df.rename(columns={col_name: file_id})

        # Add to collections
        sample_dfs[file_id] = tpm_df

        # Create metadata dictionary
        donor_id = sample_info.get("donor", sample_info.get("subject_id"))

        metadata_dict = {
            "original_sample_id": file_id,
            "tissue": sample_info.get("tissue", ""),
            "subject_id": donor_id,
            "assay_type": sample_info.get("assay_type", "total"),
            "genome_annotation": sample_info.get("genome_annotation", ""),
            "genome_assembly": sample_info.get("genome_assembly", ""),
            "data_source": "ENTEx",
            "data_type": "RNA-seq",
            "expression_unit": "TPM",
        }

        # Add donor information if available
        donor_map = metadata.get("donor_map", {})
        if donor_id and donor_id in donor_map:
            donor_info = donor_map[donor_id]
            metadata_dict["sex"] = donor_info.get("sex", "")
            metadata_dict["age"] = donor_info.get("age", "")

        # Store metadata
        sample_metadata[file_id] = metadata_dict

        # Track gene IDs
        all_gene_ids.update(tpm_df.index)
        processed_count += 1

        # Log progress
        if processed_count % 10 == 0:
            logger.info(f"Processed {processed_count} ENTEx samples")

    if not sample_dfs:
        logger.error("No valid ENTEx TPM data found")
        return None

    logger.info(f"Successfully processed {processed_count} ENTEx samples")

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
            if std_id in std_id_to_idx:
                gene_idx = std_id_to_idx[std_id]
                expr_matrix[gene_idx, sample_idx] = value

    # Create the DataFrame from the matrix
    unified_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)

    # Create observation metadata DataFrame
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            metadata_copy = sample_metadata[sample_id].copy()
            metadata_copy["_sample_id"] = sample_id
            obs_data.append(metadata_copy)
        else:
            # If we somehow have a sample in the matrix but not in metadata,
            # create minimal metadata for it
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({"_sample_id": sample_id, "dataset": "ENTEx"})

    # Create obs_df with explicit index
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index("_sample_id", inplace=True)

    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        # Fix alignment
        obs_df = obs_df.loc[unified_df.columns]

    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df["gene_id"] = var_df.index

    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add original IDs to variable metadata
    var_df["original_ids"] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])))

    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize metadata
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "ENTEx", mappings)

    # Dataset info for uns
    dataset_info = {
        "source": "ENTEx",
        "gencode_version": 24,
        "data_type": "RNA-seq",
        "expression_unit": "TPM",
        "samples": len(obs_df),
        "genes": len(var_df),
        "tissues": len(set(obs_df["tissue"])) if "tissue" in obs_df.columns else 0,
        "donors": len(set(obs_df["subject_id"])) if "subject_id" in obs_df.columns else 0,
    }

    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        entex_metadata = load_dataset_specific_metadata(metadata_dir, "entex")
        if entex_metadata:
            adata = apply_dataset_specific_metadata(adata, entex_metadata)

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(
            f"Successfully processed ENTEx data with {adata.n_obs} samples and {adata.n_vars} genes"
        )

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
        sample_attrs = pd.read_csv(GTEx_SAMPLE_ATTRS, sep="\t")

        # Load subject attributes
        logger.info(f"Loading subject attributes from {GTEx_SUBJECT_ATTRS}")
        subject_attrs = pd.read_csv(GTEx_SUBJECT_ATTRS, sep="\t")

        # Extract SUBJID from SAMPID
        sample_attrs["SUBJID"] = sample_attrs["SAMPID"].str.extract(r"(GTEX-[A-Z0-9]+)")

        # Merge sample and subject attributes
        merged_attrs = pd.merge(sample_attrs, subject_attrs, on="SUBJID", how="left")

        # Standardize column names
        column_mapping = {
            "SAMPID": "sample_id",
            "SUBJID": "subject_id",
            "SMTSD": "tissue",  # Detailed tissue site
            "SMTS": "broad_tissue",  # Broad tissue type
            "SEX": "sex",
            "AGE": "age",
            "SMGEBTCH": "batch",
            "SMRIN": "rna_integrity_number",
            "SMTSISCH": "ischemic_time",
            "SMNABTCH": "array_batch",
        }

        # Rename columns that exist in the DataFrame
        existing_cols = [col for col in column_mapping.keys() if col in merged_attrs.columns]
        merged_attrs = merged_attrs.rename(
            columns={col: column_mapping[col] for col in existing_cols}
        )

        # Ensure we have a tissue column
        if "tissue" not in merged_attrs.columns and "broad_tissue" in merged_attrs.columns:
            merged_attrs["tissue"] = merged_attrs["broad_tissue"]
            logger.info("Using broad_tissue as tissue column")

        # Set sample ID as index
        merged_attrs.set_index("sample_id", inplace=True)

        logger.info(f"Loaded metadata for {len(merged_attrs)} GTEx samples")

        # Extract tissue statistics
        if "tissue" in merged_attrs.columns:
            tissue_counts = merged_attrs["tissue"].value_counts()
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
        with gzip.open(file_path, "rt") as f:
            # Skip version line
            next(f)
            # Read dimensions line
            dims = next(f).strip().split()
            n_genes, n_samples = int(dims[0]), int(dims[1])
            logger.info(f"GTEx file dimensions: {n_genes} genes, {n_samples} samples")

            # Read column headers
            header_line = next(f).strip().split("\t")
            sample_ids = header_line[2:]  # Skip Name and Description columns

        # Now read the data in chunks
        chunk_size = 1000  # Adjust based on memory constraints

        logger.info(f"Reading GTEx data in chunks of {chunk_size} genes")

        # Initialize the gene list and expression list
        all_genes = []
        all_data = []

        with gzip.open(file_path, "rt") as f:
            # Skip header lines
            next(f)  # Version
            next(f)  # Dimensions
            next(f)  # Column headers

            # Read data line by line
            line_count = 0
            for line in f:
                line_count += 1

                fields = line.strip().split("\t")
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
        logger.info(
            f"Loaded expression data with {data_df.shape[0]} genes and {data_df.shape[1]} samples"
        )

        return data_df

    except Exception as e:
        logger.error(f"Error loading GTEx expression data: {e}")
        import traceback

        logger.error(traceback.format_exc())
        return pd.DataFrame()


def process_gtex_data(input_file, output_file, metadata_dir=None):
    """
    Process GTEx RNA-seq data into standardized AnnData - enhanced version
    with support for JSON configuration files.

    Parameters:
    -----------
    input_file : str
        Path to GTEx expression file (GCT format)
    output_file : str
        Path to save the standardized AnnData
    metadata_dir : str, optional
        Directory containing metadata JSON files

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
        metadata_df["sample_id"] = metadata_df.index

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

    # Create vectorized mapping arrays
    logger.info("Creating vectorized mappings for faster processing")
    orig_to_std_idx = {}  # Map from original gene index to standardized gene index
    
    for gene_idx, gene_id in enumerate(expr_df.index):
        std_id = gene_id_mapping[gene_id]
        if std_id in std_id_to_idx:
            std_idx = std_id_to_idx[std_id]
            orig_to_std_idx[gene_idx] = std_idx
    
    # Process all samples at once using numpy operations
    logger.info("Processing all samples with vectorized operations")
    
    # Convert to numpy array for faster processing
    expr_values = expr_df.values
    
    # Calculate batch size based on memory constraints
    # Process ~1000 samples at a time to avoid memory issues
    batch_size = 1000
    batches = (sample_count + batch_size - 1) // batch_size
    
    for batch in range(batches):
        start_idx = batch * batch_size
        end_idx = min((batch + 1) * batch_size, sample_count)
        logger.info(f"Processing sample batch {batch+1}/{batches} (samples {start_idx+1}-{end_idx})")
        
        # Extract batch of samples
        batch_samples = expr_values[:, start_idx:end_idx]
        
        # Update expression matrix for each gene mapping
        for orig_idx, std_idx in orig_to_std_idx.items():
            # For each original gene, update the corresponding standardized gene
            expr_matrix[std_idx, start_idx:end_idx] = np.maximum(
                expr_matrix[std_idx, start_idx:end_idx], 
                batch_samples[orig_idx, :]
            )
    

    # Create DataFrame with the results
    std_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=common_samples)

    # Free memory
    del expr_df
    del expr_matrix

    # Create variable (gene) metadata DataFrame
    var_df = pd.DataFrame(index=unique_std_ids)
    var_df["gene_id"] = var_df.index

    # Add original ID mapping information
    var_df["original_ids"] = ""

    # Create a reverse mapping from std_id to original ids
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add the original IDs as a concatenated string
    for std_id in var_df.index:
        if std_id in reverse_mapping:
            var_df.loc[std_id, "original_ids"] = ";".join(reverse_mapping[std_id])

    # Add GENCODE annotations
    logger.info("Adding GENCODE annotations")
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize observation metadata
    logger.info("Standardizing metadata")
    mappings = load_mappings()
    obs_df = standardize_metadata(metadata_df, "GTEx", mappings)

    # Fix index/column name conflict - rename the column if it exists
    if "sample_id" in obs_df.columns:
        logger.info(
            "Renaming 'sample_id' column to 'original_sample_id' to avoid conflict with index"
        )
        obs_df = obs_df.rename(columns={"sample_id": "original_sample_id"})

    # Add missing standard fields
    if "data_type" not in obs_df.columns:
        obs_df["data_type"] = "RNA-seq"
    if "expression_unit" not in obs_df.columns:
        obs_df["expression_unit"] = "TPM"

    # Dataset info for uns
    dataset_info = {
        "source": "GTEx",
        "version": "v10",  # Updated to v10
        "gencode_version": 24,
        "data_type": "RNA-seq",
        "expression_unit": "TPM",
        "samples": len(obs_df),
        "genes": len(var_df),
        "tissue_count": len(obs_df["tissue"].unique()) if "tissue" in obs_df.columns else 0,
        "subject_count": (
            len(obs_df["subject_id"].unique()) if "subject_id" in obs_df.columns else 0
        ),
    }

    # Create AnnData object
    logger.info("Creating AnnData object")
    logger.info(f"var_df shape before AnnData creation: {var_df.shape}")
    logger.info(f"var_df columns: {list(var_df.columns)}")

    adata = create_standard_anndata(std_expr_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        gtex_metadata = load_dataset_specific_metadata(metadata_dir, "gtex")
        if gtex_metadata:
            adata = apply_dataset_specific_metadata(adata, gtex_metadata)

    logger.info(f"adata.var shape after creation: {adata.var.shape}")
    logger.info(f"adata.var columns after creation: {list(adata.var.columns)}")

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        # Verify one more time after loading saved file
        verif_adata = ad.read_h5ad(output_file)
        logger.info(f"Verification - var shape after reload: {verif_adata.var.shape}")
        logger.info(f"Verification - var columns after reload: {list(verif_adata.var.columns)}")

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
            with open(file_path, "r") as f:
                first_line = f.readline().strip()

            # Determine delimiter based on content
            if "\t" in first_line:
                delimiter = "\t"
            elif "," in first_line:
                delimiter = ","
            else:
                # Default to tab if can't determine
                delimiter = "\t"

            logger.debug(f"Using delimiter '{delimiter}' for file {file_path}")

            # Now read with the appropriate delimiter
            df = pd.read_csv(file_path, delimiter=delimiter)

            # Check if empty after reading
            if df.empty:
                logger.warning(f"Empty data after reading: {file_path}")
                return False

            logger.debug(f"Read data with columns: {df.columns.tolist()}")

            # Check for expected columns
            if "gene_id" in df.columns and "TPM" in df.columns:
                # File has the expected structure, but TPM is actually normalized microarray intensity
                # Rename TPM to normalized_intensity if processing ADNI data
                if "ADNI" in file_path:
                    df = df.rename(columns={"TPM": "normalized_intensity"})
                    expr_col = "normalized_intensity"
                    expr_type = "normalized_intensity"
                else:  # For other datasets like MAGE, keep TPM
                    expr_col = "TPM"
                    expr_type = "TPM"

                # Remove FPKM column if it exists (as per your requirement)
                if "FPKM" in df.columns:
                    df = df.drop("FPKM", axis=1)

                # Create simplified DataFrame with gene_id and expression values
                simple_df = pd.DataFrame()
                simple_df[expr_col] = df[expr_col]
                simple_df.index = df["gene_id"]

                # Store the data
                sample_dfs[sample_id] = simple_df

                # Collect gene IDs
                all_gene_ids.update(df["gene_id"])

                # Check if gene IDs are Ensembl IDs
                is_ensembl = (
                    all(str(gene_id).startswith("ENSG") for gene_id in df["gene_id"].iloc[:5])
                    if len(df) >= 5
                    else False
                )

                # Prepare metadata
                metadata = {
                    "sample_id": sample_id,
                    "donor_id": donor_id,
                    "subject_id": donor_id,
                    "tissue": tissue,
                    "dataset": "MAGE",
                    "data_type": "RNA-seq",
                    "expression_unit": expr_type,  # This will be TPM for RNA-seq or normalized_intensity for microarray
                    "is_gencode": "gencode" in file_name.lower(),
                    "is_ensembl": is_ensembl,
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
                    if "gene" in col.lower() and "id" in col.lower():
                        gene_id_col = col
                        break

                # Look for expression column - try common names
                expr_col = None
                for col in df.columns:
                    if col.upper() in [
                        "TPM",
                        "FPKM",
                        "COUNTS",
                        "EXPRESSION",
                        "NORMALIZED_INTENSITY",
                    ]:
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
                    is_ensembl = (
                        all(str(gene_id).startswith("ENSG") for gene_id in df[gene_id_col].iloc[:5])
                        if len(df) >= 5
                        else False
                    )

                    # Determine expression type based on column name
                    if "TPM" in expr_col.upper():
                        expr_type = "TPM"
                    elif "FPKM" in expr_col.upper():
                        expr_type = "FPKM"
                    elif "COUNT" in expr_col.upper():
                        expr_type = "counts"
                    elif "NORMALIZED" in expr_col.upper() or "INTENSITY" in expr_col.upper():
                        expr_type = "normalized_intensity"
                    else:
                        expr_type = "TPM"  # Default

                    # Prepare metadata
                    metadata = {
                        "sample_id": sample_id,
                        "donor_id": donor_id,
                        "subject_id": donor_id,
                        "tissue": tissue,
                        "dataset": "MAGE",
                        "data_type": "RNA-seq",
                        "expression_unit": expr_type,
                        "is_gencode": "gencode" in file_name.lower(),
                        "is_ensembl": is_ensembl,
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


def process_mage_data(input_dir, output_file, metadata_dir=None):
    """
    Process MAGE RNA-seq data into standardized AnnData - updated to use
    JSON configurations.

    Parameters:
    -----------
    input_dir : str
        Directory containing MAGE expression files
    output_file : str
        Path to save the standardized AnnData
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
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
    donor_dirs = glob.glob(os.path.join(input_dir, "NA*")) + glob.glob(
        os.path.join(input_dir, "HG*")
    )

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
            csv_files = [f for f in os.listdir(lymphoblast_dir) if f.endswith(".csv")]

            # Separate original and gencode files
            original_files = [f for f in csv_files if "_gencode_" not in f]
            gencode_files = [f for f in csv_files if "_gencode_" in f]

            # Process original files first (prioritize them)
            for file_name in original_files:
                file_path = os.path.join(lymphoblast_dir, file_name)
                if process_mage_dir(
                    file_path, donor_id, tissue_name, sample_dfs, sample_metadata, all_gene_ids
                ):
                    processed_samples += 1

            # Only process gencode files if no original files were successfully processed
            if not any(donor_id in sample_id for sample_id in sample_dfs.keys()):
                for file_name in gencode_files:
                    file_path = os.path.join(lymphoblast_dir, file_name)
                    if process_mage_dir(
                        file_path, donor_id, tissue_name, sample_dfs, sample_metadata, all_gene_ids
                    ):
                        processed_samples += 1
        else:
            # Try to find any tissue directories
            tissue_dirs = [d for d in glob.glob(os.path.join(donor_dir, "*")) if os.path.isdir(d)]

            for tissue_dir in tissue_dirs:
                tissue_name = os.path.basename(tissue_dir)

                # Find all CSV files in the directory
                csv_files = [f for f in os.listdir(tissue_dir) if f.endswith(".csv")]

                # Separate original and gencode files
                original_files = [f for f in csv_files if "_gencode_" not in f]
                gencode_files = [f for f in csv_files if "_gencode_" in f]

                # Process original files first
                for file_name in original_files:
                    file_path = os.path.join(tissue_dir, file_name)
                    if process_mage_dir(
                        file_path, donor_id, tissue_name, sample_dfs, sample_metadata, all_gene_ids
                    ):
                        processed_samples += 1

                # Only process gencode files if no original files were successfully processed for this tissue
                tissue_sample_ids = [
                    sample_id
                    for sample_id in sample_dfs.keys()
                    if donor_id in sample_id and tissue_name in sample_id
                ]
                if not tissue_sample_ids:
                    for file_name in gencode_files:
                        file_path = os.path.join(tissue_dir, file_name)
                        if process_mage_dir(
                            file_path,
                            donor_id,
                            tissue_name,
                            sample_dfs,
                            sample_metadata,
                            all_gene_ids,
                        ):
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
    unified_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)

    # Create observation metadata DataFrame
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            metadata_copy = sample_metadata[sample_id].copy()
            metadata_copy["_sample_id"] = sample_id
            obs_data.append(metadata_copy)
        else:
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({"_sample_id": sample_id, "dataset": "MAGE"})

    # Create obs_df with explicit index
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index("_sample_id", inplace=True)

    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]

    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df["gene_id"] = var_df.index

    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add original IDs to variable metadata
    var_df["original_ids"] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])))

    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize metadata
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "MAGE", mappings)

    # Dataset info for uns
    dataset_info = {
        "source": "MAGE",
        "gencode_version": 24,  # Assuming same version as other datasets
        "data_type": "RNA-seq",
        "expression_unit": "TPM",  # Update if different
        "samples": len(obs_df),
        "genes": len(var_df),
        "donor_count": processed_donors,
        "tissue_count": len(set(obs_df["tissue"] if "tissue" in obs_df.columns else [])),
    }

    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        mage_metadata = load_dataset_specific_metadata(metadata_dir, "mage")
        if mage_metadata:
            adata = apply_dataset_specific_metadata(adata, mage_metadata)

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(
            f"Successfully processed MAGE data with {adata.n_obs} samples and {adata.n_vars} genes"
        )

    return adata


# ====================================================================================
# ADNI-specific Processing
# ====================================================================================

def preprocess_adni_file(file_path):
    """Preprocess ADNI file to fix escaped tab characters."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Check if file has escaped tabs
    if '\\t' in content:
        logger.info(f"Fixing escaped tabs in {file_path}")
        # Replace escaped tabs with actual tabs
        fixed_content = content.replace('\\t', '\t')
        
        # Create a temporary fixed file
        fixed_path = file_path + '.fixed'
        with open(fixed_path, 'w') as f:
            f.write(fixed_content)
        
        return fixed_path
    
    return file_path


def preprocess_adni_file(file_path):
    """Preprocess ADNI file to fix escaped tabs."""
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check if file has escaped tabs
        if '\t' in content:
            logger.info(f"Fixing escaped tabs in {file_path}")
            # Replace escaped tabs with actual tabs
            fixed_content = content.replace('\t', '	')
            
            # Create a temporary fixed file
            fixed_path = file_path + '.fixed'
            with open(fixed_path, 'w') as f:
                f.write(fixed_content)
            
            return fixed_path
    except Exception as e:
        logger.warning(f"Error preprocessing file {file_path}: {e}")
    
    return file_path

def process_adni_data(input_dir, output_file, dict_file=None, metadata_dir=None):
    """
    Process ADNI gene expression data into standardized AnnData with JSON configuration support.

    Parameters:
    -----------
    input_dir : str
        Directory containing ADNI sample directories with expression files
    output_file : str
        Path to save the standardized AnnData
    dict_file : str, optional
        Path to gene ID mapping dictionary file
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
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

                # Preprocess the file to fix escaped tabs
                fixed_file_path = preprocess_adni_file(file_path)
                try:
                    # Try reading with tab delimiter first
                    df = pd.read_csv(fixed_file_path, sep='\t')
                    # If successful, clean up temporary file if it was created
                    if fixed_file_path != file_path and os.path.exists(fixed_file_path):
                        os.remove(fixed_file_path)
                except Exception as e:
                    # If tab delimiter fails, try comma
                    logger.warning(f"Failed to parse with tab delimiter: {e}")
                    df = pd.read_csv(fixed_file_path)
                    # Clean up temporary file
                    if fixed_file_path != file_path and os.path.exists(fixed_file_path):
                        os.remove(fixed_file_path)

                gene_id_col, tpm_col = identify_columns(df)                
                

                if gene_id_col is None or tpm_col is None:
                    logger.warning(f"Could not identify gene_id and TPM columns in {file_path}")
                    logger.warning(f"Available columns: {list(df.columns)}")
                    continue  # Skip this file

                # Create expression DataFrame using the identified columns
                expr_df = pd.DataFrame(df[tpm_col].values, index=df[gene_id_col])

                # Store the data
                sample_dfs[sample_id] = expr_df

                # Collect gene IDs
                all_gene_ids.update(df[gene_id_col])

                metadata = {
                    "sample_id": sample_id,
                    "subject_id": subject_id,
                    "dataset": "ADNI",
                    "data_type": "Microarray",
                    "expression_unit": "Normalized intensity",
                    "tissue": "blood",
                    "platform": "Affymetrix Human Genome U219 Array",
                    "processing": "gencode_v24" if "gencode_v24" in file_name else "unknown",
                }
                # Store metadata
                sample_metadata[sample_id] = metadata
                processed_samples += 1

                logger.debug(f"Successfully processed {file_path} with {len(expr_df)} genes")

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
    unified_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)

    # Create observation metadata DataFrame
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            metadata_copy = sample_metadata[sample_id].copy()
            metadata_copy["_sample_id"] = sample_id
            obs_data.append(metadata_copy)
        else:
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({"_sample_id": sample_id, "dataset": "ADNI"})

    # Create obs_df with explicit index
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index("_sample_id", inplace=True)

    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]

    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df["gene_id"] = var_df.index

    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add original IDs to variable metadata
    var_df["original_ids"] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])))

    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize metadata
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "ADNI", mappings)

    # Dataset info for uns
    dataset_info = {
        "source": "ADNI",
        "data_type": "RNA-seq",  # Updated from microarray since these appear to be RNA-seq TPM values
        "gencode_version": 24,
        "expression_unit": "TPM",
        "samples": len(obs_df),
        "genes": len(var_df),
        "subject_count": processed_subjects,
    }

    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        adni_metadata = load_dataset_specific_metadata(metadata_dir, "adni")
        if adni_metadata:
            adata = apply_dataset_specific_metadata(adata, adni_metadata)

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(
            f"Successfully processed ADNI data with {adata.n_obs} samples and {adata.n_vars} genes"
        )

    return adata


# ====================================================================================
# Main Processing Pipeline
# ====================================================================================
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Multi-Dataset Standardization Pipeline")
    parser.add_argument("--encode-dir", help="Directory containing ENCODE cell line TPM files")
    parser.add_argument("--encode-entex-dir", help="Directory containing ENCODE ENTEx TPM files")
    parser.add_argument("--gtex-file", help="Path to GTEx expression file (GCT format)")
    parser.add_argument("--mage-dir", help="Directory containing MAGE expression files")
    parser.add_argument(
        "--adni-dir", help="Directory containing ADNI sample directories with expression files"
    )
    parser.add_argument(
        "--metadata-dir",
        default=str(DEFAULT_METADATA_DIR),
        help="Directory containing metadata JSON files",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Output directory for standardized files",
    )

    # Add ENTEx-specific arguments
    parser.add_argument("--entex-metadata-file", help="Path to ENTEx metadata JSON file")
    parser.add_argument(
        "--dataset-type",
        choices=["encode", "entex", "gtex", "mage", "adni"],
        help="Type of dataset to process when using metadata",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    start_time = time.time()

    # Process ENTEx data if metadata file is provided and dataset-type is 'entex'
    if args.entex_metadata_file and args.dataset_type == "entex":
        entex_output = os.path.join(args.output_dir, "entex_standardized_v1.h5ad")
        logger.info(f"Processing ENTEx data from metadata file: {args.entex_metadata_file}")
        entex_data = process_entex_data(args.entex_metadata_file, entex_output, args.metadata_dir)
        if entex_data is not None:
            logger.info(
                f"Successfully processed ENTEx data: {entex_data.n_obs} samples and {entex_data.n_vars} genes"
            )
            # Check gene statistics
            logger.info("ENTEx gene statistics:")
            gencode_stats = entex_data.var["mapping_source"].value_counts()
            for source, count in gencode_stats.items():
                logger.info(f"  {source}: {count} genes ({count/entex_data.n_vars:.1%})")

    # Process standard datasets if their respective inputs are provided
    elif args.encode_dir:
        encode_output = os.path.join(args.output_dir, "encode_standardized_v1.h5ad")
        logger.info(f"Processing ENCODE data from {args.encode_dir}")
        encode_data = process_encode_data(
            args.encode_dir,
            args.encode_entex_dir,
            encode_output,
            args.entex_metadata_file,
            args.metadata_dir,
        )
        if encode_data is not None:
            logger.info(
                f"Successfully processed ENCODE data: {encode_data.n_obs} samples and {encode_data.n_vars} genes"
            )
            # Check gene statistics
            logger.info("ENCODE gene statistics:")
            gencode_stats = encode_data.var["mapping_source"].value_counts()
            for source, count in gencode_stats.items():
                logger.info(f"  {source}: {count} genes ({count/encode_data.n_vars:.1%})")

    if args.gtex_file:
        gtex_output = os.path.join(args.output_dir, "gtex_standardized_v1.h5ad")
        logger.info(f"Processing GTEx data from {args.gtex_file}")
        gtex_data = process_gtex_data(args.gtex_file, gtex_output, args.metadata_dir)
        if gtex_data is not None:
            logger.info(
                f"Successfully processed GTEx data: {gtex_data.n_obs} samples, {gtex_data.n_vars} genes"
            )
            # Check gene statistics
            logger.info("GTEx gene statistics:")
            gencode_stats = gtex_data.var["mapping_source"].value_counts()
            for source, count in gencode_stats.items():
                logger.info(f"  {source}: {count} genes ({count/gtex_data.n_vars:.1%})")

    if args.mage_dir:
        mage_output = os.path.join(args.output_dir, "mage_standardized_v1.h5ad")
        logger.info(f"Processing MAGE data from {args.mage_dir}")
        mage_data = process_mage_data(args.mage_dir, mage_output, args.metadata_dir)
        if mage_data is not None:
            logger.info(
                f"Successfully processed MAGE data: {mage_data.n_obs} samples, {mage_data.n_vars} genes"
            )

    # Process ADNI data
    if args.adni_dir:
        adni_output = os.path.join(args.output_dir, "adni_standardized_v1.h5ad")
        logger.info(f"Processing ADNI data from directory: {args.adni_dir}")
        adni_data = process_adni_data(args.adni_dir, adni_output, None, args.metadata_dir)
        if adni_data is not None:
            logger.info(
                f"Successfully processed ADNI data: {adni_data.n_obs} samples, {adni_data.n_vars} genes"
            )

    # Check for common genes if multiple datasets were processed
    processed_datasets = []
    if args.encode_dir and "encode_data" in locals() and encode_data is not None:
        processed_datasets.append(("ENCODE", encode_data))
    if args.gtex_file and "gtex_data" in locals() and gtex_data is not None:
        processed_datasets.append(("GTEx", gtex_data))
    if args.mage_dir and "mage_data" in locals() and mage_data is not None:
        processed_datasets.append(("MAGE", mage_data))
    if args.adni_dir and "adni_data" in locals() and adni_data is not None:
        processed_datasets.append(("ADNI", adni_data))

    if len(processed_datasets) > 1:
        logger.info("Analyzing gene overlap between datasets:")

        # Compare each pair of datasets
        for i in range(len(processed_datasets)):
            for j in range(i + 1, len(processed_datasets)):
                name1, data1 = processed_datasets[i]
                name2, data2 = processed_datasets[j]

                genes1 = set(data1.var_names)
                genes2 = set(data2.var_names)
                common_genes = genes1.intersection(genes2)

                logger.info(f"  {name1} ({len(genes1)} genes) vs {name2} ({len(genes2)} genes):")
                logger.info(f"    Common genes: {len(common_genes)}")
                logger.info(
                    f"    Overlap: {len(common_genes)/len(genes1):.1%} of {name1}, {len(common_genes)/len(genes2):.1%} of {name2}"
                )

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

                logger.info(
                    f"    Top 10% expressed genes: {name1} ({len(top_genes1)}), {name2} ({len(top_genes2)})"
                )
                logger.info(f"    Common top genes: {len(common_top)}")
                logger.info(
                    f"    Top gene overlap: {len(common_top)/len(top_genes1):.1%} of top {name1}, {len(common_top)/len(top_genes2):.1%} of top {name2}"
                )

        # If both ENCODE and GTEx were processed, create a gene compatibility report
        if (
            "encode_data" in locals()
            and "gtex_data" in locals()
            and encode_data is not None
            and gtex_data is not None
        ):
            report_file = os.path.join(args.output_dir, "gene_compatibility_report.csv")
            logger.info(f"Creating gene compatibility report: {report_file}")

            # Create DataFrame with gene information from both datasets
            encode_genes = pd.DataFrame(index=encode_data.var_names)
            encode_genes["gene_name"] = encode_data.var["gene_name"]
            encode_genes["gene_type"] = encode_data.var["gene_type"]
            encode_genes["in_encode"] = True

            gtex_genes = pd.DataFrame(index=gtex_data.var_names)
            gtex_genes["gene_name"] = gtex_data.var["gene_name"]
            gtex_genes["gene_type"] = gtex_data.var["gene_type"]
            gtex_genes["in_gtex"] = True

            # Merge the DataFrames
            all_genes = encode_genes.join(
                gtex_genes, how="outer", lsuffix="_encode", rsuffix="_gtex"
            )
            all_genes["in_encode"] = all_genes["in_encode"].fillna(False)
            all_genes["in_gtex"] = all_genes["in_gtex"].fillna(False)

            # Fill in missing gene names using the other dataset
            all_genes["gene_name"] = all_genes["gene_name_encode"].combine_first(
                all_genes["gene_name_gtex"]
            )
            all_genes["gene_type"] = all_genes["gene_type_encode"].combine_first(
                all_genes["gene_type_gtex"]
            )

            # Calculate average expression in each dataset
            all_genes["mean_expr_encode"] = 0.0
            all_genes["mean_expr_gtex"] = 0.0

            for gene in all_genes.index:
                if gene in encode_data.var_names:
                    gene_idx = encode_data.var_names.get_loc(gene)
                    all_genes.loc[gene, "mean_expr_encode"] = encode_data.X[:, gene_idx].mean()
                if gene in gtex_data.var_names:
                    gene_idx = gtex_data.var_names.get_loc(gene)
                    all_genes.loc[gene, "mean_expr_gtex"] = gtex_data.X[:, gene_idx].mean()

            # Save the report
            all_genes.to_csv(report_file)
            logger.info(f"Saved gene compatibility report with {len(all_genes)} genes")

    total_time = time.time() - start_time
    logger.info(f"Total processing time: {total_time:.2f} seconds")

    # Print validation summary for all processed datasets
    validation_summary = {}
    if processed_datasets:
        logger.info("=== Metadata Validation Summary ===")
        for name, data in processed_datasets:
            if "metadata_validation" in data.uns:
                report = data.uns["metadata_validation"]
                validation_summary[name] = report

                # Print summary
                logger.info(f"Dataset: {name}")
                logger.info(f"  Samples: {report['sample_count']}")
                logger.info(f"  Status: {report['validation_status']}")

                if report["unmapped_tissues"]:
                    logger.info(f"  Unmapped tissues: {len(report['unmapped_tissues'])}")

                if report["unmapped_assays"]:
                    logger.info(f"  Unmapped assays: {len(report['unmapped_assays'])}")

                if report["unmapped_ages"]:
                    logger.info(f"  Unmapped ages: {len(report['unmapped_ages'])}")

                logger.info(
                    f"  Missing fields: {', '.join(report['missing_fields']) if report['missing_fields'] else 'None'}"
                )

        logger.info("================================")


if __name__ == "__main__":
    main()
