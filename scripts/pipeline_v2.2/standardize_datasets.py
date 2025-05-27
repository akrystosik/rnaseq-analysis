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
    ensure_serializable,
    standardize_ensembl_id,
    load_gencode_mapping, # Still import here, but call once in main
    add_gencode_annotations,
    standardize_metadata,
    validate_metadata,
    prepare_for_anndata,
    load_mappings,
    load_json_mapping,
    load_dataset_specific_metadata,
    apply_dataset_specific_metadata,
    NIH_RACE_CATEGORIES_LOWER,
    NIH_ETHNICITY_CATEGORIES_LOWER,
    SRE_BASE_CATEGORIES_LOWER
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


def identify_columns(df):
    """Identify gene_id and expression value columns across different file formats"""
    gene_id_col = None
    expr_val_col = None

    columns_lower = [str(col).lower() for col in df.columns]

    # Priority for gene identifier column
    gene_id_candidates = ["gene_id", "target_id", "gene", "name", "id"] # Ordered by preference
    for candidate in gene_id_candidates:
        if candidate in columns_lower:
            gene_id_col = df.columns[columns_lower.index(candidate)]
            logger.debug(f"Identified gene ID column: '{gene_id_col}' (matched '{candidate}')")
            break

    if gene_id_col is None and len(df.columns) >=1: # Fallback if no specific name matches
        gene_id_col = df.columns[0]
        logger.debug(f"Using fallback first column '{gene_id_col}' as gene ID column.")


    # Priority for expression value column
    expr_val_candidates = ["tpm", "fpkm", "counts", "normalized_intensity", "expression", "value"]
    for candidate in expr_val_candidates:
        if candidate in columns_lower:
            expr_val_col = df.columns[columns_lower.index(candidate)]
            logger.debug(f"Identified expression value column: '{expr_val_col}' (matched '{candidate}')")
            break

    if expr_val_col is None and len(df.columns) >= 2: # Fallback if no specific name matches
        # Try to find a numeric column if gene_id_col is the first one
        if gene_id_col == df.columns[0] and pd.api.types.is_numeric_dtype(df[df.columns[1]]):
            expr_val_col = df.columns[1]
            logger.debug(f"Using fallback second column '{expr_val_col}' as expression value column.")
        elif len(df.select_dtypes(include=np.number).columns) > 0 : # Generic numeric column
            expr_val_col = df.select_dtypes(include=np.number).columns[0]
            logger.debug(f"Using first numeric column '{expr_val_col}' as expression value column.")


    return gene_id_col, expr_val_col

def get_encode_cell_info(metadata_dir=None):
    """Get ENCODE cell line info from JSON."""
    if not metadata_dir:
        metadata_dir = DEFAULT_METADATA_DIR

    json_path = Path(metadata_dir) / "encode_metadata.json"
    if json_path.exists():
        try:
            with open(json_path, "r") as f:
                metadata_config = json.load(f)
            if "dataset_info" in metadata_config and "cell_lines" in metadata_config["dataset_info"]:
                logger.info(f"Successfully loaded ENCODE cell line info from {json_path}")
                return metadata_config["dataset_info"]["cell_lines"]
            else:
                logger.warning(f"'dataset_info.cell_lines' not found in {json_path}. No ENCODE cell line specific info will be used beyond defaults.")
        except Exception as e:
            logger.error(f"Error loading ENCODE cell info from {json_path}: {e}")
    else:
        logger.error(f"ENCODE metadata file not found: {json_path}. Cannot retrieve cell line info.")
    return {}




def create_standard_anndata(data_df, obs_df, var_df, dataset_specific_cfg, dataset_name_for_log):
    """
    Create standardized AnnData object with consistent structure.
    Relies on dataset_specific_cfg for .uns population.
    """
    logger.info(f"Creating AnnData for {dataset_name_for_log}")

    obs_df, var_df, data_df = prepare_for_anndata(obs_df.copy(), var_df.copy(), data_df.copy())

    logger.debug(
        f"Creating AnnData with X shape={data_df.shape}, obs shape={obs_df.shape}, var shape={var_df.shape}"
    )
    try:
        # Try sparse matrix creation first for X, especially for RNA-seq TPM/counts
        if not sp.issparse(data_df.values): # Check if not already sparse
             adata = ad.AnnData(X=sp.csr_matrix(data_df.values), obs=obs_df, var=var_df)
        else:
             adata = ad.AnnData(X=data_df.values, obs=obs_df, var=var_df) # If already sparse or small enough
    except Exception as e:
        logger.error(f"Error creating AnnData object for {dataset_name_for_log}: {e}")
        try: # Fallback to dense if sparse failed for some reason (should be rare for typical numeric data)
            logger.warning("Retrying AnnData creation with dense X matrix due to previous error.")
            adata = ad.AnnData(X=data_df.values, obs=obs_df, var=var_df)
        except Exception as e_dense:
            logger.error(f"Dense AnnData creation also failed for {dataset_name_for_log}: {e_dense}")
            raise e_dense

    # Populate adata.uns using dataset_specific_cfg
    if dataset_specific_cfg:
        for key, value in dataset_specific_cfg.items():
            if key == "obs_columns": # These are defaults for obs_df, not for direct .uns storage here.
                continue
            adata.uns[key] = value # Will be serialized by save_anndata -> ensure_serializable

    # Ensure a basic dataset_info structure exists in uns, even if cfg is minimal
    if 'dataset_info' not in adata.uns:
        adata.uns['dataset_info'] = {}

    # Populate calculated/summary fields into dataset_info
    adata.uns['dataset_info']['source_dataset_label'] = dataset_specific_cfg.get("label", dataset_name_for_log)
    adata.uns['dataset_info']['processed_version'] = dataset_specific_cfg.get("version", "1.0") # Our internal pipeline version
    adata.uns['dataset_info']['original_dataset_version'] = dataset_specific_cfg.get("dataset_version", "unknown")
    adata.uns['dataset_info']['sample_count'] = adata.n_obs
    adata.uns['dataset_info']['gene_count'] = adata.n_vars
    adata.uns['dataset_info']['data_type'] = obs_df['data_type'].unique()[0] if 'data_type' in obs_df and obs_df['data_type'].nunique() == 1 else dataset_specific_cfg.get('data_type', 'unknown')
    adata.uns['dataset_info']['expression_unit'] = obs_df['expression_unit'].unique()[0] if 'expression_unit' in obs_df and obs_df['expression_unit'].nunique() == 1 else dataset_specific_cfg.get('obs_columns',{}).get('expression_unit','unknown')


    # Set project-wide harmonized standards directly in .uns
    adata.uns["harmonized_reference_genome"] = "hg38"
    adata.uns["harmonized_gencode_version"] = "v24"

    # Original versions from config are already top-level in adata.uns if they existed in dataset_specific_cfg
    # If not, set them as 'unknown' or based on config.
    adata.uns.setdefault("original_reference_genome", dataset_specific_cfg.get("reference_genome", "unknown"))
    adata.uns.setdefault("original_gencode_version", dataset_specific_cfg.get("gencode_version", "unknown"))

    # Validation report will be added in Stage 2 by standardize_metadata.py script
    # Ontology mappings also added in Stage 2
    adata.uns["stage1_processing_date"] = pd.Timestamp.now().strftime("%Y-%m-%d")

    return adata



def save_anndata(adata, file_path):
    """Save AnnData object to file with validation and improved serialization."""
    try:
        logger.info(f"Preparing to save AnnData to {file_path}")

        if adata.n_vars == 0:
            logger.error("Cannot save AnnData: var DataFrame is empty!")
            return False
        if 'gene_id' not in adata.var.columns and len(adata.var_names) > 0:
             logger.warning("Reconstructing potentially missing 'gene_id' column in var from index.")
             adata.var['gene_id'] = adata.var_names.astype(str)
        elif adata.var.shape[1] == 0 and len(adata.var_names) > 0: # No columns but index exists
             logger.error("var DataFrame has no columns! Attempting reconstruction with 'gene_id'.")
             adata.var['gene_id'] = adata.var_names.astype(str)

        if 'gene_id' in adata.var.columns: # Ensure string type
            adata.var['gene_id'] = adata.var['gene_id'].astype(str)

        if adata.obs.index.name is not None and adata.obs.index.name in adata.obs.columns:
            new_name = f"original_obs_index_{adata.obs.index.name}"
            logger.warning(f"Fixing index/column name conflict in obs: Renaming column '{adata.obs.index.name}' to '{new_name}'")
            adata.obs = adata.obs.rename(columns={adata.obs.index.name: new_name})

        if adata.var.index.name is not None and adata.var.index.name in adata.var.columns:
            new_name = f"original_var_index_{adata.var.index.name}"
            logger.warning(f"Fixing index/column name conflict in var: Renaming column '{adata.var.index.name}' to '{new_name}'")
            adata.var = adata.var.rename(columns={adata.var.index.name: new_name})

        adata.obs.index = adata.obs.index.astype(str)
        adata.var.index = adata.var.index.astype(str)

        logger.info("Converting object columns in obs/var to string before saving...")
        for df_attr_name in ['obs', 'var']:
            df = getattr(adata, df_attr_name)
            # Create a new DataFrame for modifications to avoid SettingWithCopyWarning on slices
            new_df_data = {}
            object_cols_converted = False
            for col in df.columns:
                if df[col].dtype == 'object':
                    # Safely convert dicts/lists to JSON strings, then everything to string
                    def safe_stringify(x):
                        if isinstance(x, (dict, list)):
                            try:
                                return json.dumps(x, default=str)
                            except TypeError:
                                return str(x) # Fallback for unjsonable complex objects
                        return str(x) if pd.notna(x) else '' # Convert others to string, NaN to empty string

                    new_df_data[col] = df[col].apply(safe_stringify).astype(str)
                    logger.debug(f"Converted {df_attr_name} object column '{col}' to string.")
                    object_cols_converted = True
                else:
                    new_df_data[col] = df[col] # Keep non-object columns as is

            if object_cols_converted:
                setattr(adata, df_attr_name, pd.DataFrame(new_df_data, index=df.index))
        logger.info("Object column conversion complete.")

        logger.info("Ensuring adata.uns is serializable...")
        # Create a new dictionary for the serialized .uns to avoid modifying during iteration
        serializable_uns_dict = {}
        for key, value in adata.uns.items():
            try:
                serializable_uns_dict[key] = ensure_serializable(value)
            except Exception as e_uns_item:
                logger.error(f"Could not serialize adata.uns['{key}'] (type: {type(value)}). Error: {e_uns_item}. Storing as string representation.")
                serializable_uns_dict[key] = str(value)
        adata.uns.clear()
        adata.uns.update(serializable_uns_dict)
        logger.info("adata.uns converted to serializable format.")

        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        logger.info(f"Writing AnnData to {file_path}...")
        adata.write_h5ad(file_path)
        logger.info(f"Successfully wrote AnnData to {file_path}")

        logger.info(f"Verifying saved file: {file_path}")
        test_load = ad.read_h5ad(file_path)
        logger.info(
            f"Verification successful: Loaded AnnData with shape={test_load.shape}, "
            f"var columns ({len(test_load.var.columns)}): {list(test_load.var.columns)}, "
            f"obs columns ({len(test_load.obs.columns)}): {list(test_load.obs.columns)}"
        )
        if test_load.n_vars > 0 and test_load.var.shape[1] == 0 :
             logger.error("Verification FAILED: Loaded AnnData has empty var DataFrame columns despite having genes!")
             return False
        if test_load.n_obs > 0 and test_load.obs.shape[1] == 0:
             logger.error("Verification FAILED: Loaded AnnData has empty obs DataFrame columns despite having samples!")
             return False
        return True

    except Exception as e:
        logger.error(f"Unhandled error in save_anndata for {file_path}: {e}", exc_info=True)
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
    """Read an ENCODE TPM file - uses identify_columns."""
    try:
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        gene_id_col, expr_val_col = identify_columns(df) # expr_val_col is the identified expression column

        if gene_id_col is not None and expr_val_col is not None:
            df[gene_id_col] = df[gene_id_col].astype(str)
            # Keep only gene ID and the identified expression column
            result_df = df[[gene_id_col, expr_val_col]].set_index(gene_id_col)
            # RENAME the identified expression column to a STANDARD intermediate name "expression_value"
            result_df = result_df.rename(columns={expr_val_col: "expression_value"})
            logger.debug(f"Read {file_path}, gene_col='{gene_id_col}', expr_col='{expr_val_col}', renamed to 'expression_value'. Shape: {result_df.shape}")
            return result_df
        else:
            logger.warning(f"Could not identify gene_id and/or expression value columns in {file_path}. Gene: {gene_id_col}, Expr: {expr_val_col}")
            logger.warning(f"Columns present: {df.columns.tolist()}")
            return None
    except Exception as e:
        logger.warning(f"Error reading TPM file {file_path}: {e}")
        return None


def extract_encode_metadata(file_path, entex_full_metadata=None, encode_specific_cfg=None):
    """
    Extract metadata for an ENCODE or ENTEx file.
    Prioritizes specific metadata from entex_full_metadata if it's an ENTEx sample.
    Uses encode_specific_cfg for cell line info and general ENCODE defaults.
    """
    metadata = {}
    file_name = os.path.basename(file_path)
    file_id = file_name.split(".")[0]
    metadata["original_sample_id"] = file_id

    encode_cell_line_info_from_json = (encode_specific_cfg or {}).get("dataset_info", {}).get("cell_lines", {})

    is_entex_sample = "entex" in str(file_path).lower() # Check path first
    if entex_full_metadata and entex_full_metadata.get("sample_lookup") and file_id in entex_full_metadata["sample_lookup"]:
        is_entex_sample = True # Confirm if found in ENTEx lookup

    if is_entex_sample:
        metadata["data_source"] = "ENTEx"
        sample_info = (entex_full_metadata.get("sample_lookup", {})).get(file_id, {})
        if sample_info:
            logger.debug(f"Using specific ENTEx metadata for {file_id}")
            donor_id = sample_info.get("donor_id", sample_info.get("donor", sample_info.get("subject_id")))
            metadata["tissue"] = sample_info.get("tissue", "")
            metadata["subject_id"] = donor_id
            metadata["assay_type"] = sample_info.get("assay_type", sample_info.get("library_preparation", "unknown"))
            metadata["genome_annotation_source"] = sample_info.get("genome_annotation", "unknown")
            metadata["genome_assembly_source"] = sample_info.get("genome_assembly", "unknown")
            if donor_id and donor_id in entex_full_metadata.get("donor_map", {}):
                donor_info_map = entex_full_metadata["donor_map"][donor_id]
                metadata["sex"] = donor_info_map.get("sex", "")
                metadata["age"] = str(donor_info_map.get("age", ""))
        else:
            logger.warning(f"File {file_id} appears to be ENTEx by path, but no specific metadata in lookup. Using defaults.")
    else: # ENCODE Cell Line
        metadata["data_source"] = "ENCODE"
        path_parts = Path(file_path).parts
        cell_line_name_from_path = None
        for part in reversed(path_parts):
            for cl_key in encode_cell_line_info_from_json.keys():
                if cl_key.lower() in part.lower():
                    cell_line_name_from_path = cl_key
                    break
            if cell_line_name_from_path: break

        if cell_line_name_from_path:
            metadata["cell_line"] = cell_line_name_from_path
            # Set cell_type to the cell line name for Cell Ontology mapping
            metadata["cell_type"] = cell_line_name_from_path
            cell_specific_info = encode_cell_line_info_from_json.get(cell_line_name_from_path, {})
            # Update metadata with specific cell line info from JSON if available
            for key_to_update in ['tissue', 'sex', 'age', 'organism', 'subject_id', 'disease', 'ethnicity']:
                if key_to_update in cell_specific_info:
                    metadata[key_to_update] = cell_specific_info[key_to_update]

            relevant_path_part = next((p for p in path_parts if cell_line_name_from_path.lower() in p.lower()), "")
            if "_polya_plus" in relevant_path_part.lower() or "_polya+" in relevant_path_part.lower():
                metadata["extraction_method"] = "polyA_plus"
            elif "_polya_minus" in relevant_path_part.lower() or "_polya-" in relevant_path_part.lower():
                metadata["extraction_method"] = "polyA_minus"
            elif "_total" in relevant_path_part.lower():
                metadata["extraction_method"] = "total"
        else:
            logger.warning(f"Could not determine cell line for ENCODE sample {file_id} from path. Defaults will be used.")

    # Apply general ENCODE defaults from encode_specific_cfg if not already set
    if encode_specific_cfg and "obs_columns" in encode_specific_cfg:
        for key, default_val in encode_specific_cfg["obs_columns"].items():
            metadata.setdefault(key, default_val)

    # Final fallback for core fields if still missing
    for fld in ['tissue', 'subject_id', 'sex', 'age', 'data_type', 'expression_unit']:
        metadata.setdefault(fld, "unknown" if fld != 'age' else "")
    metadata['age'] = str(metadata['age']) # Ensure age is string

    return metadata

def process_encode_data(
    input_dir, entex_dir=None, output_file=None, entex_metadata_file=None, metadata_dir=None, gencode_map=None
): # Added gencode_map parameter
    logger.info(f"Processing ENCODE & ENTEx data. Cell line dir: {input_dir}, ENTEx dir: {entex_dir}")
    dataset_name_for_log = "ENCODE_ENTEx_Combined"

    encode_specific_cfg = load_dataset_specific_metadata(metadata_dir, "encode") or {}
    entex_full_metadata = load_entex_metadata(entex_metadata_file) if (entex_metadata_file or entex_dir) else {}

    all_tsv_files_paths = []
    if input_dir and os.path.exists(input_dir):
        all_tsv_files_paths.extend(glob.glob(os.path.join(input_dir, "**", "*.tsv"), recursive=True))
    if entex_dir and os.path.exists(entex_dir) and entex_dir != input_dir :
        all_tsv_files_paths.extend(glob.glob(os.path.join(entex_dir, "**", "*.tsv"), recursive=True))

    all_tsv_files_paths = sorted(list(set(all_tsv_files_paths)))

    final_tsv_files_to_process = []
    for file_path_str in all_tsv_files_paths:
        file_id = Path(file_path_str).name.split(".")[0]
        is_potentially_entex = "entex" in file_path_str.lower() or (entex_full_metadata.get("sample_lookup") and file_id in entex_full_metadata["sample_lookup"])

        if is_potentially_entex:
            sample_info = (entex_full_metadata.get("sample_lookup", {})).get(file_id)
            if sample_info:
                quant_method = sample_info.get("quantification_method", "").lower()
                genome_annotation = sample_info.get("genome_annotation", "").upper()
                is_preferred_quant = "gene" in quant_method and "transcript" not in quant_method
                is_preferred_annotation = genome_annotation in ["V24", "V29", ""] or not genome_annotation
                if is_preferred_quant and is_preferred_annotation:
                    final_tsv_files_to_process.append(file_path_str)
                else:
                    logger.debug(f"Skipping ENTEx file {file_id} (from {file_path_str}) due to non-preferred quant/annotation: {quant_method}, {genome_annotation}")
            else:
                logger.warning(f"File {file_id} (from {file_path_str}) seems like ENTEx by path but no specific metadata. Including.")
                final_tsv_files_to_process.append(file_path_str)
        else:
            final_tsv_files_to_process.append(file_path_str)

    logger.info(f"Total files after initial scan: {len(all_tsv_files_paths)}. After filtering for ENCODE/ENTEx: {len(final_tsv_files_to_process)}")

    if not final_tsv_files_to_process:
        logger.error(f"No TSV files to process for ENCODE/ENTEx.")
        return None

    sample_dfs = {}
    sample_metadata_list = []
    all_gene_ids = set()

    for file_path_str in final_tsv_files_to_process:
        tpm_df = read_encode_tpm_file(file_path_str)
        if tpm_df is None:
            continue
        file_name = os.path.basename(file_path_str)
        sample_id = file_name.split(".")[0]

        if "expression_value" in tpm_df.columns:
            tpm_df = tpm_df.rename(columns={"expression_value": sample_id})
        else: # Should ideally not happen if read_encode_tpm_file worked
            logger.error(f"For {sample_id}, 'expression_value' column was expected but not found after read_encode_tpm_file. Columns: {tpm_df.columns.tolist()}")
            continue

        sample_dfs[sample_id] = tpm_df
        current_metadata = extract_encode_metadata(file_path_str, entex_full_metadata, encode_specific_cfg)
        current_metadata["sample_id_from_file"] = sample_id
        sample_metadata_list.append(current_metadata)
        all_gene_ids.update(tpm_df.index)

    if not sample_dfs:
        logger.error(f"No valid TPM data found for {dataset_name_for_log} (ENCODE/ENTEx)")
        return None

    num_processed_encode_entex_files = len(sample_dfs)
    logger.info(f"Successfully read data from {num_processed_encode_entex_files} ENCODE/ENTEx files.")

    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs for {dataset_name_for_log}")
    gene_id_mapping = {gid: standardize_ensembl_id(gid) for gid in all_gene_ids}

    valid_mapped_values_encode = [str(gid) for gid in gene_id_mapping.values() if gid is not None]
    if not valid_mapped_values_encode:
        logger.error(f"{dataset_name_for_log}: No gene IDs could be meaningfully standardized. Aborting processing for this dataset.")
        return None
    unique_std_ids = sorted(list(set(valid_mapped_values_encode)))

    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs for {dataset_name_for_log}")

    std_id_to_idx = {gid: i for i, gid in enumerate(unique_std_ids)}

    obs_df_encode = pd.DataFrame(sample_metadata_list)
    if 'sample_id_from_file' not in obs_df_encode.columns or obs_df_encode.empty:
        logger.error(f"{dataset_name_for_log}: 'sample_id_from_file' key missing or obs_df is empty.")
        return None

    obs_df_encode = obs_df_encode.set_index('sample_id_from_file')

    final_sample_ids_encode = obs_df_encode.index.tolist()
    if obs_df_encode.index.has_duplicates:
        logger.warning(f"{dataset_name_for_log}: Correcting duplicate sample IDs in obs_df index. Original count: {len(obs_df_encode)}")
        obs_df_encode = obs_df_encode[~obs_df_encode.index.duplicated(keep='first')]
        final_sample_ids_encode = obs_df_encode.index.tolist()

    sample_dfs = {sid: df_val for sid, df_val in sample_dfs.items() if sid in final_sample_ids_encode}

    num_genes_encode = len(unique_std_ids)
    num_samples_encode = len(final_sample_ids_encode)

    if num_samples_encode == 0:
        logger.error(f"{dataset_name_for_log}: No samples remaining after metadata processing and filtering. Aborting.")
        return None
    if num_genes_encode == 0:
        logger.error(f"{dataset_name_for_log}: No genes remaining after standardization. Aborting.")
        return None

    logger.info(f"Creating unified expression matrix for {dataset_name_for_log} ({num_genes_encode} genes x {num_samples_encode} samples)")

    expr_matrix = np.full((num_genes_encode, num_samples_encode), np.nan, dtype=np.float32)
    sample_id_to_col_idx = {sid: i for i, sid in enumerate(final_sample_ids_encode)}

    logger.info("DEBUG: Checking sample_dfs structure before ENCODE/ENTEx aggregation loop...")
    for temp_sample_id_debug, temp_df_debug in sample_dfs.items():
        if temp_df_debug is None:
            logger.info(f"DEBUG: sample_dfs['{temp_sample_id_debug}'] is None")
            continue
        # Ensure the column name is indeed the sample_id (which it should be after rename)
        if temp_sample_id_debug not in temp_df_debug.columns:
            logger.warning(f"DEBUG WARNING: Column '{temp_sample_id_debug}' NOT FOUND in sample_dfs['{temp_sample_id_debug}']. Columns are: {temp_df_debug.columns.tolist()}")

    std_to_originals_map_encode = {}
    for original_gid, mapped_std_id_val in gene_id_mapping.items():
        if mapped_std_id_val is not None:
            str_mapped_std_id_val = str(mapped_std_id_val)
            std_to_originals_map_encode.setdefault(str_mapped_std_id_val, []).append(str(original_gid))

    for std_id, original_ids_for_this_std in std_to_originals_map_encode.items():
        if std_id not in std_id_to_idx:
            logger.debug(f"Standardized ID '{std_id}' not in final index for {dataset_name_for_log}. Skipping.")
            continue
        gene_idx = std_id_to_idx[std_id]
        for sample_id, col_idx in sample_id_to_col_idx.items():
            sample_df = sample_dfs.get(sample_id)
            if sample_df is None: continue
            if sample_id not in sample_df.columns: # This check is crucial
                logger.error(f"CRITICAL DEBUG in ENCODE: Column '{sample_id}' missing in its DataFrame. Columns: {sample_df.columns.tolist()}")
                continue

            relevant_orig_ids_for_this_std_id = [orig_id for orig_id in original_ids_for_this_std if orig_id in sample_df.index]
            if relevant_orig_ids_for_this_std_id:
                values = sample_df.loc[relevant_orig_ids_for_this_std_id, sample_id]
                if not values.empty and not values.isna().all():
                    expr_matrix[gene_idx, col_idx] = values.max()

    expr_matrix = np.nan_to_num(expr_matrix, nan=0.0)
    unified_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=final_sample_ids_encode)

    var_df_encode = pd.DataFrame(index=unique_std_ids)
    var_df_encode["gene_id"] = var_df_encode.index.astype(str)
    var_df_encode["original_ids"] = var_df_encode.index.map(lambda x: ";".join(std_to_originals_map_encode.get(x, [])))

    # gencode_map = load_gencode_mapping() # Removed: use passed gencode_map
    var_df_encode = add_gencode_annotations(var_df_encode, gencode_map)

    mappings = load_mappings()
    obs_df_encode = standardize_metadata(obs_df_encode, "ENCODE", mappings)
    if obs_df_encode is None: # Check if standardize_metadata failed
        logger.error("ENCODE: standardize_metadata returned None. Aborting ENCODE/ENTEx processing.")
        return None


    # ---- ADD DEBUG LOGS BEFORE CREATING ANNDATA FOR ENCODE ---
    logger.info(f"ENCODE: Type of unified_expr_df.T: {type(unified_expr_df.T)}, Shape: {unified_expr_df.T.shape if hasattr(unified_expr_df, 'T') else 'N/A'}")
    logger.info(f"ENCODE: Type of obs_df_encode: {type(obs_df_encode)}, Shape: {obs_df_encode.shape if hasattr(obs_df_encode, 'shape') else 'N/A'}")
    logger.info(f"ENCODE: Type of var_df_encode: {type(var_df_encode)}, Shape: {var_df_encode.shape if hasattr(var_df_encode, 'shape') else 'N/A'}")

    if unified_expr_df is None: logger.error("CRITICAL ENCODE: unified_expr_df is None!")
    if obs_df_encode is None: logger.error("CRITICAL ENCODE: obs_df_encode is None!")
    if var_df_encode is None: logger.error("CRITICAL ENCODE: var_df_encode is None!")
    # ---- END DEBUG LOGS ---

    adata = create_standard_anndata(unified_expr_df.T, obs_df_encode, var_df_encode, encode_specific_cfg, dataset_name_for_log)

    if encode_specific_cfg:
        adata = apply_dataset_specific_metadata(adata, encode_specific_cfg)

    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed {dataset_name_for_log} data: {adata.n_obs} samples and {adata.n_vars} genes")
    else:
        logger.error(f"Failed to save AnnData for {dataset_name_for_log}")
        return None
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
        sample_attrs = pd.read_csv(GTEx_SAMPLE_ATTRS, sep="\t", low_memory=False)

        # Load subject attributes
        logger.info(f"Loading subject attributes from {GTEx_SUBJECT_ATTRS}")
        subject_attrs = pd.read_csv(GTEx_SUBJECT_ATTRS, sep="\t", low_memory=False)

        # Extract SUBJID from SAMPID
        sample_attrs["SUBJID"] = sample_attrs["SAMPID"].str.extract(r"(GTEX-[A-Z0-9]+)")

        # Merge sample and subject attributes
        merged_attrs = pd.merge(sample_attrs, subject_attrs, on="SUBJID", how="left")

        # Standardize column names based on a mapping
        column_mapping = {
            "SAMPID": "sample_id", # Will become index
            "SUBJID": "subject_id",
            "SMTSD": "tissue",  # Detailed tissue site
            "SMTS": "broad_tissue",  # Broad tissue type
            "SEX": "sex_code",
            "AGE": "age_bracket",
            "RACE": "race_code_source",
            "ETHNCTY": "ethnicity_code_source",
        }

        cols_to_rename = {k: v for k, v in column_mapping.items() if k in merged_attrs.columns}
        merged_attrs.rename(columns=cols_to_rename, inplace=True)

        gtex_cfg = load_dataset_specific_metadata(DEFAULT_METADATA_DIR, "gtex") or {}

        if 'race_code_source' in merged_attrs.columns and gtex_cfg.get('source_race_to_standard_label_map'):
            race_map = {str(k): v for k,v in gtex_cfg['source_race_to_standard_label_map'].items()}
            merged_attrs['race'] = merged_attrs['race_code_source'].astype(str).map(race_map).fillna('unknown or not reported')
            logger.info("Mapped GTEx 'race_code_source' to 'race' column.")
        else:
            merged_attrs['race'] = 'unknown or not reported'
            logger.warning(f"GTEx: Could not map race codes. Source column 'race_code_source' or mapping in JSON missing.")

        if 'ethnicity_code_source' in merged_attrs.columns and gtex_cfg.get('source_ethnicity_to_standard_label_map'):
            eth_map = {str(k): v for k,v in gtex_cfg['source_ethnicity_to_standard_label_map'].items()}
            merged_attrs['is_hispanic_or_latino'] = merged_attrs['ethnicity_code_source'].astype(str).map(eth_map).fillna('unknown or not reported')
            logger.info("Mapped GTEx 'ethnicity_code_source' to 'is_hispanic_or_latino' column.")
        else:
            merged_attrs['is_hispanic_or_latino'] = 'unknown or not reported'
            logger.warning(f"GTEx: Could not map ethnicity codes. Source column 'ethnicity_code_source' or mapping in JSON missing.")

        if 'sex_code' in merged_attrs.columns:
            sex_map = {'1': 'male', '2': 'female'}
            merged_attrs['sex'] = merged_attrs['sex_code'].astype(str).map(sex_map).fillna('unknown')
        elif 'sex' not in merged_attrs.columns:
            merged_attrs['sex'] = 'unknown'

        if 'age_bracket' in merged_attrs.columns:
            merged_attrs['age'] = merged_attrs['age_bracket'].astype(str)
        elif 'age' not in merged_attrs.columns:
            merged_attrs['age'] = ''

        if "tissue" not in merged_attrs.columns and "broad_tissue" in merged_attrs.columns:
            merged_attrs["tissue"] = merged_attrs["broad_tissue"]
            logger.info("Using 'broad_tissue' as 'tissue' column for GTEx.")
        elif "tissue" not in merged_attrs.columns:
             merged_attrs["tissue"] = "unknown"

        if "sample_id" in merged_attrs.columns:
            merged_attrs.set_index("sample_id", inplace=True)
        else:
            logger.error("GTEx metadata merging failed to produce 'sample_id' column from SAMPID.")
            return pd.DataFrame()

        logger.info(f"Loaded and initially processed metadata for {len(merged_attrs)} GTEx samples")
        return merged_attrs

    except Exception as e:
        logger.error(f"Error loading GTEx metadata: {e}", exc_info=True)
        return pd.DataFrame()


def process_gtex_data(input_file, output_file, metadata_dir=None, gencode_map=None): # Added gencode_map parameter
    logger.info(f"Processing GTEx data from {input_file} (Memory Optimized V2)")
    dataset_name_for_log = "GTEx"

    metadata_df = load_gtex_metadata()
    if metadata_df.empty:
        logger.error("Failed to load GTEx metadata. Cannot proceed.")
        return None
    logger.info(f"Loaded metadata for {len(metadata_df)} potential GTEx samples.")

    sample_ids_in_gct = []
    n_genes_in_gct = 0
    try:
        with gzip.open(input_file, "rt") as f:
            next(f)
            dims = next(f).strip().split()
            n_genes_in_gct = int(dims[0])
            header_line = next(f).strip().split("\t")
            sample_ids_in_gct = header_line[2:]
        logger.info(f"GCT header indicates {n_genes_in_gct} genes and {len(sample_ids_in_gct)} samples.")
    except Exception as e:
        logger.error(f"Error reading GCT header: {e}")
        return None

    common_samples_set = set(metadata_df.index).intersection(set(sample_ids_in_gct))
    if not common_samples_set:
        logger.error("No common samples found between metadata and GCT file.")
        return None

    final_sample_ids = sorted(list(common_samples_set))
    n_final_samples = len(final_sample_ids)
    logger.info(f"Processing {n_final_samples} common samples.")
    obs_df = metadata_df.loc[final_sample_ids].copy()

    obs_df['expression_unit'] = 'TPM'
    obs_df['data_type'] = 'RNA-seq'

    logger.info("Standardizing observation metadata for GTEx...")
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, dataset_name_for_log, mappings)
    if obs_df is None: # Check if standardize_metadata failed
        logger.error("GTEx: standardize_metadata returned None. Aborting GTEx processing.")
        return None
    obs_df['sample_id'] = obs_df.index.astype(str)

    gct_sample_indices = {sid: i for i, sid in enumerate(sample_ids_in_gct)}
    sample_indices_to_keep = [gct_sample_indices[sid] for sid in final_sample_ids]

    logger.info("First pass: Identifying unique standardized genes...")
    std_gene_id_map = {}
    all_std_gene_ids = set()
    original_id_to_std_map = {}

    line_count = 0
    with gzip.open(input_file, "rt") as f:
        next(f); next(f); next(f)
        for line in f:
            line_count += 1
            try:
                original_gene_id = line.strip().split("\t", 1)[0]
                std_gene_id = standardize_ensembl_id(original_gene_id)
                if std_gene_id:
                    std_gene_id_map[original_gene_id] = std_gene_id
                    all_std_gene_ids.add(std_gene_id)
                    if std_gene_id not in original_id_to_std_map:
                        original_id_to_std_map[std_gene_id] = []
                    original_id_to_std_map[std_gene_id].append(original_gene_id)
            except IndexError:
                logger.warning(f"Skipping malformed line {line_count}")
                continue
            if line_count % 10000 == 0:
                logger.info(f"First pass processed {line_count}/{n_genes_in_gct} genes...")

    final_std_gene_ids = sorted(list(all_std_gene_ids))
    if not final_std_gene_ids:
        logger.error("GTEx: No genes could be standardized. Aborting GTEx processing.")
        return None
    n_final_genes = len(final_std_gene_ids)
    std_gene_id_to_idx = {gid: i for i, gid in enumerate(final_std_gene_ids)}
    logger.info(f"Identified {n_final_genes} unique standardized genes.")

    logger.info("Second pass: Populating expression matrix...")
    final_expr_matrix = np.full((n_final_genes, n_final_samples), -1.0, dtype=np.float32)

    line_count = 0
    with gzip.open(input_file, "rt") as f:
        next(f); next(f); next(f)
        for line in f:
            line_count += 1
            try:
                fields = line.strip().split("\t")
                original_gene_id = fields[0]
                std_gene_id = std_gene_id_map.get(original_gene_id)

                if not std_gene_id: continue

                std_idx = std_gene_id_to_idx.get(std_gene_id)
                if std_idx is None: continue

                expression_values_all = fields[2:]
                expression_values_filtered = np.array([float(expression_values_all[i]) for i in sample_indices_to_keep], dtype=np.float32)

                current_row = final_expr_matrix[std_idx, :]
                update_mask = (expression_values_filtered > current_row) | (current_row == -1.0)
                final_expr_matrix[std_idx, update_mask] = expression_values_filtered[update_mask]

            except Exception as line_e:
                logger.warning(f"Skipping line {line_count} during second pass due to error: {line_e}")
                continue
            if line_count % 5000 == 0:
                logger.info(f"Second pass processed {line_count}/{n_genes_in_gct} genes...")

    final_expr_matrix[final_expr_matrix == -1.0] = 0.0
    logger.info("Finished populating expression matrix.")

    var_df = pd.DataFrame(index=final_std_gene_ids)
    var_df['gene_id'] = var_df.index.astype(str)
    var_df['original_ids'] = [";".join(original_id_to_std_map.get(gid, [])) for gid in final_std_gene_ids]
    logger.info("Adding GENCODE annotations...")
    # gencode_mapping = load_gencode_mapping() # Removed: use passed gencode_map
    var_df = add_gencode_annotations(var_df, gencode_map)

    gtex_specific_cfg = load_dataset_specific_metadata(metadata_dir, "gtex") or {}
    adata = create_standard_anndata(pd.DataFrame(final_expr_matrix.T, index=final_sample_ids, columns=final_std_gene_ids),
                                    obs_df, var_df, gtex_specific_cfg, dataset_name_for_log)
    logger.info(f"Created AnnData object with shape {adata.shape}")

    if gtex_specific_cfg:
        adata = apply_dataset_specific_metadata(adata, gtex_specific_cfg)

    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed GTEx data with {adata.n_obs} samples and {adata.n_vars} genes (Memory Optimized V2)")
    else:
        logger.error("Failed to save memory-optimized GTEx AnnData V2.")
        return None

    return adata

# ====================================================================================
# MAGE-specific Processing
# ====================================================================================

def process_mage_dir(file_path, donor_id, tissue, sample_dfs, sample_metadata_list, all_gene_ids):
    try:
        file_name = os.path.basename(file_path)
        sample_id = f"{donor_id}_{os.path.splitext(file_name)[0]}"
        logger.debug(f"Processing MAGE file: {file_path} as sample_id: {sample_id}")

        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            logger.warning(f"File does not exist or is empty: {file_path}")
            return False

        with open(file_path, "r") as f:
            first_line = f.readline().strip()
        delimiter = "\t" if "\t" in first_line else ","

        df = pd.read_csv(file_path, delimiter=delimiter, low_memory=False)
        if df.empty:
            logger.warning(f"Empty DataFrame after reading MAGE file: {file_path}")
            return False

        gene_id_col_name, expr_col_name = identify_columns(df)
        expr_type = "unknown_expression"
        if expr_col_name and "TPM" in expr_col_name.upper():
            expr_type = "TPM"

        if not (gene_id_col_name and expr_col_name and gene_id_col_name in df.columns and expr_col_name in df.columns):
            logger.warning(f"Could not reliably identify gene ID ('{gene_id_col_name}') and/or expression ('{expr_col_name}') columns in {file_path}. Columns: {df.columns.tolist()}")
            return False

        gene_id_series = df[gene_id_col_name].astype(str)
        expr_series = pd.to_numeric(df[expr_col_name], errors='coerce')

        valid_expr_mask = expr_series.notna()
        if not valid_expr_mask.any():
            logger.warning(f"No valid numeric expression values in column '{expr_col_name}' for {file_path}.")
            return False

        gene_id_series = gene_id_series[valid_expr_mask]
        expr_series = expr_series[valid_expr_mask]

        if gene_id_series.empty:
            logger.warning(f"Gene ID series became empty after filtering for valid expression in {file_path}.")
            return False

        simple_df = pd.DataFrame({sample_id: expr_series.values}, index=gene_id_series.values)

        if simple_df.index.has_duplicates:
            logger.warning(f"Duplicate gene IDs found in {file_path}. Aggregating by max.")
            simple_df = simple_df.groupby(simple_df.index).max()

        sample_dfs[sample_id] = simple_df
        all_gene_ids.update(simple_df.index)

        metadata = {
            "sample_id_from_file": sample_id,
            "donor_id": donor_id, "subject_id": donor_id, "tissue": tissue,
            "dataset": "MAGE", "data_type": "RNA-seq", "expression_unit": expr_type,
            "is_gencode": "gencode" in file_name.lower(),
            # Set cell_type for Cell Ontology mapping (MAGE uses lymphoblastoid cell lines)
            "cell_type": "lymphoblastoid cell line",
            # Add placeholders for sex and ethnicity, to be filled from 1000G data
            "sex": "unknown",
            "race": "unknown or not reported",
            "is_hispanic_or_latino": "unknown or not reported",
            "population_code_1000g": "unknown"
        }
        sample_metadata_list.append(metadata)
        logger.debug(f"Successfully processed {file_path} with {len(simple_df)} genes for sample {sample_id}")
        return True

    except Exception as e:
        logger.error(f"Error processing MAGE file {file_path}: {e}", exc_info=True)
        return False


def process_mage_data(input_dir, output_file, metadata_dir=None, gencode_map=None): # Added gencode_map
    logger.info(f"Processing MAGE data from {input_dir}")
    dataset_name_for_log = "MAGE"
    mage_specific_cfg = load_dataset_specific_metadata(metadata_dir, "mage") or {}

    # --- Load 1000 Genomes Metadata ---
    genomes_metadata_file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/project_gene_regulation/data/MAGE/WGS/20130606_g1k_3202_samples_ped_population.txt'
    genomes_metadata_map = {}
    if os.path.exists(genomes_metadata_file_path):
        try:
            logger.info(f"MAGE: Loading 1000 Genomes metadata from {genomes_metadata_file_path}")
            ped_df = pd.read_csv(genomes_metadata_file_path, sep='\s+', dtype=str)
            ped_df.columns = ped_df.columns.str.strip()

            required_ped_cols_for_mapping = ['SampleID', 'Sex', 'Population']
            missing_cols = [col for col in required_ped_cols_for_mapping if col not in ped_df.columns]

            if missing_cols:
                logger.error(f"MAGE: 1000G PED file missing required columns for mapping after stripping: {missing_cols}. Columns found: {ped_df.columns.tolist()}")
            else:
                for _, row in ped_df.iterrows():
                    sample_id_1000g = str(row['SampleID']).strip().upper()
                    sex_code = str(row['Sex']).strip()
                    population_code = str(row['Population']).strip()
                    genomes_metadata_map[sample_id_1000g] = {
                        'sex_code_1000g': sex_code,
                        'population_code_1000g': population_code
                    }
                logger.info(f"MAGE: Loaded metadata for {len(genomes_metadata_map)} samples from 1000 Genomes PED file.")
                # ---- MAGE DEBUGGING ----
                if genomes_metadata_map:
                    logger.info(f"MAGE DEBUG: First 5 keys in genomes_metadata_map: {list(genomes_metadata_map.keys())[:5]}")
                # ---- END MAGE DEBUGGING ----

        except Exception as e:
            logger.error(f"MAGE: Error loading or processing 1000 Genomes metadata file {genomes_metadata_file_path}: {e}", exc_info=True)
    else:
        logger.warning(f"MAGE: 1000 Genomes metadata file not found at {genomes_metadata_file_path}. Sex and ethnicity will be 'unknown'.")




    sample_dfs = {}
    sample_metadata_list = []
    all_gene_ids = set()
    processed_samples_count = 0
    donor_dirs = glob.glob(os.path.join(input_dir, "NA*")) + glob.glob(os.path.join(input_dir, "HG*"))

    if not donor_dirs:
        logger.error(f"No donor directories (NA* or HG*) found in MAGE input directory: {input_dir}")
        return None

    for donor_dir_path in donor_dirs:
        donor_id_original_case = os.path.basename(donor_dir_path)
        donor_id_lookup_key = donor_id_original_case.upper() # Normalize for lookup
        # ---- MAGE DEBUGGING ----
        if processed_samples_count < 5: # Log for the first few samples
            logger.info(f"MAGE DEBUG: Processing donor_id_original_case: '{donor_id_original_case}', lookup_key: '{donor_id_lookup_key}'")
        # ---- END MAGE DEBUGGING ----


        lymphoblast_dir = os.path.join(donor_dir_path, "lymphoblast")
        if not os.path.isdir(lymphoblast_dir):
            logger.warning(f"No 'lymphoblast' subdirectory for MAGE donor {donor_id_original_case}, skipping.")
            continue

        all_csv_files_in_dir = glob.glob(os.path.join(lymphoblast_dir, "*.csv"))
        preferred_file = None
        gencode_v24_pruned_files = [f for f in all_csv_files_in_dir if Path(f).name.endswith("_gencode_V24_pruned.csv")]

        if gencode_v24_pruned_files:
            preferred_file = gencode_v24_pruned_files[0]
            if len(gencode_v24_pruned_files) > 1:
                logger.warning(f"Multiple '*_gencode_V24_pruned.csv' files for MAGE donor {donor_id_original_case} in {lymphoblast_dir}. Using {Path(preferred_file).name}")
        else:
            simple_name_files = [f for f in all_csv_files_in_dir if Path(f).name == f"{donor_id_original_case}.csv"]
            if simple_name_files:
                preferred_file = simple_name_files[0]
                logger.info(f"Using '{Path(preferred_file).name}' for MAGE donor {donor_id_original_case} as no '*_gencode_V24_pruned.csv' was found.")
                if len(simple_name_files) > 1:
                     logger.warning(f"Multiple '{donor_id_original_case}.csv' files found (unexpected). Using the first one.")
            elif all_csv_files_in_dir:
                preferred_file = all_csv_files_in_dir[0]
                logger.warning(f"No '*_gencode_V24_pruned.csv' or '{donor_id_original_case}.csv' file for MAGE donor {donor_id_original_case}. Using first available CSV: {Path(preferred_file).name}.")

        if preferred_file:
            logger.info(f"Selected file for MAGE donor {donor_id_original_case}: {Path(preferred_file).name}")
            if process_mage_dir(preferred_file, donor_id_original_case, "lymphoblast", sample_dfs, sample_metadata_list, all_gene_ids):
                processed_samples_count +=1
                if sample_metadata_list:
                    entry_to_update = next((item for item in reversed(sample_metadata_list) if item.get('donor_id') == donor_id_original_case), None)
                    if entry_to_update:
                        sample_1000g_info = genomes_metadata_map.get(donor_id_lookup_key) # Use normalized key for lookup
                        if sample_1000g_info:
                            logger.info(f"MAGE SUCCESS: Found 1000G metadata for {donor_id_original_case} using key {donor_id_lookup_key}") # Add success log
                            sex_code_1000g = sample_1000g_info.get('sex_code_1000g')
                            pop_code_1000g = sample_1000g_info.get('population_code_1000g')
                            sex_map_1000g = {'1': 'male', '2': 'female'}
                            entry_to_update['sex'] = sex_map_1000g.get(sex_code_1000g, 'unknown')
                            pop_to_race_map = mage_specific_cfg.get('pop_to_race_map', {'ACB':'black or african american', 'ASW':'black or african american', 'ESN':'black or african american', 'GWD':'black or african american', 'LWK':'black or african american', 'MSL':'black or african american', 'YRI':'black or african american', 'CLM':'multiethnic', 'MXL':'multiethnic', 'PEL':'multiethnic', 'PUR':'multiethnic', 'CDX':'asian', 'CHB':'asian', 'CHS':'asian', 'JPT':'asian', 'KHV':'asian', 'CEU':'white', 'FIN':'white', 'GBR':'white', 'IBS':'white', 'TSI':'white', 'BEB':'asian', 'GIH':'asian', 'ITU':'asian', 'PJL':'asian', 'STU':'asian'})
                            pop_to_hispanic_map = mage_specific_cfg.get('pop_to_hispanic_map', {'CLM':'hispanic or latino', 'MXL':'hispanic or latino', 'PEL':'hispanic or latino', 'PUR':'hispanic or latino'})
                            default_race = mage_specific_cfg.get('default_race_for_unmapped_pop', 'unknown or not reported')
                            derived_race = pop_to_race_map.get(pop_code_1000g, default_race)
                            derived_is_hispanic = pop_to_hispanic_map.get(pop_code_1000g, 'not hispanic or latino')
                            if derived_race == default_race and derived_is_hispanic == 'hispanic or latino': derived_race = 'multiethnic'
                            entry_to_update['race'] = derived_race
                            entry_to_update['is_hispanic_or_latino'] = derived_is_hispanic
                            entry_to_update['population_code_1000g'] = pop_code_1000g
                        else:
                            logger.warning(f"MAGE: No 1000G metadata found for sample {donor_id_original_case} (lookup key: {donor_id_lookup_key}). Sex/Ethnicity remain defaults.")
                    else:
                        logger.error(f"MAGE: Could not find metadata entry for {donor_id_original_case} to enrich with 1000G data.")
        else:
            logger.warning(f"No suitable CSV file found for MAGE donor {donor_id_original_case} in {lymphoblast_dir}")

    if processed_samples_count == 0:
        logger.error(f"No MAGE samples were successfully processed from {input_dir}.")
        return None
    logger.info(f"Processed {processed_samples_count} MAGE samples.")

    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs for {dataset_name_for_log}")
    gene_id_mapping = {gid: standardize_ensembl_id(gid) for gid in all_gene_ids}

    valid_mapped_values_mage = [str(gid) for gid in gene_id_mapping.values() if gid is not None]
    if not valid_mapped_values_mage:
        logger.error(f"{dataset_name_for_log}: No gene IDs could be meaningfully standardized. Aborting.")
        return None
    unique_std_ids = sorted(list(set(valid_mapped_values_mage)))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs for {dataset_name_for_log}")

    std_id_to_idx = {gid: i for i, gid in enumerate(unique_std_ids)}

    obs_df_mage = pd.DataFrame(sample_metadata_list)
    if 'sample_id_from_file' not in obs_df_mage.columns or obs_df_mage.empty:
        logger.error(f"{dataset_name_for_log}: 'sample_id_from_file' key missing or obs_df is empty.")
        return None
    obs_df_mage = obs_df_mage.set_index('sample_id_from_file')
    obs_df_mage.index.name = "sample_id"

    final_sample_ids_mage = obs_df_mage.index.tolist()
    if obs_df_mage.index.has_duplicates:
        logger.warning(f"{dataset_name_for_log}: Correcting duplicate sample IDs in obs_df index.")
        obs_df_mage = obs_df_mage[~obs_df_mage.index.duplicated(keep='first')]
        final_sample_ids_mage = obs_df_mage.index.tolist()

    sample_dfs = {sid: df_val for sid, df_val in sample_dfs.items() if sid in final_sample_ids_mage}

    num_genes_mage = len(unique_std_ids)
    num_samples_mage = len(final_sample_ids_mage)

    if num_samples_mage == 0 or num_genes_mage == 0:
        logger.error(f"{dataset_name_for_log}: No samples ({num_samples_mage}) or genes ({num_genes_mage}) remaining. Aborting.")
        return None

    logger.info(f"Creating unified expression matrix for {dataset_name_for_log} ({num_genes_mage} genes x {num_samples_mage} samples) using optimized method.")

    # Optimized Matrix Creation for MAGE
    all_series = []
    for sample_id, s_df in sample_dfs.items():
        series = s_df[sample_id].rename(sample_id)
        all_series.append(series)

    if not all_series:
        logger.error(f"No data series collected for {dataset_name_for_log}. Aborting matrix creation.")
        unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_mage).astype(np.float32)

    else:
        logger.info(f"Concatenating {len(all_series)} sample series for {dataset_name_for_log}...")
        try:
            concat_df = pd.concat(all_series, axis=1, join='outer').astype(np.float32)
            logger.info(f"Concatenated DataFrame shape: {concat_df.shape}")

            std_to_originals_map_mage = {}
            for original_gid, mapped_std_id_val in gene_id_mapping.items():
                if mapped_std_id_val is not None:
                    std_to_originals_map_mage.setdefault(str(mapped_std_id_val), []).append(str(original_gid))

            original_to_std_series_map = {orig_id: std_id_val for std_id_val, orig_ids_list in std_to_originals_map_mage.items() for orig_id in orig_ids_list}
            mappable_original_ids = [idx for idx in concat_df.index if idx in original_to_std_series_map]

            if not mappable_original_ids:
                logger.warning(f"No mappable original gene IDs in concat_df for {dataset_name_for_log}. Resulting matrix might be mostly zeros or have unexpected shape.")
                unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_mage).astype(np.float32)
            else:
                filtered_concat_df = concat_df.loc[mappable_original_ids]
                std_ids_for_aggregation = filtered_concat_df.index.map(original_to_std_series_map)

                if std_ids_for_aggregation.isna().all():
                    logger.error(f"All standardized IDs for aggregation are NaN for {dataset_name_for_log}. Check mapping.")
                    unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_mage).astype(np.float32)
                else:
                    filtered_concat_df.index = std_ids_for_aggregation
                    logger.info(f"Grouping by standardized gene IDs and aggregating for {dataset_name_for_log}...")
                    unified_expr_df_grouped = filtered_concat_df.groupby(level=0).max()
                    logger.info(f"Aggregation complete. Shape: {unified_expr_df_grouped.shape}")
                    unified_expr_df = unified_expr_df_grouped.reindex(index=unique_std_ids, columns=final_sample_ids_mage, fill_value=0.0).astype(np.float32)
        except Exception as e_concat:
            logger.error(f"Error during optimized matrix creation for {dataset_name_for_log}: {e_concat}", exc_info=True)
            logger.warning("Falling back to previous matrix creation method for MAGE due to error.")
            # Fallback to original loop method if concat/groupby fails
            expr_matrix = np.full((num_genes_mage, num_samples_mage), np.nan, dtype=np.float32)
            sample_id_to_col_idx = {sid: i for i, sid in enumerate(final_sample_ids_mage)}
            std_to_originals_map_mage = {} # Re-create as it might not be in scope if optimized failed early
            for original_gid, mapped_std_id_val in gene_id_mapping.items():
                if mapped_std_id_val is not None:
                    std_to_originals_map_mage.setdefault(str(mapped_std_id_val), []).append(str(original_gid))

            for std_id, original_ids_for_this_std in std_to_originals_map_mage.items():
                if std_id not in std_id_to_idx:
                    logger.warning(f"Standardized ID '{std_id}' not in final index for MAGE (fallback). Skipping.")
                    continue
                gene_idx = std_id_to_idx[std_id]
                for sample_id, col_idx in sample_id_to_col_idx.items():
                    sample_df = sample_dfs.get(sample_id)
                    if sample_df is None: continue
                    if sample_id not in sample_df.columns:
                        logger.error(f"MAGE Aggregation Error (fallback): Column '{sample_id}' not found in sample_df. Columns: {sample_df.columns.tolist()}")
                        continue
                    relevant_orig_ids_in_sample = [orig_id for orig_id in original_ids_for_this_std if orig_id in sample_df.index]
                    if relevant_orig_ids_in_sample:
                        values = pd.to_numeric(sample_df.loc[relevant_orig_ids_in_sample, sample_id], errors='coerce')
                        if not values.empty and not values.isna().all():
                            expr_matrix[gene_idx, col_idx] = values.max()
            expr_matrix = np.nan_to_num(expr_matrix, nan=0.0)
            unified_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=final_sample_ids_mage)
    # End of Optimized Matrix Creation

    var_df_mage = pd.DataFrame(index=unique_std_ids)
    var_df_mage["gene_id"] = var_df_mage.index.astype(str)
    var_df_mage["original_ids"] = var_df_mage.index.map(lambda x: ";".join(std_to_originals_map_mage.get(x, [])))

    # gencode_map passed as parameter
    var_df_mage = add_gencode_annotations(var_df_mage, gencode_map)

    mappings = load_mappings()
    obs_df_mage = standardize_metadata(obs_df_mage, dataset_name_for_log, mappings)
    if obs_df_mage is None: # Check if standardize_metadata failed
        logger.error("MAGE: standardize_metadata returned None. Aborting MAGE processing.")
        return None


    adata = create_standard_anndata(unified_expr_df.T, obs_df_mage, var_df_mage, mage_specific_cfg, dataset_name_for_log)

    if mage_specific_cfg:
        adata = apply_dataset_specific_metadata(adata, mage_specific_cfg)

    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed MAGE data: {adata.n_obs} samples and {adata.n_vars} genes")
    else:
        logger.error(f"Failed to save AnnData for MAGE.")
        return None
    return adata


# ====================================================================================
# ADNI-specific Processing
# ====================================================================================

def preprocess_adni_file(file_path):
    """Preprocess ADNI file to fix escaped tab characters if present."""
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        if '\\t' in content:
            logger.info(f"Fixing escaped tabs '\\t' in {file_path}")
            fixed_content = content.replace('\\t', '\t')

            fixed_file_path = file_path + '.fixed_tabs'
            with open(fixed_file_path, 'w', encoding='utf-8') as f_out:
                f_out.write(fixed_content)
            return fixed_file_path
    except Exception as e:
        logger.warning(f"Error during ADNI file preprocessing (tab fixing) for {file_path}: {e}")

    return file_path

def process_adni_data(input_dir, output_file, adni_demographics_file=None, adni_diagnosis_file=None, dict_file=None, metadata_dir=None, gencode_map=None): # Added gencode_map
    logger.info(f"Processing ADNI data from {input_dir} to capture longitudinal info.")
    dataset_name_for_log = "ADNI"
    adni_specific_cfg = load_dataset_specific_metadata(metadata_dir, "adni") or {}

    sample_dfs = {}
    sample_metadata_list = []
    all_gene_ids = set()
    processed_expression_files_count = 0

    demographics_df = None
    if adni_demographics_file and os.path.exists(adni_demographics_file):
        logger.info(f"ADNI: Loading demographics metadata from {adni_demographics_file}")
        try:
            demographics_df = pd.read_csv(adni_demographics_file, low_memory=False)
            if 'RID' in demographics_df.columns:
                demographics_df['RID'] = pd.to_numeric(demographics_df['RID'], errors='coerce').astype('Int64')
            if 'PTID' in demographics_df.columns:
                demographics_df['PTID'] = demographics_df['PTID'].astype(str)
            demographics_df = demographics_df.drop_duplicates(subset=['RID' if 'RID' in demographics_df.columns else 'PTID'], keep='first')
        except Exception as e:
            logger.error(f"ADNI: Error loading demographics file {adni_demographics_file}: {e}", exc_info=True)
            demographics_df = None

    diagnosis_long_df = None
    all_subject_diagnoses_for_uns = {}

    if adni_diagnosis_file and os.path.exists(adni_diagnosis_file):
        logger.info(f"ADNI: Loading longitudinal diagnosis data from {adni_diagnosis_file}")
        try:
            dx_df_raw = pd.read_csv(adni_diagnosis_file, low_memory=False)
            if 'RID' not in dx_df_raw.columns:
                logger.error("ADNI Diagnosis file must contain 'RID' column.")
                all_subject_diagnoses_for_uns = {}
                diagnosis_long_df = None
            else:
                dx_df_raw['EXAMDATE'] = pd.to_datetime(dx_df_raw['EXAMDATE'], errors='coerce')
                diagnosis_long_df = dx_df_raw.copy()
                diagnosis_long_df['RID'] = pd.to_numeric(diagnosis_long_df['RID'], errors='coerce').astype('Int64')
                diagnosis_long_df['DIAGNOSIS'] = pd.to_numeric(diagnosis_long_df['DIAGNOSIS'], errors='coerce')
                if 'VISCODE2' in diagnosis_long_df.columns:
                    diagnosis_long_df['VISCODE2'] = diagnosis_long_df['VISCODE2'].astype(str).str.lower().str.strip()
                logger.info(f"ADNI: Loaded {len(diagnosis_long_df)} total diagnosis entries.")
                diag_map_codes = adni_specific_cfg.get("diagnosis_mapping", {
                    1: "Cognitively Normal", 2: "Mild Cognitive Impairment", 3: "Alzheimer's Disease",
                    -4: "Missing/Unknown"
                })
                dx_df_raw_cleaned = dx_df_raw.dropna(subset=['RID', 'EXAMDATE', 'DIAGNOSIS'])
                dx_df_raw_cleaned.loc[:, 'RID'] = pd.to_numeric(dx_df_raw_cleaned['RID'], errors='coerce').astype('Int64')
                for rid_val, group in dx_df_raw_cleaned.groupby('RID'):
                    subject_diagnoses = []
                    for _, row in group.sort_values(by='EXAMDATE').iterrows():
                        diag_code = row['DIAGNOSIS']
                        diag_label = diag_map_codes.get(diag_code, f"Unknown code: {diag_code}")
                        exam_date_val = row['EXAMDATE']
                        subject_diagnoses.append({
                            "visit_code": str(row.get('VISCODE2', 'N/A')),
                            "exam_date": str(exam_date_val.date()) if pd.notna(exam_date_val) and hasattr(exam_date_val, 'date') else 'N/A',
                            "diagnosis_code": int(diag_code) if pd.notna(diag_code) else None,
                            "diagnosis_label": diag_label
                        })
                    if subject_diagnoses:
                        all_subject_diagnoses_for_uns[str(int(rid_val))] = subject_diagnoses
                logger.info(f"ADNI: Processed full diagnostic history for {len(all_subject_diagnoses_for_uns)} subjects for adata.uns.")
        except Exception as e:
            logger.error(f"ADNI: Error loading or processing diagnosis file {adni_diagnosis_file}: {e}", exc_info=True)
            diagnosis_long_df = None
            all_subject_diagnoses_for_uns = {}

    subject_dirs = glob.glob(os.path.join(input_dir, "*_S_*"))
    expression_files = []
    for subject_dir_path in subject_dirs:
        expression_files.extend(glob.glob(os.path.join(subject_dir_path, "*.txt")))
        expression_files.extend(glob.glob(os.path.join(subject_dir_path, "*.csv")))
        expression_files.extend(glob.glob(os.path.join(subject_dir_path, "**", "*.txt"), recursive=True))
        expression_files.extend(glob.glob(os.path.join(subject_dir_path, "**", "*.csv"), recursive=True))

    expression_files.extend(glob.glob(os.path.join(input_dir, "*.txt")))
    expression_files.extend(glob.glob(os.path.join(input_dir, "*.csv")))
    expression_files = sorted(list(set(expression_files)))

    if not expression_files:
        logger.error(f"No expression files (.txt or .csv) found for ADNI in {input_dir} or its subdirectories.")
        return None

    for file_path_str in expression_files:
        if "ADNI_Gene_Expression_Profile.csv" in os.path.basename(file_path_str) or \
           "ADNI_Gene_Expression_Profile_DICT.csv" in os.path.basename(file_path_str):
            logger.info(f"ADNI: Skipping summary/dictionary file {file_path_str} in per-sample processing loop.")
            continue

        subject_id_from_path = "unknown_subject"
        path_obj = Path(file_path_str)
        parent_dir_name = path_obj.parent.name
        if re.match(r"^\d{3}_S_\d{4,}$", parent_dir_name):
            subject_id_from_path = parent_dir_name
        else:
            match_fn = re.search(r"(\d{3}_S_\d{4,})", path_obj.name)
            if match_fn:
                subject_id_from_path = match_fn.group(1)
            else:
                logger.warning(f"Could not determine ADNI subject ID for file: {file_path_str}. Using 'unknown_subject'.")

        file_name = os.path.basename(file_path_str)
        sample_id_suffix = Path(file_name).stem
        sample_id = f"{subject_id_from_path}_{sample_id_suffix}"

        processed_file_path = preprocess_adni_file(file_path_str)

        try:
            df = None
            file_basename_lower = os.path.basename(processed_file_path).lower()
            if "_gencode_v24_pruned.csv" in file_basename_lower:
                logger.debug(f"ADNI: Attempting tab delimiter for pruned CSV: {processed_file_path}")
                try: df = pd.read_csv(processed_file_path, sep='\t', comment='#', low_memory=False)
                except Exception:
                    logger.warning(f"Tab delimiter failed for pruned CSV {processed_file_path}. Trying comma.")
                    try: df = pd.read_csv(processed_file_path, sep=',', comment='#', low_memory=False)
                    except Exception as e_comma:
                        logger.error(f"Comma delimiter also failed for pruned CSV {processed_file_path}: {e_comma}. Skipping.")
                        if processed_file_path != file_path_str: os.remove(processed_file_path)
                        continue
            elif file_basename_lower.endswith('.csv'):
                try: df = pd.read_csv(processed_file_path, sep=',', comment='#', low_memory=False)
                except Exception:
                    logger.warning(f"Comma delimiter failed for CSV {processed_file_path}. Trying tab.")
                    try: df = pd.read_csv(processed_file_path, sep='\t', comment='#', low_memory=False)
                    except Exception as e_tab:
                        logger.error(f"Tab delimiter also failed for CSV {processed_file_path}: {e_tab}. Skipping.")
                        if processed_file_path != file_path_str: os.remove(processed_file_path)
                        continue
            else:
                df = pd.read_csv(processed_file_path, sep='\t', comment='#', low_memory=False)

            if df is None:
                if processed_file_path != file_path_str: os.remove(processed_file_path)
                continue
            if processed_file_path != file_path_str and os.path.exists(processed_file_path):
                os.remove(processed_file_path)

            if df.empty:
                logger.warning(f"Empty data after reading ADNI file: {file_path_str}")
                continue

            gene_id_col, expr_val_col = identify_columns(df)
            if gene_id_col is None or expr_val_col is None:
                logger.warning(f"Could not identify required columns in ADNI file {file_path_str}. Gene: {gene_id_col}, Expr: {expr_val_col}. Columns: {df.columns.tolist()}")
                continue

            df[gene_id_col] = df[gene_id_col].astype(str)
            df[expr_val_col] = pd.to_numeric(df[expr_val_col], errors='coerce')
            simple_df = df.set_index(gene_id_col)[[expr_val_col]]
            simple_df = simple_df.rename(columns={expr_val_col: sample_id})
            sample_dfs[sample_id] = simple_df
            all_gene_ids.update(simple_df.index)

            current_metadata = {
                "sample_id_from_file": sample_id, "subject_id": subject_id_from_path,
                "dataset": adni_specific_cfg.get("label", "ADNI"),
                "data_type": adni_specific_cfg.get("data_type", "Microarray"),
                "expression_unit": adni_specific_cfg.get("obs_columns", {}).get("expression_unit", "Normalized intensity"),
                "tissue": adni_specific_cfg.get("obs_columns", {}).get("tissue", "blood"),
                "platform": adni_specific_cfg.get("platform", "Affymetrix Human Genome U219 Array"),
                "expression_sample_visit_code": "unknown", "age_at_expression_sampling": "",
                "diagnosis_at_expression_sampling_code": pd.NA,
                "diagnosis_at_expression_sampling_label": "Unknown",
                "date_of_expression_sampling": ""
            }
            parsed_visit_code = "unknown"
            filename_lower = sample_id_suffix.lower()
            visit_pattern = r"_(bl|sc|scrn|screening|m\d{1,3}|v\d{1,3}|aim\d{1,3})[._]"
            visit_match = re.search(visit_pattern, filename_lower)
            if visit_match:
                parsed_visit_code = visit_match.group(1)
                if parsed_visit_code.startswith("aim"): parsed_visit_code = "m" + parsed_visit_code[3:].zfill(2)
                elif parsed_visit_code in ["scrn", "screening"]: parsed_visit_code = "sc"
            current_metadata["expression_sample_visit_code"] = parsed_visit_code
            subject_rid = None
            try:
                if "_S_" in subject_id_from_path: subject_rid = int(subject_id_from_path.split('_S_')[-1])
            except ValueError: logger.warning(f"Could not parse RID for {subject_id_from_path}")
            if demographics_df is not None and subject_rid is not None:
                demog_subject_row = demographics_df[demographics_df['RID'] == subject_rid]
                if not demog_subject_row.empty:
                    demog_data = demog_subject_row.iloc[0]
                    sex_col_demog = next((c for c in ['PTGENDER', 'SEX'] if c in demog_data.index), None)
                    if sex_col_demog and pd.notna(demog_data[sex_col_demog]):
                        sex_val = demog_data[sex_col_demog]
                        sex_map = {1:'male', 2:'female', '1':'male', '2':'female', 'Male':'male', 'Female':'female'}
                        current_metadata['sex'] = sex_map.get(str(sex_val).lower(), sex_map.get(int(sex_val) if isinstance(sex_val, (int,float)) or str(sex_val).isdigit() else None, 'unknown'))
                    current_metadata['PTDOBYY_temp'] = pd.to_numeric(demog_data.get('PTDOBYY'), errors='coerce')
                    if adni_specific_cfg.get('race_source_column') in demog_data:
                        race_val = demog_data[adni_specific_cfg['race_source_column']]
                        race_map = adni_specific_cfg.get('source_race_to_standard_label_map',{})
                        current_metadata['race'] = race_map.get(str(race_val), race_map.get(int(race_val) if pd.notna(race_val) and str(race_val).isdigit() else race_val, 'unknown or not reported'))
                    else: current_metadata['race'] = 'unknown or not reported'
                    if adni_specific_cfg.get('ethnicity_source_column') in demog_data:
                        eth_val = demog_data[adni_specific_cfg['ethnicity_source_column']]
                        eth_map = adni_specific_cfg.get('source_ethnicity_to_standard_label_map',{})
                        current_metadata['is_hispanic_or_latino'] = eth_map.get(str(eth_val), eth_map.get(int(eth_val) if pd.notna(eth_val) and str(eth_val).isdigit() else eth_val, 'unknown or not reported'))
                    else: current_metadata['is_hispanic_or_latino'] = 'unknown or not reported'
                else: logger.warning(f"ADNI: No demographics for RID {subject_rid}")
            exam_date_for_this_visit = None
            if diagnosis_long_df is not None and subject_rid is not None and current_metadata["expression_sample_visit_code"] != "unknown":
                visit_dx_entry = diagnosis_long_df[(diagnosis_long_df['RID'] == subject_rid) & (diagnosis_long_df['VISCODE2'] == current_metadata["expression_sample_visit_code"])]
                if not visit_dx_entry.empty:
                    diag_data = visit_dx_entry.iloc[0]
                    diag_code = diag_data.get('DIAGNOSIS'); exam_date_for_this_visit = diag_data.get('EXAMDATE')
                    diag_map = adni_specific_cfg.get("diagnosis_mapping", {1:"CN", 2:"MCI", 3:"AD", -4:"Unknown"})
                    current_metadata['diagnosis_at_expression_sampling_label'] = diag_map.get(diag_code, "Unknown")
                    current_metadata['diagnosis_at_expression_sampling_code'] = diag_code if pd.notna(diag_code) else pd.NA
                    current_metadata['date_of_expression_sampling'] = str(exam_date_for_this_visit.date()) if pd.notna(exam_date_for_this_visit) and hasattr(exam_date_for_this_visit, 'date') else ""
                else: logger.warning(f"ADNI: No DXSUM for RID {subject_rid}, VISCODE2 {current_metadata['expression_sample_visit_code']}")
            elif diagnosis_long_df is None: logger.debug("ADNI: DXSUM not loaded/RID missing.")
            if pd.notna(current_metadata.get('PTDOBYY_temp')) and pd.notna(exam_date_for_this_visit) and hasattr(exam_date_for_this_visit, 'year'):
                current_metadata['age_at_expression_sampling'] = str(exam_date_for_this_visit.year - int(current_metadata['PTDOBYY_temp']))
            elif demographics_df is not None and subject_rid is not None:
                 demog_row = demographics_df[demographics_df['RID'] == subject_rid]
                 if not demog_row.empty and 'AGE' in demog_row.columns and pd.notna(demog_row['AGE'].iloc[0]):
                    age_baseline = pd.to_numeric(demog_row['AGE'].iloc[0], errors='coerce')
                    if pd.notna(age_baseline): current_metadata['age_at_expression_sampling'] = str(int(age_baseline))
            current_metadata.pop('PTDOBYY_temp', None)
            sample_metadata_list.append(current_metadata)
            processed_expression_files_count += 1
        except Exception as e:
            logger.error(f"Error processing ADNI file {file_path_str} for sample {sample_id}: {e}", exc_info=True)
            if processed_file_path != file_path_str and os.path.exists(processed_file_path):
                os.remove(processed_file_path)

    if processed_expression_files_count == 0:
        logger.error(f"No valid ADNI expression data processed from {input_dir}")
        return None
    logger.info(f"Processed {processed_expression_files_count} ADNI expression files.")

    obs_df_adni = pd.DataFrame(sample_metadata_list)
    if 'sample_id_from_file' not in obs_df_adni.columns or obs_df_adni.empty:
        logger.error(f"{dataset_name_for_log}: 'sample_id_from_file' key missing or obs_df is empty.")
        return None
    obs_df_adni = obs_df_adni.set_index('sample_id_from_file')
    obs_df_adni.index.name = "sample_id"
    final_sample_ids_adni = obs_df_adni.index.tolist()
    if obs_df_adni.index.has_duplicates:
        logger.warning(f"{dataset_name_for_log}: Correcting duplicate sample IDs in obs_df index.")
        obs_df_adni = obs_df_adni[~obs_df_adni.index.duplicated(keep='first')]
        final_sample_ids_adni = obs_df_adni.index.tolist()
    sample_dfs = {sid: df_val for sid, df_val in sample_dfs.items() if sid in final_sample_ids_adni}

    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs for {dataset_name_for_log}")
    gene_id_mapping = {gid: standardize_ensembl_id(gid) for gid in all_gene_ids}
    valid_mapped_values_adni = [str(gid) for gid in gene_id_mapping.values() if gid is not None]
    if not valid_mapped_values_adni:
        logger.error(f"{dataset_name_for_log}: No gene IDs could be meaningfully standardized. Aborting.")
        return None
    unique_std_ids = sorted(list(set(valid_mapped_values_adni)))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs for {dataset_name_for_log}.")

    # Optimized Matrix Creation for ADNI
    logger.info(f"Creating unified expression matrix for {dataset_name_for_log} ({len(unique_std_ids)} genes x {len(final_sample_ids_adni)} samples) using optimized method.")
    all_series = [s_df[sample_id].rename(sample_id) for sample_id, s_df in sample_dfs.items()]
    if not all_series:
        logger.error(f"No data series collected for {dataset_name_for_log}. Aborting matrix creation.")
        unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_adni).astype(np.float32) # Create empty if no series
    else:
        try:
            concat_df = pd.concat(all_series, axis=1, join='outer').astype(np.float32)
            logger.info(f"Concatenated DataFrame shape: {concat_df.shape}")

            std_to_originals_map_adni = {}
            for original_gid, mapped_std_id_val in gene_id_mapping.items():
                if mapped_std_id_val is not None:
                    std_to_originals_map_adni.setdefault(str(mapped_std_id_val), []).append(str(original_gid))

            original_to_std_series_map = {orig_id: std_id_val for std_id_val, orig_ids_list in std_to_originals_map_adni.items() for orig_id in orig_ids_list}
            mappable_original_ids = [idx for idx in concat_df.index if idx in original_to_std_series_map]

            if not mappable_original_ids:
                logger.warning(f"No mappable original gene IDs in concat_df for {dataset_name_for_log}. Resulting matrix might be all zeros or have unexpected shape.")
                unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_adni).astype(np.float32)
            else:
                filtered_concat_df = concat_df.loc[mappable_original_ids]
                std_ids_for_aggregation = filtered_concat_df.index.map(original_to_std_series_map)
                if std_ids_for_aggregation.isna().all():
                    logger.error(f"All standardized IDs for aggregation are NaN for {dataset_name_for_log}. Check mapping.")
                    unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_adni).astype(np.float32)
                else:
                    filtered_concat_df.index = std_ids_for_aggregation
                    unified_expr_df_grouped = filtered_concat_df.groupby(level=0).max()
                    unified_expr_df = unified_expr_df_grouped.reindex(index=unique_std_ids, columns=final_sample_ids_adni, fill_value=0.0).astype(np.float32)
        except Exception as e_concat:
            logger.error(f"Error during optimized matrix creation for {dataset_name_for_log}: {e_concat}", exc_info=True)
            logger.warning("Falling back to previous (slower) matrix creation method for ADNI due to error.")
            # Fallback to original loop method
            expr_matrix = np.full((len(unique_std_ids), len(final_sample_ids_adni)), np.nan, dtype=np.float32)
            sample_id_to_col_idx = {sid: i for i, sid in enumerate(final_sample_ids_adni)}
            std_id_to_idx = {gid: i for i, gid in enumerate(unique_std_ids)} # Re-define if not in scope
            std_to_originals_map_adni = {} # Re-create as it might not be in scope if optimized failed early
            for original_gid, mapped_std_id_val in gene_id_mapping.items():
                if mapped_std_id_val is not None: std_to_originals_map_adni.setdefault(str(mapped_std_id_val), []).append(str(original_gid))

            for std_id, original_ids_for_this_std in std_to_originals_map_adni.items():
                if std_id not in std_id_to_idx: continue
                gene_idx = std_id_to_idx[std_id]
                for sample_id, col_idx in sample_id_to_col_idx.items():
                    sample_df = sample_dfs.get(sample_id)
                    if sample_df is None: continue
                    if sample_id not in sample_df.columns: continue
                    relevant_orig_ids = [orig_id for orig_id in original_ids_for_this_std if orig_id in sample_df.index]
                    if relevant_orig_ids:
                        values = pd.to_numeric(sample_df.loc[relevant_orig_ids, sample_id], errors='coerce')
                        if not values.empty and not values.isna().all(): expr_matrix[gene_idx, col_idx] = values.max()
            expr_matrix = np.nan_to_num(expr_matrix, nan=0.0)
            unified_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=final_sample_ids_adni)

    logger.info(f"Optimized ADNI matrix created. Shape: {unified_expr_df.shape}")
    # End of Optimized Matrix Creation

    var_df_adni = pd.DataFrame(index=unique_std_ids)
    var_df_adni["gene_id"] = var_df_adni.index.astype(str)
    var_df_adni["original_ids"] = var_df_adni.index.map(lambda x: ";".join(std_to_originals_map_adni.get(x, [])))

    # gencode_map passed as parameter
    var_df_adni = add_gencode_annotations(var_df_adni, gencode_map)

    mappings = load_mappings()
    # Add worst diagnosis over time before metadata standardization
    if all_subject_diagnoses_for_uns:
        logger.info("ADNI: Computing worst diagnosis over time for each subject...")
        # Define diagnosis severity order based on ADNI progression model:
        # 1 = CN (Cognitively Normal) - least severe
        # 2 = MCI (Mild Cognitive Impairment) - moderate severity  
        # 3 = AD (Alzheimer's Disease) - most severe
        # -4 = Unknown/Missing - lowest priority (severity 0)
        diagnosis_severity_order = {1: 1, 2: 2, 3: 3, -4: 0}
        
        worst_diagnosis_map = {}
        for subject_rid_str, diagnosis_history in all_subject_diagnoses_for_uns.items():
            subject_rid = int(subject_rid_str)
            valid_diagnoses = [entry for entry in diagnosis_history if entry.get('diagnosis_code') is not None]
            
            if valid_diagnoses:
                # Find diagnosis with highest severity
                worst_entry = max(valid_diagnoses, 
                                key=lambda x: diagnosis_severity_order.get(x['diagnosis_code'], 0))
                worst_diagnosis_map[subject_rid] = {
                    'worst_diagnosis_code': worst_entry['diagnosis_code'],
                    'worst_diagnosis_label': worst_entry['diagnosis_label'],
                    'worst_diagnosis_visit': worst_entry['visit_code'],
                    'worst_diagnosis_date': worst_entry['exam_date']
                }
            else:
                worst_diagnosis_map[subject_rid] = {
                    'worst_diagnosis_code': -4,
                    'worst_diagnosis_label': 'Unknown',
                    'worst_diagnosis_visit': 'N/A',
                    'worst_diagnosis_date': 'N/A'
                }
        
        # Add worst diagnosis fields to obs_df_adni
        for col in ['worst_diagnosis_code', 'worst_diagnosis_label', 'worst_diagnosis_visit', 'worst_diagnosis_date']:
            obs_df_adni[col] = pd.NA
        
        for idx, row in obs_df_adni.iterrows():
            subject_id = row['subject_id']
            if subject_id and "_S_" in subject_id:
                try:
                    subject_rid = int(subject_id.split('_S_')[-1])
                    if subject_rid in worst_diagnosis_map:
                        worst_data = worst_diagnosis_map[subject_rid]
                        obs_df_adni.loc[idx, 'worst_diagnosis_code'] = worst_data['worst_diagnosis_code']
                        obs_df_adni.loc[idx, 'worst_diagnosis_label'] = worst_data['worst_diagnosis_label']
                        obs_df_adni.loc[idx, 'worst_diagnosis_visit'] = worst_data['worst_diagnosis_visit']
                        obs_df_adni.loc[idx, 'worst_diagnosis_date'] = worst_data['worst_diagnosis_date']
                except ValueError:
                    continue
        
        logger.info(f"ADNI: Added worst diagnosis over time for {len(worst_diagnosis_map)} subjects")
    
    obs_df_adni = standardize_metadata(obs_df_adni, dataset_name_for_log, mappings)
    if obs_df_adni is None:
        logger.error("ADNI: standardize_metadata returned None. Aborting ADNI processing.")
        return None

    obs_df_adni.rename(columns={
        'age_at_expression_sampling': 'age',
        'diagnosis_at_expression_sampling_code': 'diagnosis_code',
        'diagnosis_at_expression_sampling_label': 'diagnosis',
        'date_of_expression_sampling': 'visit_date'
    }, inplace=True)

    for col, dtype, default in [
        ('age', str, ''), ('sex', 'category', 'unknown'),
        ('race', 'category', 'unknown or not reported'),
        ('is_hispanic_or_latino', 'category', 'unknown or not reported'),
        ('self_reported_ethnicity', 'category', 'unknown or not reported'),
        ('self_reported_ethnicity_ontology_term_id', 'category', 'unknown'),
        ('diagnosis', 'category', 'Unknown'),
        ('diagnosis_code', 'Int64', pd.NA),
        ('worst_diagnosis_code', 'Int64', pd.NA),
        ('worst_diagnosis_label', 'category', 'Unknown'),
        ('worst_diagnosis_visit', str, 'N/A'),
        ('worst_diagnosis_date', str, 'N/A'),
        ('visit_date', str, ''),
        ('expression_sample_visit_code', str, 'unknown')
    ]:
        if col not in obs_df_adni.columns:
            obs_df_adni[col] = default
        if dtype == 'category':
            current_categories = None
            if col == 'race': current_categories = NIH_RACE_CATEGORIES_LOWER
            elif col == 'is_hispanic_or_latino': current_categories = NIH_ETHNICITY_CATEGORIES_LOWER
            elif col == 'self_reported_ethnicity': current_categories = SRE_BASE_CATEGORIES_LOWER
            elif col == 'self_reported_ethnicity_ontology_term_id': pass # Handled by standardize_metadata

            if current_categories:
                 obs_df_adni[col] = pd.Categorical(
                     obs_df_adni[col].astype(str).fillna(str(default) if pd.notna(default) else ''),
                     categories=current_categories, ordered=False)
            elif col != 'self_reported_ethnicity_ontology_term_id':
                 obs_df_adni[col] = obs_df_adni[col].astype(str).fillna(str(default) if pd.notna(default) else '').astype('category')
        elif dtype == 'Int64':
            obs_df_adni[col] = pd.to_numeric(obs_df_adni[col], errors='coerce').astype('Int64')
        else:
            obs_df_adni[col] = obs_df_adni[col].astype(str).fillna(str(default) if pd.notna(default) else '')

    # ---- Moved Deduplication Logic Here (ADNI) ----
    if obs_df_adni.columns.has_duplicates:
        duplicated_cols_obs = obs_df_adni.columns[obs_df_adni.columns.duplicated(keep=False)].unique().tolist()
        logger.warning(f"ADNI: obs_df_adni has duplicate columns before creating AnnData. Duplicates: {duplicated_cols_obs}")
        obs_df_adni = obs_df_adni.loc[:, ~obs_df_adni.columns.duplicated(keep='first')]
        logger.info(f"ADNI: Deduplicated obs_df_adni columns. New shape: {obs_df_adni.shape}, Columns: {obs_df_adni.columns.tolist()}")

    if var_df_adni.columns.has_duplicates:
        duplicated_cols_var = var_df_adni.columns[var_df_adni.columns.duplicated(keep=False)].unique().tolist()
        logger.warning(f"ADNI: var_df_adni has duplicate columns before creating AnnData. Duplicates: {duplicated_cols_var}")
        var_df_adni = var_df_adni.loc[:, ~var_df_adni.columns.duplicated(keep='first')]
        logger.info(f"ADNI: Deduplicated var_df_adni columns. New shape: {var_df_adni.shape}, Columns: {var_df_adni.columns.tolist()}")
    # ---- End Moved Deduplication Logic ----

    logger.info(f"ADNI: Columns in obs_df_adni before AnnData creation: {obs_df_adni.columns.tolist()}")
    logger.info(f"ADNI: dtypes in obs_df_adni before AnnData creation:\n{obs_df_adni.dtypes}")
    for col_debug in obs_df_adni.columns:
        if isinstance(obs_df_adni[col_debug], pd.DataFrame):
            logger.error(f"ADNI ERROR (FINAL CHECK): obs_df_adni['{col_debug}'] is still a DataFrame!")
    for col_debug_var in var_df_adni.columns:
        if isinstance(var_df_adni[col_debug_var], pd.DataFrame):
            logger.error(f"ADNI ERROR (FINAL CHECK): var_df_adni['{col_debug_var}'] is still a DataFrame!")

    adata = create_standard_anndata(unified_expr_df.T, obs_df_adni, var_df_adni, adni_specific_cfg, dataset_name_for_log)

    if all_subject_diagnoses_for_uns:
        adata.uns['longitudinal_diagnoses_adni'] = all_subject_diagnoses_for_uns
        logger.info("ADNI: Added full longitudinal diagnosis history to adata.uns['longitudinal_diagnoses_adni']")

    if adni_specific_cfg:
        adata = apply_dataset_specific_metadata(adata, adni_specific_cfg)

    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed ADNI data: {adata.n_obs} samples and {adata.n_vars} genes")
    else:
        logger.error(f"Failed to save AnnData for ADNI.")
        return None
    return adata


# ====================================================================================
# Main Processing Pipeline
# ====================================================================================
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Multi-Dataset Standardization Pipeline (Stage 1)")
    parser.add_argument("--encode-dir", help="Directory containing ENCODE cell line TPM files")
    parser.add_argument("--encode-entex-dir", help="Directory containing ENCODE ENTEx TPM files (can be same as encode-dir if co-located)")
    parser.add_argument("--gtex-file", help="Path to GTEx expression file (GCT format)")
    parser.add_argument("--mage-dir", help="Directory containing MAGE expression files")
    parser.add_argument(
        "--adni-dir", help="Directory containing ADNI sample directories with expression files"
    )
    parser.add_argument(
        "--metadata-dir",
        default=str(DEFAULT_METADATA_DIR),
        help=f"Directory containing metadata JSON files (default: {DEFAULT_METADATA_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"Output directory for standardized files (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument("--entex-metadata-file", help="Path to ENTEx metadata JSON file (can override default path if entex-dir is used)")
    parser.add_argument("--adni-demographics-file", help="Path to ADNI patient demographics CSV (e.g., PTDEMOG_*.csv)")
    parser.add_argument("--adni-diagnosis-file", help="Path to ADNI diagnosis summary CSV (e.g., DXSUM_*.csv)")
    return parser.parse_args()

def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    start_time = time.time()

    # Load GENCODE mapping once
    logger.info("Loading GENCODE mapping for the pipeline run...")
    gencode_map = load_gencode_mapping()
    if not gencode_map:
        logger.error("Failed to load GENCODE mapping. This is critical for gene annotations. Exiting.")
        sys.exit(1)
    logger.info("GENCODE mapping loaded successfully.")


    metadata_dir = Path(args.metadata_dir if args.metadata_dir else DEFAULT_METADATA_DIR)
    processed_datasets_summary = []

    # Process ENCODE (which can include ENTEx based on its args)
    if args.encode_dir or args.encode_entex_dir or args.entex_metadata_file:
        encode_output_file = Path(args.output_dir) / "encode_standardized_v1.h5ad"
        logger.info("Starting ENCODE/ENTEx processing...")
        encode_data = process_encode_data(
            args.encode_dir,
            args.encode_entex_dir,
            str(encode_output_file),
            args.entex_metadata_file,
            str(metadata_dir),
            gencode_map=gencode_map # Pass loaded map
        )
        if encode_data:
            processed_label = "ENCODE"
            if args.encode_entex_dir or args.entex_metadata_file:
                processed_label = "ENCODE_ENTEx_Combined" if args.encode_dir else "ENTEx"

            processed_datasets_summary.append((processed_label, encode_data.n_obs, encode_data.n_vars))
            logger.info(
                f"Successfully processed {processed_label} data: {encode_data.n_obs} samples and {encode_data.n_vars} genes"
            )
            if hasattr(encode_data, 'var') and 'mapping_source' in encode_data.var.columns:
                logger.info(f"{processed_label} gene statistics (mapping source):")
                gencode_stats = encode_data.var["mapping_source"].value_counts()
                for source, count in gencode_stats.items():
                    logger.info(f"  {source}: {count} genes ({count/encode_data.n_vars:.1%})")


    # Process GTEx
    if args.gtex_file:
        gtex_output_file = Path(args.output_dir) / "gtex_standardized_v1.h5ad"
        logger.info(f"Processing GTEx data from {args.gtex_file}")
        gtex_data = process_gtex_data(args.gtex_file, str(gtex_output_file), str(metadata_dir), gencode_map=gencode_map) # Pass loaded map
        if gtex_data:
            processed_datasets_summary.append(("GTEx", gtex_data.n_obs, gtex_data.n_vars))
            logger.info(
                f"Successfully processed GTEx data: {gtex_data.n_obs} samples, {gtex_data.n_vars} genes"
            )
            if hasattr(gtex_data, 'var') and 'mapping_source' in gtex_data.var.columns:
                logger.info("GTEx gene statistics (mapping source):")
                gencode_stats = gtex_data.var["mapping_source"].value_counts()
                for source, count in gencode_stats.items():
                    logger.info(f"  {source}: {count} genes ({count/gtex_data.n_vars:.1%})")

    # Process MAGE
    if args.mage_dir:
        mage_output_file = Path(args.output_dir) / "mage_standardized_v1.h5ad"
        logger.info(f"Processing MAGE data from {args.mage_dir}")
        mage_data = process_mage_data(args.mage_dir, str(mage_output_file), str(metadata_dir), gencode_map=gencode_map) # Pass loaded map
        if mage_data:
            processed_datasets_summary.append(("MAGE", mage_data.n_obs, mage_data.n_vars))
            logger.info(
                f"Successfully processed MAGE data: {mage_data.n_obs} samples, {mage_data.n_vars} genes"
            )

    # Process ADNI
    if args.adni_dir:
        adni_output_file = Path(args.output_dir) / "adni_standardized_v1.h5ad"
        logger.info(f"Processing ADNI data from directory: {args.adni_dir}")
        adni_data = process_adni_data(
            args.adni_dir,
            str(adni_output_file),
            adni_demographics_file=args.adni_demographics_file,
            adni_diagnosis_file=args.adni_diagnosis_file,
            dict_file=None,
            metadata_dir=str(metadata_dir),
            gencode_map=gencode_map # Pass loaded map
        )
        if adni_data:
            processed_datasets_summary.append(("ADNI", adni_data.n_obs, adni_data.n_vars))
            logger.info(
                f"Successfully processed ADNI data: {adni_data.n_obs} samples, {adni_data.n_vars} genes"
            )

    logger.info("=== Stage 1 Processing Summary ===")
    if processed_datasets_summary:
        for name, obs, var in processed_datasets_summary:
            logger.info(f"  {name}: {obs} samples, {var} genes.")
    else:
        logger.info("  No datasets were processed in this run based on provided arguments.")

    total_time = time.time() - start_time
    logger.info(f"Total Stage 1 processing time: {total_time:.2f} seconds")


if __name__ == "__main__":
    main()