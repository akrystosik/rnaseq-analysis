#!/usr/bin/env python3
"""
RNA-seq Metadata Standardization Script (V2 - JSON Prioritized)

This script standardizes metadata for RNA-seq datasets, adding consistent
annotations, ontology mappings, and applying dataset-specific metadata
from JSON configuration files. Inference is used primarily as a fallback.
"""

import os
import re
import logging
import argparse
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import numpy as np
import json # <-- Add json import

# Import shared utilities
from rnaseq_utils import (
    # We might not need standardize_metadata here anymore if Stage 1 is sufficient
    load_dataset_specific_metadata,
    apply_dataset_specific_metadata,
    load_mappings # Keep for ontology lookups if needed
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('metadata_standardization_v2')

# Define paths and constants
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
DEFAULT_METADATA_DIR = BASE_DIR / "metadata/json"
DEFAULT_OUTPUT_DIR = BASE_DIR / "standardized_data"

# --- TissueOntologyMapper Class (Keep as is from previous version) ---
class TissueOntologyMapper:
    """
    A class to handle mapping of tissue terms to standardized ontology terms.
    Loads mappings from JSON files (tissue_to_uberon.json by default).
    """
    def _load_json_mapping(self, mapping_file):
        """Load mappings from a JSON file"""
        try:
            import json
            with open(mapping_file, 'r') as f:
                mapping_data = json.load(f)
                # Handle different JSON structures (simple key:value or key:dict)
                for tissue, ontology_info in mapping_data.items():
                    if isinstance(ontology_info, dict) and 'id' in ontology_info:
                        self.tissue_mappings[tissue.lower()] = ontology_info['id'] # Store keys lowercase
                        self.mapping_confidence[tissue.lower()] = ontology_info.get('confidence', 'medium')
                    elif isinstance(ontology_info, str) and (ontology_info.startswith("UBERON:") or ontology_info.startswith("CL:")):
                        self.tissue_mappings[tissue.lower()] = ontology_info # Store keys lowercase
                        self.mapping_confidence[tissue.lower()] = 'medium' # Default confidence
                    else:
                        logger.debug(f"Skipping invalid mapping entry: {tissue}: {ontology_info}")

            logger.info(f"Loaded {len(self.tissue_mappings)} tissue mappings from {mapping_file}")
        except FileNotFoundError:
             logger.warning(f"Tissue mapping file not found: {mapping_file}. Skipping.")
        except Exception as e:
            logger.error(f"Error loading tissue mappings from JSON {mapping_file}: {e}")

    def __init__(self, default_mapping_dir=DEFAULT_METADATA_DIR):
        """
        Initialize the tissue mapper.
        Args:
            default_mapping_dir: Directory containing the default tissue_to_uberon.json
        """
        self.tissue_mappings = {}
        self.mapping_confidence = {}
        self.unmapped_tissues = set() # Track tissues we couldn't map

        # Load the default JSON mapping file
        default_mapping_file = Path(default_mapping_dir) / "tissue_to_uberon.json"
        if default_mapping_file.exists():
            self._load_json_mapping(default_mapping_file)
        else:
            logger.warning(f"Default tissue mapping file not found at {default_mapping_file}")

    def map_tissue(self, tissue_name):
        """
        Map a tissue name to an ontology ID.
        Returns (ontology_id, confidence) tuple.
        """
        # Handle edge cases
        if pd.isna(tissue_name) or not tissue_name or str(tissue_name).strip() == '':
            return '', 'none'

        # Normalize tissue name (lowercase, strip)
        tissue_lower = str(tissue_name).strip().lower()

        # Check exact lowercase match
        if tissue_lower in self.tissue_mappings:
            return (self.tissue_mappings[tissue_lower],
                    self.mapping_confidence.get(tissue_lower, 'medium'))

        # If not found directly, mark as unmapped for now
        self.unmapped_tissues.add(tissue_name) # Store original case
        return '', 'none'

    def map_tissues_in_adata(self, adata, tissue_field='tissue'):
        """Add ontology mappings to an AnnData object based on loaded mappings."""
        if tissue_field not in adata.obs.columns:
            logger.warning(f"Tissue field '{tissue_field}' not found in adata.obs. Skipping ontology mapping.")
            # Ensure columns exist even if mapping is skipped
            if 'tissue_ontology' not in adata.obs:
                adata.obs['tissue_ontology'] = ''
            if 'tissue_ontology_confidence' not in adata.obs:
                adata.obs['tissue_ontology_confidence'] = ''
            return adata

        logger.info(f"Applying tissue ontology mappings using field '{tissue_field}'")

        # Initialize columns if they don't exist, ensuring string type
        if 'tissue_ontology' not in adata.obs.columns:
            adata.obs['tissue_ontology'] = pd.Series(index=adata.obs.index, dtype='str')
        else: # Ensure existing column is string type and fill NaNs
             adata.obs['tissue_ontology'] = adata.obs['tissue_ontology'].astype(str).fillna('')

        if 'tissue_ontology_confidence' not in adata.obs.columns:
            adata.obs['tissue_ontology_confidence'] = pd.Series(index=adata.obs.index, dtype='str')
        else: # Ensure existing column is string type and fill NaNs
             adata.obs['tissue_ontology_confidence'] = adata.obs['tissue_ontology_confidence'].astype(str).fillna('none')

        # --- Efficient Mapping using pandas map ---
        # Create mapping series from dictionaries
        tissue_map_series = pd.Series(self.tissue_mappings)
        confidence_map_series = pd.Series(self.mapping_confidence)

        # Normalize the tissue column for mapping (lowercase, strip)
        normalized_tissue_col = adata.obs[tissue_field].astype(str).str.strip().str.lower()

        # Apply mapping
        adata.obs['tissue_ontology'] = normalized_tissue_col.map(tissue_map_series).fillna('')
        adata.obs['tissue_ontology_confidence'] = normalized_tissue_col.map(confidence_map_series).fillna('none')

        # Identify unmapped tissues (where ontology is empty but original tissue was not)
        unmapped_mask = (adata.obs['tissue_ontology'] == '') & (adata.obs[tissue_field].astype(str).str.strip() != '') & (adata.obs[tissue_field].notna())
        self.unmapped_tissues.update(adata.obs.loc[unmapped_mask, tissue_field].unique())

        # Log stats
        total_samples = adata.n_obs
        mapped_samples = (adata.obs['tissue_ontology'] != '').sum()
        mapping_percentage = (mapped_samples / total_samples) * 100 if total_samples > 0 else 0
        confidence_counts = adata.obs['tissue_ontology_confidence'].value_counts().to_dict()

        logger.info(f"Tissue ontology mapping complete: {mapped_samples}/{total_samples} ({mapping_percentage:.1f}%) mapped.")
        logger.info(f"Mapping confidence counts: {confidence_counts}")

        if self.unmapped_tissues:
             # Log only a sample if there are too many
            sample_unmapped = list(self.unmapped_tissues)[:20]
            logger.warning(f"Found {len(self.unmapped_tissues)} unmapped tissues. Examples: {sample_unmapped}")

        # Convert back to categorical if desired (optional)
        # adata.obs['tissue_ontology'] = pd.Categorical(adata.obs['tissue_ontology'])
        # adata.obs['tissue_ontology_confidence'] = pd.Categorical(adata.obs['tissue_ontology_confidence'])

        return adata

# --- CellTypeOntologyMapper Class ---
class CellTypeOntologyMapper:
    """
    A class to handle mapping of cell type/cell line names to Cell Ontology (CL) terms.
    Loads mappings from celltype_to_cl.json.
    """
    def _load_json_mapping(self, mapping_file):
        """Load cell type mappings from a JSON file"""
        try:
            with open(mapping_file, 'r') as f:
                mapping_data = json.load(f)
                mappings = mapping_data.get('mappings', {})
                
                for cell_type, ontology_info in mappings.items():
                    # Skip comment entries
                    if cell_type.startswith('_comment'):
                        continue
                        
                    if isinstance(ontology_info, dict):
                        cl_term_id = ontology_info.get('cl_term_id', '')
                        cl_term_name = ontology_info.get('cl_term_name', '')
                        confidence = ontology_info.get('confidence', 'medium')
                        
                        # Store the mapping (case-insensitive key)
                        self.celltype_mappings[cell_type.lower()] = cl_term_id
                        self.celltype_names[cell_type.lower()] = cl_term_name
                        self.mapping_confidence[cell_type.lower()] = confidence
                    else:
                        logger.debug(f"Skipping invalid cell type mapping entry: {cell_type}: {ontology_info}")

            logger.info(f"Loaded {len(self.celltype_mappings)} cell type mappings from {mapping_file}")
        except FileNotFoundError:
             logger.warning(f"Cell type mapping file not found: {mapping_file}. Skipping.")
        except Exception as e:
            logger.error(f"Error loading cell type mappings from JSON {mapping_file}: {e}")

    def __init__(self, default_mapping_dir=DEFAULT_METADATA_DIR):
        """
        Initialize the cell type mapper.
        Args:
            default_mapping_dir: Directory containing the celltype_to_cl.json
        """
        self.celltype_mappings = {}
        self.celltype_names = {}
        self.mapping_confidence = {}
        self.unmapped_celltypes = set()  # Track cell types we couldn't map

        # Load the default JSON mapping file
        default_mapping_file = Path(default_mapping_dir) / "celltype_to_cl.json"
        if default_mapping_file.exists():
            self._load_json_mapping(default_mapping_file)
        else:
            logger.warning(f"Default cell type mapping file not found at {default_mapping_file}")

    def map_celltype(self, celltype_name):
        """
        Map a cell type name to a Cell Ontology ID.
        Returns (cl_term_id, cl_term_name, confidence) tuple.
        """
        # Handle edge cases
        if pd.isna(celltype_name) or not celltype_name or str(celltype_name).strip() == '':
            return '', '', 'none'

        # Normalize cell type name (lowercase, strip)
        celltype_lower = str(celltype_name).strip().lower()

        # Check exact lowercase match
        if celltype_lower in self.celltype_mappings:
            return (self.celltype_mappings[celltype_lower],
                    self.celltype_names.get(celltype_lower, ''),
                    self.mapping_confidence.get(celltype_lower, 'medium'))

        # If not found directly, mark as unmapped
        self.unmapped_celltypes.add(celltype_name)  # Store original case
        return '', '', 'none'

    def map_celltypes_in_adata(self, adata, celltype_field='cell_type'):
        """Add cell type ontology mappings to an AnnData object based on loaded mappings."""
        if celltype_field not in adata.obs.columns:
            logger.warning(f"Cell type field '{celltype_field}' not found in adata.obs. Skipping cell type ontology mapping.")
            # Ensure columns exist even if mapping is skipped
            if 'cell_type_ontology_term_id' not in adata.obs:
                adata.obs['cell_type_ontology_term_id'] = ''
            if 'cell_type_ontology_confidence' not in adata.obs:
                adata.obs['cell_type_ontology_confidence'] = ''
            return adata

        logger.info(f"Applying cell type ontology mappings using field '{celltype_field}'")

        # Initialize columns if they don't exist, ensuring string type
        if 'cell_type_ontology_term_id' not in adata.obs.columns:
            adata.obs['cell_type_ontology_term_id'] = pd.Series(index=adata.obs.index, dtype='str')
        else:
             adata.obs['cell_type_ontology_term_id'] = adata.obs['cell_type_ontology_term_id'].astype(str).fillna('')

        if 'cell_type_ontology_confidence' not in adata.obs.columns:
            adata.obs['cell_type_ontology_confidence'] = pd.Series(index=adata.obs.index, dtype='str')
        else:
             adata.obs['cell_type_ontology_confidence'] = adata.obs['cell_type_ontology_confidence'].astype(str).fillna('none')

        # --- Efficient Mapping using pandas map ---
        # Create mapping series from dictionaries
        celltype_map_series = pd.Series(self.celltype_mappings)
        confidence_map_series = pd.Series(self.mapping_confidence)

        # Normalize the cell type column for mapping (lowercase, strip)
        normalized_celltype_col = adata.obs[celltype_field].astype(str).str.strip().str.lower()

        # Apply mapping
        adata.obs['cell_type_ontology_term_id'] = normalized_celltype_col.map(celltype_map_series).fillna('')
        adata.obs['cell_type_ontology_confidence'] = normalized_celltype_col.map(confidence_map_series).fillna('none')

        # Identify unmapped cell types (where ontology is empty but original cell type was not)
        unmapped_mask = (adata.obs['cell_type_ontology_term_id'] == '') & (adata.obs[celltype_field].astype(str).str.strip() != '') & (adata.obs[celltype_field].notna())
        self.unmapped_celltypes.update(adata.obs.loc[unmapped_mask, celltype_field].unique())

        # Log stats
        total_samples = adata.n_obs
        mapped_samples = (adata.obs['cell_type_ontology_term_id'] != '').sum()
        mapping_percentage = (mapped_samples / total_samples) * 100 if total_samples > 0 else 0
        confidence_counts = adata.obs['cell_type_ontology_confidence'].value_counts().to_dict()

        logger.info(f"Cell type ontology mapping complete: {mapped_samples}/{total_samples} ({mapping_percentage:.1f}%) mapped.")
        logger.info(f"Cell type mapping confidence counts: {confidence_counts}")

        if self.unmapped_celltypes:
            # Log only a sample if there are too many
            sample_unmapped = list(self.unmapped_celltypes)[:20]
            logger.warning(f"Found {len(self.unmapped_celltypes)} unmapped cell types. Examples: {sample_unmapped}")

        return adata

# --- Main Enhancement Function ---

def enhance_standardized_metadata(dataset_name, adata, output_file=None, metadata_dir=None):
    """
    Enhance metadata standardization by prioritizing dataset-specific JSON files.

    Args:
        dataset_name: Name of the dataset (e.g., 'encode', 'gtex')
        adata: AnnData object with dataset loaded from Stage 1 output.
        output_file: Optional file path to save the enhanced standardized data.
        metadata_dir: Directory containing dataset-specific metadata JSON files.

    Returns:
        adata with enhanced standardized metadata, or None if error occurs.
    """
    logger.info(f"Enhancing metadata for {dataset_name} dataset using JSON configurations.")

    # --- 1. Load and Apply Dataset-Specific JSON Metadata ---
    dataset_specific_metadata = None
    if metadata_dir:
        try:
            # Validate the function exists before calling
            if callable(load_dataset_specific_metadata):
                 dataset_specific_metadata = load_dataset_specific_metadata(metadata_dir, dataset_name)
            else:
                 logger.error("load_dataset_specific_metadata function not found in rnaseq_utils.")

            if dataset_specific_metadata:
                logger.info(f"Applying metadata from {dataset_name}_metadata.json")
                # Validate the function exists before calling
                if callable(apply_dataset_specific_metadata):
                    adata = apply_dataset_specific_metadata(adata, dataset_specific_metadata)
                else:
                    logger.error("apply_dataset_specific_metadata function not found in rnaseq_utils.")
            else:
                logger.warning(f"No specific metadata JSON found for {dataset_name}. Will rely on defaults and Stage 1 metadata.")
        except ValueError as ve:
             logger.error(f"Validation error loading metadata for {dataset_name}: {ve}. Proceeding without it.")
        except Exception as e:
             logger.error(f"Could not load or apply metadata JSON for {dataset_name}: {e}")


    # --- 2. Ensure Core Fields Exist and Set Defaults (if not set by JSON) ---
    logger.info("Ensuring core fields and harmonized versions are set...")

    # Harmonization versions (critical for validation)
    if 'harmonized_gencode_version' not in adata.uns or adata.uns.get('harmonized_gencode_version') in [None, 'None', 'unknown']:
        adata.uns['harmonized_gencode_version'] = '24' # Default
        logger.warning(f"Setting default harmonized_gencode_version to '24' for {dataset_name}")
    else: # Ensure format consistency (remove 'v')
        adata.uns['harmonized_gencode_version'] = str(adata.uns['harmonized_gencode_version']).replace('v','')

    if 'harmonized_reference_genome' not in adata.uns or adata.uns.get('harmonized_reference_genome') in [None, 'None', 'unknown']:
        adata.uns['harmonized_reference_genome'] = 'hg38' # Default
        logger.warning(f"Setting default harmonized_reference_genome to 'hg38' for {dataset_name}")
    else: # Allow GRCh38 as well
        if adata.uns['harmonized_reference_genome'] == 'GRCh38':
             adata.uns['harmonized_reference_genome'] = 'hg38' # Standardize to hg38 representation

    # Core obs fields
    if 'dataset' not in adata.obs.columns:
        adata.obs['dataset'] = dataset_name
    adata.obs['dataset'] = adata.obs['dataset'].astype('category') # Ensure categorical

    if 'sample_id' not in adata.obs.columns and adata.obs.index.name != 'sample_id':
        adata.obs['sample_id'] = adata.obs.index.astype(str)
    elif 'sample_id' in adata.obs.columns:
         adata.obs['sample_id'] = adata.obs['sample_id'].astype(str)


    # Basic biological fields (ensure they exist, even if empty or 'unknown')
    for field, default_val, dtype in [
        ('species', 'unknown', 'category'),
        ('species_ontology', '', 'str'),
        ('tissue', 'unknown', 'category'),
        ('tissue_ontology', '', 'str'), # Will be populated by mapper
        ('sex', 'unknown', 'category'),
        ('age', '', 'str'), # Keep age flexible, ontology is separate
        ('developmental_stage_ontology', '', 'str'), # Provide a default if needed
        ('data_type', 'unknown', 'category'),
        ('assay_ontology', '', 'str'),
        ('expression_unit', 'unknown', 'category') # Added expression unit
    ]:
        if field not in adata.obs.columns:
            logger.debug(f"Adding missing obs column '{field}' with default '{default_val}'")
            adata.obs[field] = default_val
        # Ensure correct dtype and handle NaNs appropriately
        if dtype == 'category':
            adata.obs[field] = adata.obs[field].astype('category')
        elif dtype == 'str':
             adata.obs[field] = adata.obs[field].astype(str).fillna(default_val if default_val != 'unknown' else '')


    # --- 3. Apply Tissue Ontology Mapping ---
    if 'tissue' in adata.obs.columns:
        logger.info(f"Applying tissue ontology mapping for {dataset_name}")
        # Assuming default mapping dir contains tissue_to_uberon.json
        tissue_mapper = TissueOntologyMapper(default_mapping_dir=metadata_dir or DEFAULT_METADATA_DIR)
        adata = tissue_mapper.map_tissues_in_adata(adata, tissue_field='tissue')
    else:
        logger.warning(f"'tissue' column not found in {dataset_name}. Skipping tissue ontology mapping.")
        if 'tissue_ontology' not in adata.obs: adata.obs['tissue_ontology'] = ''
        if 'tissue_ontology_confidence' not in adata.obs: adata.obs['tissue_ontology_confidence'] = ''

    # --- 3.5. Apply Cell Type Ontology Mapping ---
    # Check for cell_type field (for ENCODE cell lines, MAGE LCLs, GTEx single-cell)
    if 'cell_type' in adata.obs.columns:
        logger.info(f"Applying cell type ontology mapping for {dataset_name}")
        # Load celltype_to_cl.json mapping
        celltype_mapper = CellTypeOntologyMapper(default_mapping_dir=metadata_dir or DEFAULT_METADATA_DIR)
        adata = celltype_mapper.map_celltypes_in_adata(adata, celltype_field='cell_type')
    else:
        logger.info(f"'cell_type' column not found in {dataset_name}. Skipping cell type ontology mapping.")
        if 'cell_type_ontology_term_id' not in adata.obs: adata.obs['cell_type_ontology_term_id'] = ''
        if 'cell_type_ontology_confidence' not in adata.obs: adata.obs['cell_type_ontology_confidence'] = ''


    # --- 4. Final Checks and Logging ---
    logger.info(f"Final metadata checks for {dataset_name}")

    # Ensure required uns fields are strings for saving
    for key in ['harmonized_gencode_version', 'harmonized_reference_genome',
                'original_gencode_version', 'original_reference_genome',
                'gencode_mapping_notes', 'genome_mapping_notes']:
        if key in adata.uns:
            adata.uns[key] = str(adata.uns[key]) if adata.uns[key] is not None else None

    # Log final state
    log_original_gencode = adata.uns.get('original_gencode_version', 'unknown')
    log_harmonized_gencode = adata.uns.get('harmonized_gencode_version', 'unknown')
    log_original_genome = adata.uns.get('original_reference_genome', 'unknown')
    log_harmonized_genome = adata.uns.get('harmonized_reference_genome', 'unknown')
    protocol_type = adata.uns.get('rna_seq_protocol', {}).get('protocol_type', 'unknown')
    protocol_confidence = adata.uns.get('rna_seq_protocol', {}).get('protocol_confidence', 'low')

    logger.info(f"Metadata enhancement complete for {dataset_name}")
    logger.info(f"  Samples: {adata.n_obs}")
    logger.info(f"  Genes: {adata.n_vars}")
    logger.info(f"  Reference Genome: {log_original_genome} (harmonized to {log_harmonized_genome})")
    logger.info(f"  GENCODE Version: v{log_original_gencode} (harmonized to v{log_harmonized_gencode})")
    logger.info(f"  Protocol Type: {protocol_type} (confidence: {protocol_confidence})")

    # --- 5. Save Standardized Data (if output file specified) ---
    if output_file:
        logger.info(f"Attempting to save standardized {dataset_name} dataset to {output_file}")

        # Ensure necessary columns added by TissueOntologyMapper are strings
        for col in ['tissue_ontology', 'tissue_ontology_confidence']:
            if col in adata.obs.columns:
                 adata.obs[col] = adata.obs[col].astype(str)

        try:
            # Use the robust saving logic including serialization of uns
            from standardize_datasets import save_anndata # Import save function locally
            if callable(save_anndata):
                if save_anndata(adata, output_file):
                    logger.info(f"Successfully saved enhanced dataset to {output_file}")
                else:
                    logger.error(f"Failed to save enhanced dataset {output_file}")
                    return None # Indicate failure
            else:
                logger.error("save_anndata function not found in standardize_datasets.py. Cannot save.")
                return None

        except ImportError:
             logger.error("Could not import save_anndata from standardize_datasets.py. Ensure it's available.")
             return None
        except Exception as e:
            logger.error(f"Critical error during final save of {dataset_name}: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None # Indicate failure

    return adata # Return the enhanced AnnData object

# --- Keep process_all_datasets and main function ---
# (No significant changes needed here, they just call enhance_standardized_metadata)

def process_all_datasets(data_dir, output_dir=None, metadata_dir=None):
    """
    Process all datasets in a directory using enhance_standardized_metadata.

    Args:
        data_dir: Directory containing datasets (expects *_standardized_v1.h5ad)
        output_dir: Optional directory to save standardized datasets (as *_standardized_v2.h5ad)
        metadata_dir: Directory containing dataset-specific metadata JSON files
    """
    data_dir = Path(data_dir)

    if not output_dir:
        output_dir = data_dir
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)

    # Prioritize v1 files from Stage 1
    h5ad_files = list(data_dir.glob("*_standardized_v1.h5ad"))

    if not h5ad_files:
        logger.warning(f"No *_standardized_v1.h5ad files found in {data_dir}. Looking for any *.h5ad...")
        h5ad_files = list(data_dir.glob("*.h5ad"))
        # Avoid processing v2 files if they exist
        h5ad_files = [f for f in h5ad_files if '_standardized_v2' not in f.name]

    if not h5ad_files:
        logger.error(f"No suitable h5ad files found in {data_dir} to process.")
        return

    logger.info(f"Found {len(h5ad_files)} h5ad files to process in Stage 2")

    processed_count = 0
    failed_count = 0
    for h5ad_file in h5ad_files:
        # Extract dataset name (e.g., 'encode' from 'encode_standardized_v1.h5ad')
        dataset_name = h5ad_file.stem.split('_')[0].lower()

        try:
            logger.info(f"--- Processing dataset: {dataset_name} ---")
            logger.info(f"Loading {dataset_name} dataset from {h5ad_file}")
            adata = sc.read_h5ad(h5ad_file)

            # Define output file path for v2
            output_file = output_dir / f"{dataset_name}_standardized_v2.h5ad"

            # Enhance metadata standardization
            enhanced_adata = enhance_standardized_metadata(
                dataset_name,
                adata,
                output_file=output_file, # Pass output path for saving
                metadata_dir=metadata_dir
            )

            if enhanced_adata is not None:
                 processed_count += 1
            else:
                 failed_count += 1
                 logger.error(f"Failed to enhance metadata for {dataset_name}")


        except Exception as e:
            failed_count += 1
            logger.error(f"Error processing {dataset_name} from {h5ad_file}: {e}")
            import traceback
            logger.error(traceback.format_exc())

    logger.info(f"--- Stage 2 Summary ---")
    logger.info(f"Successfully processed: {processed_count} datasets")
    logger.info(f"Failed to process: {failed_count} datasets")
    logger.info("----------------------")


def main():
    parser = argparse.ArgumentParser(description="Standardize metadata for RNA-seq datasets (V2 - JSON Prioritized)")
    parser.add_argument("--data-dir", required=True, help="Directory containing *_standardized_v1.h5ad files")
    parser.add_argument("--output-dir", help="Directory to save *_standardized_v2.h5ad datasets (defaults to data-dir)")
    # Removed mapping-dir as TissueOntologyMapper now uses default path
    parser.add_argument("--metadata-dir", default=str(DEFAULT_METADATA_DIR),
                         help="Directory containing dataset-specific metadata JSON files (e.g., encode_metadata.json)")
    parser.add_argument("--dataset", help="Process only this specific dataset name (e.g., 'encode')")

    args = parser.parse_args()

    output_dir = args.output_dir if args.output_dir else args.data_dir

    if args.dataset:
        # Process a single specified dataset
        dataset_name = args.dataset.lower()
        data_dir = Path(args.data_dir)
        # Look for the v1 file first
        h5ad_file = data_dir / f"{dataset_name}_standardized_v1.h5ad"

        if not h5ad_file.exists():
            logger.warning(f"V1 file not found for {dataset_name}, looking for any {dataset_name}*.h5ad...")
            found_files = list(data_dir.glob(f"{dataset_name}*.h5ad"))
            # Exclude v2 files specifically
            found_files = [f for f in found_files if '_standardized_v2' not in f.name]
            if not found_files:
                 logger.error(f"No suitable h5ad file found for dataset {dataset_name} in {data_dir}")
                 return
            h5ad_file = found_files[0] # Process the first one found

        logger.info(f"Processing single dataset: {dataset_name} from {h5ad_file}")
        try:
            adata = sc.read_h5ad(h5ad_file)
            output_file_path = Path(output_dir) / f"{dataset_name}_standardized_v2.h5ad"
            enhance_standardized_metadata(dataset_name, adata, output_file=output_file_path, metadata_dir=args.metadata_dir)
        except Exception as e:
            logger.error(f"Error processing single dataset {dataset_name}: {e}")
            import traceback
            logger.error(traceback.format_exc())
    else:
        # Process all datasets found in the directory
        process_all_datasets(args.data_dir, output_dir, args.metadata_dir)

if __name__ == "__main__":
    main()