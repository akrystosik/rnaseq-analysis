#!/usr/bin/env python3
"""
Apply Dataset Metadata from JSON Files

This script applies dataset-specific metadata from JSON configuration files
to standardized RNA-seq datasets.
"""

import os
import sys
import json
import logging
import argparse
import scanpy as sc
import pandas as pd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('metadata_application')

def load_dataset_metadata(metadata_dir, dataset_name):
    """
    Load dataset-specific metadata from JSON file.
    
    Args:
        metadata_dir: Directory containing metadata JSON files
        dataset_name: Name of the dataset (e.g., 'encode', 'gtex')
    
    Returns:
        Dictionary with dataset metadata or None if not found
    """
    json_path = os.path.join(metadata_dir, f"{dataset_name}_metadata.json")
    
    if not os.path.exists(json_path):
        logger.warning(f"Metadata file not found: {json_path}")
        return None
    
    try:
        with open(json_path, 'r') as f:
            metadata = json.load(f)
        logger.info(f"Loaded metadata for {dataset_name} from {json_path}")
        return metadata
    except Exception as e:
        logger.error(f"Error loading metadata file {json_path}: {e}")
        return None

def apply_dataset_metadata(adata, metadata):
    """
    Apply dataset-specific metadata to an AnnData object.
    
    Args:
        adata: AnnData object
        metadata: Dictionary with dataset metadata
    
    Returns:
        Modified AnnData object
    """
    if metadata is None:
        logger.warning("No metadata provided, nothing to apply")
        return adata
    
    # Apply dataset-level metadata from JSON
    if 'dataset_info' in metadata:
        # Update uns with dataset info
        for key, value in metadata['dataset_info'].items():
            adata.uns[key] = value
            logger.info(f"Updated uns[{key}]")
    
    # Apply obs columns if defined
    if 'obs_columns' in metadata:
        for col, value in metadata['obs_columns'].items():
            # If value is a simple string/number, apply to all cells
            if not isinstance(value, dict):
                adata.obs[col] = value
                logger.info(f"Updated obs column {col} with value {value}")
            else:
                # If value is a dict, interpret as field mappings
                logger.info(f"Applying mappings for obs column {col}")
    
    # Apply ontology fields if defined
    if 'ontology_mappings' in metadata:
        for field, mapping in metadata['ontology_mappings'].items():
            source_field = mapping.get('source_field', field)
            target_field = mapping.get('target_field', f"{field}_ontology")
            
            if source_field in adata.obs.columns:
                logger.info(f"Applying ontology mapping for {source_field} -> {target_field}")
                # Create or update ontology field
                adata.obs[target_field] = adata.obs[source_field].map(
                    lambda x: mapping.get('mappings', {}).get(str(x), '')
                )
            else:
                logger.warning(f"Source field {source_field} not found in obs columns")
    
    # Special handling for common required fields
    required_fields = {
        'dataset': metadata.get('dataset_name', os.path.basename(metadata.get('dataset_info', {}).get('source', ''))),
        'data_type': metadata.get('dataset_info', {}).get('data_type', 'RNA-seq'),
        'expression_unit': metadata.get('dataset_info', {}).get('expression_unit', 'TPM'),
        'species': metadata.get('dataset_info', {}).get('species', 'Homo sapiens'),
        'species_ontology': metadata.get('dataset_info', {}).get('species_ontology', 'NCBITaxon:9606')
    }
    
    # Add missing required fields
    for field, value in required_fields.items():
        if field not in adata.obs.columns:
            adata.obs[field] = value
            logger.info(f"Added missing required field {field} with value {value}")
    
    # Handle assay ontology
    if 'assay_ontology' not in adata.obs.columns:
        assay_value = metadata.get('dataset_info', {}).get('assay_ontology', 'EFO:0009922')  # Default to RNA-seq
        adata.obs['assay_ontology'] = assay_value
        logger.info(f"Added assay_ontology with value {assay_value}")
    
    # Handle developmental stage ontology
    if 'developmental_stage_ontology' not in adata.obs.columns:
        dev_stage = metadata.get('dataset_info', {}).get('developmental_stage_ontology', 'HsapDv:0000087')  # Human adult
        adata.obs['developmental_stage_ontology'] = dev_stage
        logger.info(f"Added developmental_stage_ontology with value {dev_stage}")
    
    logger.info(f"Successfully applied metadata to dataset")
    return adata

def process_file(input_file, output_file, metadata_dir):
    """
    Process a single AnnData file to apply metadata.
    
    Args:
        input_file: Path to input AnnData file
        output_file: Path to output AnnData file
        metadata_dir: Directory containing metadata JSON files
    
    Returns:
        Path to the output file
    """
    logger.info(f"Processing {input_file}")
    
    # Extract dataset name from filename
    filename = os.path.basename(input_file)
    dataset_name = filename.split('_')[0].lower()
    
    # Load AnnData
    adata = sc.read_h5ad(input_file)
    logger.info(f"Loaded AnnData with {adata.n_obs} observations and {adata.n_vars} variables")
    
    # Load and apply metadata
    metadata = load_dataset_metadata(metadata_dir, dataset_name)
    if metadata:
        adata = apply_dataset_metadata(adata, metadata)
    
    # Save modified AnnData
    logger.info(f"Saving to {output_file}")
    adata.write_h5ad(output_file)
    logger.info(f"Saved updated AnnData to {output_file}")
    
    return output_file

def process_directory(input_dir, output_dir, metadata_dir, file_pattern='*.h5ad'):
    """
    Process all AnnData files in a directory.
    
    Args:
        input_dir: Path to input directory
        output_dir: Path to output directory
        metadata_dir: Directory containing metadata JSON files
        file_pattern: Pattern to match AnnData files
    
    Returns:
        Number of files processed
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all matching files
    import glob
    h5ad_files = glob.glob(os.path.join(input_dir, file_pattern))
    
    if not h5ad_files:
        logger.warning(f"No files found matching {file_pattern} in {input_dir}")
        return 0
    
    logger.info(f"Found {len(h5ad_files)} files to process")
    
    # Process each file
    processed_files = 0
    for input_file in h5ad_files:
        output_file = os.path.join(output_dir, os.path.basename(input_file))
        try:
            process_file(input_file, output_file, metadata_dir)
            processed_files += 1
        except Exception as e:
            logger.error(f"Error processing {input_file}: {e}")
            import traceback
            logger.error(traceback.format_exc())
    
    logger.info(f"Processed {processed_files} files")
    return processed_files

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Apply dataset metadata from JSON files')
    parser.add_argument('--input', required=True, help='Input AnnData file or directory')
    parser.add_argument('--output', required=True, help='Output AnnData file or directory')
    parser.add_argument('--metadata-dir', required=True, help='Directory containing metadata JSON files')
    parser.add_argument('--pattern', default='*.h5ad', help='File pattern when processing directories')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Check if input is a file or directory
    if os.path.isfile(args.input):
        # Process single file
        process_file(args.input, args.output, args.metadata_dir)
    elif os.path.isdir(args.input):
        # Process directory
        process_directory(args.input, args.output, args.metadata_dir, args.pattern)
    else:
        logger.error(f"Input not found: {args.input}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
