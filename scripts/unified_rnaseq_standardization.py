#!/usr/bin/env python3
"""
Unified RNA-seq Standardization Pipeline

This script combines the functionality of standardize_datasets.py and standardize_metadata.py
to provide a complete pipeline for standardizing RNA-seq data from multiple sources.

The pipeline processes data in two stages:
1. Initial data conversion (from standardize_datasets.py)
2. Enhanced metadata standardization (from standardize_metadata.py)

Usage:
  python unified_rnaseq_standardization.py --encode-dir /path/to/encode/data \
                                           --gtex-file /path/to/gtex.gct.gz \
                                           --mage-dir /path/to/mage/data \
                                           --adni-dir /path/to/adni/data \
                                           --metadata-dir /path/to/metadata/json \
                                           --output-dir /path/to/output
"""

import os
import sys
import json
import argparse
import logging
import time
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from pathlib import Path
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"rnaseq_standardization_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('unified_pipeline')

# Define paths and constants
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
GENCODE_MAPPING_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gencode_v24_complete_mapping.csv")
DEFAULT_OUTPUT_DIR = BASE_DIR / "standardized_data"

# Import functions from standardize_datasets.py
sys.path.append(str(BASE_DIR / "scripts"))
try:
    from standardize_datasets import (
        # Utility functions
        standardize_ensembl_id, load_gencode_mapping, add_gencode_annotations,
        standardize_metadata as std_ds_standardize_metadata, 
        validate_metadata, prepare_for_anndata, create_standard_anndata, save_anndata,
        
        # Dataset-specific processing functions
        process_encode_data, process_entex_data, process_gtex_data, 
        process_mage_data, process_adni_data,
        
        # Standardization mappings (to be replaced with JSON loading)
        TISSUE_TO_UBERON, ASSAY_TO_EFO, AGE_TO_HSAPDV, SEX_STANDARDIZATION,
        SPECIES_TO_NCBI_TAXON, ETHNICITY_STANDARDIZATION
    )
    logger.info("Successfully imported functions from standardize_datasets.py")
except ImportError as e:
    logger.error(f"Error importing from standardize_datasets.py: {e}")
    logger.error("Make sure the script is in the expected location")
    sys.exit(1)

# Import functions from standardize_metadata.py
try:
    from standardize_metadata import (
        TissueOntologyMapper, infer_rnaseq_protocol, infer_gene_ids_format,
        infer_species, infer_genome_version, load_dataset_specific_metadata,
        apply_dataset_specific_metadata, standardize_metadata,
        validate_dataset_metadata
    )
    logger.info("Successfully imported functions from standardize_metadata.py")
except ImportError as e:
    logger.error(f"Error importing from standardize_metadata.py: {e}")
    logger.error("Make sure the script is in the expected location")
    sys.exit(1)

def load_metadata_mappings(metadata_dir):
    """
    Load metadata mappings from JSON files.
    
    Args:
        metadata_dir: Directory containing JSON metadata files
        
    Returns:
        Dictionary of mappings
    """
    mappings = {}
    metadata_dir = Path(metadata_dir)
    
    # Define the files to load
    metadata_files = {
        'ENCODE_METADATA': 'encode_metadata.json',
        'GTEX_METADATA': 'gtex_metadata.json',
        'MAGE_METADATA': 'mage_metadata.json',
        'ADNI_METADATA': 'adni_metadata.json',
        'ENTEX_METADATA': 'entex_metadata.json',
        # Add other mapping files here
    }
    
    # Load each file
    for key, filename in metadata_files.items():
        file_path = metadata_dir / filename
        if file_path.exists():
            try:
                with open(file_path, 'r') as f:
                    mappings[key] = json.load(f)
                logger.info(f"Loaded metadata from {file_path}")
            except Exception as e:
                logger.error(f"Error loading metadata from {file_path}: {e}")
                mappings[key] = {}
        else:
            logger.warning(f"Metadata file not found: {file_path}")
            mappings[key] = {}
    
    # Load tissue mapping
    tissue_mapping_file = metadata_dir / 'tissue_to_uberon.json'
    if tissue_mapping_file.exists():
        try:
            with open(tissue_mapping_file, 'r') as f:
                mappings['TISSUE_TO_UBERON'] = json.load(f)
            logger.info(f"Loaded tissue mapping from {tissue_mapping_file}")
        except Exception as e:
            logger.error(f"Error loading tissue mapping: {e}")
            # Fall back to hardcoded mapping
            mappings['TISSUE_TO_UBERON'] = TISSUE_TO_UBERON
    else:
        logger.warning(f"Tissue mapping file not found, using default mapping")
        mappings['TISSUE_TO_UBERON'] = TISSUE_TO_UBERON
    
    return mappings

def stage1_process_dataset(dataset_type, args, mappings):
    """
    Stage 1: Initial data conversion using functions from standardize_datasets.py
    
    Args:
        dataset_type: Type of dataset to process (encode, gtex, mage, adni, entex)
        args: Command-line arguments
        mappings: Dictionary of metadata mappings from JSON files
        
    Returns:
        Path to the output file if successful, None otherwise
    """
    try:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Define stage 1 output filename
        output_file = output_dir / f"{dataset_type}_standardized_v1.h5ad"
        
        logger.info(f"=== Stage 1: Processing {dataset_type.upper()} dataset ===")
        
        # Process based on dataset type
        if dataset_type == 'encode':
            if not args.encode_dir:
                logger.warning("No ENCODE directory specified, skipping")
                return None
                
            result = process_encode_data(
                args.encode_dir, 
                args.encode_entex_dir, 
                output_file,
                args.entex_metadata_file
            )
            
        elif dataset_type == 'entex':
            if not args.entex_metadata_file:
                logger.warning("No ENTEx metadata file specified, skipping")
                return None
                
            result = process_entex_data(
                args.entex_metadata_file, 
                output_file
            )
            
        elif dataset_type == 'gtex':
            if not args.gtex_file:
                logger.warning("No GTEx file specified, skipping")
                return None
                
            result = process_gtex_data(
                args.gtex_file, 
                output_file
            )
            
        elif dataset_type == 'mage':
            if not args.mage_dir:
                logger.warning("No MAGE directory specified, skipping")
                return None
                
            result = process_mage_data(
                args.mage_dir, 
                output_file
            )
            
        elif dataset_type == 'adni':
            if not args.adni_dir:
                logger.warning("No ADNI directory specified, skipping")
                return None
                
            result = process_adni_data(
                args.adni_dir, 
                output_file
            )
            
        else:
            logger.error(f"Unknown dataset type: {dataset_type}")
            return None
        
        if result is not None:
            logger.info(f"Stage 1 processing for {dataset_type} completed successfully")
            logger.info(f"Output saved to {output_file}")
            return output_file
        else:
            logger.error(f"Stage 1 processing for {dataset_type} failed")
            return None
            
    except Exception as e:
        logger.error(f"Error in Stage 1 processing for {dataset_type}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def stage2_enhance_metadata(dataset_type, input_file, args, mappings):
    """
    Stage 2: Enhanced metadata standardization using functions from standardize_metadata.py
    
    Args:
        dataset_type: Type of dataset (encode, gtex, mage, adni, entex)
        input_file: Path to the input file from Stage 1
        args: Command-line arguments
        mappings: Dictionary of metadata mappings from JSON files
        
    Returns:
        Path to the output file if successful, None otherwise
    """
    try:
        output_dir = Path(args.output_dir)
        input_file = Path(input_file)
        
        # Define stage 2 output filename
        output_file = output_dir / f"{dataset_type}_standardized_v2.h5ad"
        
        logger.info(f"=== Stage 2: Enhancing metadata for {dataset_type.upper()} dataset ===")
        
        # Load the dataset
        logger.info(f"Loading {dataset_type} dataset from {input_file}")
        adata = sc.read_h5ad(input_file)
        
        # Load dataset-specific metadata
        metadata_key = f"{dataset_type.upper()}_METADATA"
        dataset_metadata = mappings.get(metadata_key, {})
        
        if not dataset_metadata:
            logger.warning(f"No metadata found for {dataset_type} in mappings")
            # Try to load from file
            metadata_file = Path(args.metadata_dir) / f"{dataset_type}_metadata.json"
            if metadata_file.exists():
                try:
                    with open(metadata_file, 'r') as f:
                        dataset_metadata = json.load(f)
                    logger.info(f"Loaded metadata from {metadata_file}")
                except Exception as e:
                    logger.error(f"Error loading metadata from {metadata_file}: {e}")
        
        # Apply dataset-specific metadata from JSON
        if dataset_metadata:
            adata = apply_dataset_specific_metadata(adata, dataset_metadata)
            logger.info(f"Applied dataset-specific metadata for {dataset_type}")
        
        # Create tissue ontology mapper
        mapping_dir = Path(args.metadata_dir) if args.metadata_dir else None
        tissue_mapper = TissueOntologyMapper()
        
        # Add tissue ontology mappings if 'tissue' is in obs
        if 'tissue' in adata.obs.columns:
            adata = tissue_mapper.map_tissues_in_adata(adata, tissue_field='tissue')
        
        # Enhance metadata using functions from standardize_metadata.py
        adata = standardize_metadata(
            dataset_type, 
            adata, 
            output_file=output_file,
            mapping_dir=args.metadata_dir, 
            metadata_dir=args.metadata_dir
        )
        
        logger.info(f"Stage 2 processing for {dataset_type} completed successfully")
        logger.info(f"Output saved to {output_file}")
        return output_file
        
    except Exception as e:
        logger.error(f"Error in Stage 2 processing for {dataset_type}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def validate_standardized_dataset(dataset_type, file_path, mappings):
    """
    Validate a standardized dataset against requirements.
    
    Args:
        dataset_type: Type of dataset (encode, gtex, mage, adni, entex)
        file_path: Path to the standardized dataset file
        mappings: Dictionary of metadata mappings
        
    Returns:
        Dictionary with validation results
    """
    try:
        logger.info(f"Validating standardized {dataset_type} dataset: {file_path}")
        file_path = Path(file_path)
        
        if not file_path.exists():
            logger.error(f"File not found: {file_path}")
            return {
                'dataset': dataset_type,
                'status': 'error',
                'message': 'File not found',
                'validations': {}
            }
        
        # Load the dataset
        adata = sc.read_h5ad(file_path)
        
        # Initialize validation results
        validation_results = {
            'dataset': dataset_type,
            'status': 'passed',
            'n_samples': adata.n_obs,
            'n_genes': adata.n_vars,
            'validations': {}
        }
        
        # Check harmonized GENCODE version
        gencode_version = adata.uns.get('harmonized_gencode_version', None)
        if gencode_version is not None:
            gencode_version = str(gencode_version).replace('v', '')
            if gencode_version == '24':
                validation_results['validations']['gencode_version'] = {
                    'status': 'passed',
                    'value': gencode_version
                }
            else:
                validation_results['validations']['gencode_version'] = {
                    'status': 'failed',
                    'value': gencode_version,
                    'expected': '24'
                }
                validation_results['status'] = 'failed'
        else:
            validation_results['validations']['gencode_version'] = {
                'status': 'missing',
                'expected': '24'
            }
            validation_results['status'] = 'failed'
        
        # Check harmonized reference genome
        genome_version = adata.uns.get('harmonized_reference_genome', None)
        if genome_version is not None:
            if genome_version in ['hg38', 'GRCh38']:
                validation_results['validations']['reference_genome'] = {
                    'status': 'passed',
                    'value': genome_version
                }
            else:
                validation_results['validations']['reference_genome'] = {
                    'status': 'failed',
                    'value': genome_version,
                    'expected': 'hg38/GRCh38'
                }
                validation_results['status'] = 'failed'
        else:
            validation_results['validations']['reference_genome'] = {
                'status': 'missing',
                'expected': 'hg38/GRCh38'
            }
            validation_results['status'] = 'failed'
        
        # Check observation metadata fields
        metadata_fields = {
            'tissue': {
                'ontology_field': 'tissue_ontology',
                'ontology_prefix': 'UBERON:',
                'importance': 'critical'
            },
            'sex': {
                'values': ['male', 'female', 'unknown'],
                'importance': 'important'
            },
            'species': {
                'ontology_field': 'species_ontology',
                'ontology_prefix': 'NCBITaxon:',
                'importance': 'important'
            },
            'data_type': {
                'values': ['RNA-seq', 'microarray'],
                'importance': 'important'
            },
            'assay_ontology': {
                'ontology_prefix': 'EFO:',
                'importance': 'important'
            }
        }
        
        for field, config in metadata_fields.items():
            if field in adata.obs.columns:
                # Check for missing values
                missing_count = adata.obs[field].isna().sum()
                missing_percentage = (missing_count / adata.n_obs) * 100
                
                # Check for ontology field
                if 'ontology_field' in config and config['ontology_field'] in adata.obs.columns:
                    # Check if ontology values are valid
                    ontology_field = config['ontology_field']
                    ontology_prefix = config.get('ontology_prefix', '')
                    
                    # Count values with correct prefix
                    valid_values = adata.obs[ontology_field].astype(str).str.startswith(ontology_prefix)
                    valid_count = valid_values.sum()
                    valid_percentage = (valid_count / adata.n_obs) * 100
                    
                    validation_results['validations'][field] = {
                        'status': 'passed' if valid_percentage >= 90 else 'warning' if valid_percentage >= 70 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'valid_percentage': float(valid_percentage),
                        'has_ontology': True
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
                
                # Check for enumerated values
                elif 'values' in config:
                    valid_values = adata.obs[field].isin(config['values'])
                    valid_count = valid_values.sum()
                    valid_percentage = (valid_count / adata.n_obs) * 100
                    
                    validation_results['validations'][field] = {
                        'status': 'passed' if valid_percentage >= 90 else 'warning' if valid_percentage >= 70 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'valid_percentage': float(valid_percentage),
                        'has_ontology': False
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
                
                # Simple presence check
                else:
                    validation_results['validations'][field] = {
                        'status': 'passed' if missing_percentage <= 10 else 'warning' if missing_percentage <= 30 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'has_ontology': False
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
            else:
                validation_results['validations'][field] = {
                    'status': 'missing',
                    'importance': config['importance']
                }
                
                if config['importance'] == 'critical':
                    validation_results['status'] = 'failed'
        
        logger.info(f"Validation for {dataset_type} completed: {validation_results['status']}")
        return validation_results
        
    except Exception as e:
        logger.error(f"Error validating {dataset_type} dataset: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {
            'dataset': dataset_type,
            'status': 'error',
            'message': str(e),
            'validations': {}
        }

def generate_summary_report(validation_results, output_dir):
    """
    Generate a summary report of validation results.
    
    Args:
        validation_results: List of validation result dictionaries
        output_dir: Directory to save the report
        
    Returns:
        Path to the report file
    """
    try:
        output_dir = Path(output_dir)
        report_file = output_dir / f"standardization_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        
        # Create summary statistics
        summary = {
            'timestamp': datetime.now().isoformat(),
            'datasets_processed': len(validation_results),
            'datasets_passed': sum(1 for r in validation_results if r['status'] == 'passed'),
            'datasets_warning': sum(1 for r in validation_results if r['status'] == 'warning'),
            'datasets_failed': sum(1 for r in validation_results if r['status'] == 'failed'),
            'datasets_error': sum(1 for r in validation_results if r['status'] == 'error'),
            'total_samples': sum(r.get('n_samples', 0) for r in validation_results),
            'total_genes': sum(r.get('n_genes', 0) for r in validation_results),
            'dataset_results': validation_results
        }
        
        # Write report to file
        with open(report_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Also log summary
        logger.info("=== Standardization Summary ===")
        logger.info(f"Datasets processed: {summary['datasets_processed']}")
        logger.info(f"Datasets passed: {summary['datasets_passed']}")
        logger.info(f"Datasets with warnings: {summary['datasets_warning']}")
        logger.info(f"Datasets failed: {summary['datasets_failed']}")
        logger.info(f"Datasets with errors: {summary['datasets_error']}")
        logger.info(f"Total samples processed: {summary['total_samples']}")
        logger.info(f"Total genes processed: {summary['total_genes']}")
        logger.info(f"Detailed report saved to {report_file}")
        
        return report_file
        
    except Exception as e:
        logger.error(f"Error generating summary report: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Unified RNA-seq Standardization Pipeline')
    
    parser.add_argument('--encode-dir', help='Directory containing ENCODE cell line TPM files')
    parser.add_argument('--encode-entex-dir', help='Directory containing ENCODE ENTEx TPM files')
    parser.add_argument('--entex-metadata-file', help='Path to ENTEx metadata JSON file')
    parser.add_argument('--gtex-file', help='Path to GTEx expression file (GCT format)')
    parser.add_argument('--mage-dir', help='Directory containing MAGE expression files')
    parser.add_argument('--adni-dir', help='Directory containing ADNI sample directories with expression files')
    parser.add_argument('--metadata-dir', default=str(BASE_DIR / "metadata/json"), 
                        help='Directory containing metadata JSON files')
    parser.add_argument('--output-dir', default=str(DEFAULT_OUTPUT_DIR), 
                        help='Output directory for standardized files')
    parser.add_argument('--datasets', nargs='+', 
                        choices=['encode', 'entex', 'gtex', 'mage', 'adni', 'all'],
                        default=['all'], help='Datasets to process')
    parser.add_argument('--skip-stage1', action='store_true', 
                        help='Skip Stage 1 processing (use existing v1 files)')
    parser.add_argument('--skip-stage2', action='store_true', 
                        help='Skip Stage 2 processing (only perform Stage 1)')
    parser.add_argument('--skip-validation', action='store_true', 
                        help='Skip validation of standardized datasets')
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Start time for performance tracking
    start_time = time.time()
    
    # Determine datasets to process
    datasets_to_process = []
    if 'all' in args.datasets:
        datasets_to_process = ['encode', 'entex', 'gtex', 'mage', 'adni']
    else:
        datasets_to_process = args.datasets
    
    logger.info(f"=== Starting Unified RNA-seq Standardization Pipeline ===")
    logger.info(f"Datasets to process: {', '.join(datasets_to_process)}")
    
    # Load metadata mappings from JSON files
    mappings = load_metadata_mappings(args.metadata_dir)
    
    # Track processing results
    stage1_results = {}
    stage2_results = {}
    validation_results = []
    
    # Process each dataset
    for dataset_type in datasets_to_process:
        logger.info(f"Processing {dataset_type} dataset")
        
        # Stage 1: Initial data conversion
        if not args.skip_stage1:
            stage1_output = stage1_process_dataset(dataset_type, args, mappings)
            stage1_results[dataset_type] = stage1_output
        else:
            # Check for existing v1 file
            v1_file = Path(args.output_dir) / f"{dataset_type}_standardized_v1.h5ad"
            if v1_file.exists():
                logger.info(f"Using existing Stage 1 output for {dataset_type}: {v1_file}")
                stage1_results[dataset_type] = v1_file
            else:
                logger.warning(f"Stage 1 output not found for {dataset_type} and --skip-stage1 specified")
                continue
        
        # Stage 2: Enhanced metadata standardization
        if not args.skip_stage2 and stage1_results.get(dataset_type) is not None:
            stage2_output = stage2_enhance_metadata(
                dataset_type, 
                stage1_results[dataset_type], 
                args, 
                mappings
            )
            stage2_results[dataset_type] = stage2_output
        elif args.skip_stage2:
            logger.info(f"Skipping Stage 2 for {dataset_type} (--skip-stage2 specified)")
        elif stage1_results.get(dataset_type) is None:
            logger.warning(f"Skipping Stage 2 for {dataset_type} (Stage 1 failed)")
        
        # Validation
        if not args.skip_validation:
            # Determine which file to validate
            if not args.skip_stage2 and dataset_type in stage2_results and stage2_results[dataset_type] is not None:
                validation_file = stage2_results[dataset_type]
            elif dataset_type in stage1_results and stage1_results[dataset_type] is not None:
                validation_file = stage1_results[dataset_type]
            else:
                logger.warning(f"No output file available for {dataset_type} validation")
                continue
            
            # Validate the dataset
            validation_result = validate_standardized_dataset(
                dataset_type, 
                validation_file, 
                mappings
            )
            validation_results.append(validation_result)
    
    # Generate summary report
    if validation_results:
        report_file = generate_summary_report(validation_results, args.output_dir)
        if report_file:
            logger.info(f"Summary report saved to {report_file}")
    
    # Log total execution time
    total_time = time.time() - start_time
    logger.info(f"=== Unified Pipeline Completed in {total_time:.2f} seconds ===")

if __name__ == '__main__':
    main()