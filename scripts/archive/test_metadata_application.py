#!/usr/bin/env python3
"""
Test applying metadata from JSON files to AnnData objects
"""
import scanpy as sc
import pandas as pd
import numpy as np
import logging
import sys
import os
import json

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('metadata_test')

def create_test_anndata(name="adni"):
    """Create a test AnnData object"""
    # Create a simple AnnData object
    x = np.random.rand(5, 10)  # 5 cells, 10 genes
    obs = pd.DataFrame(index=[f'cell_{i}' for i in range(5)])
    var = pd.DataFrame(index=[f'gene_{i}' for i in range(10)])
    
    # Add some basic gene information
    var['gene_id'] = [f'ENSG{i:08d}' for i in range(10)]
    var['ensembl_id'] = [f'ENSG{i:08d}' if i < 8 else f'ENTREZ:{i}' for i in range(10)]
    
    # Create AnnData
    adata = sc.AnnData(X=x, obs=obs, var=var)
    
    # Add minimal dataset info
    adata.uns['dataset'] = name
    
    return adata

def load_metadata_json(metadata_dir, dataset_name):
    """Load dataset metadata from JSON file"""
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

def apply_metadata(adata, metadata):
    """Apply metadata from JSON to AnnData object"""
    if metadata is None:
        logger.warning("No metadata provided, nothing to apply")
        return adata
    
    # Apply dataset-level metadata from JSON
    if 'dataset_info' in metadata:
        # Update uns with dataset info
        for key, value in metadata['dataset_info'].items():
            adata.uns[key] = value
            logger.info(f"Updated uns[{key}] = {value}")
    
    # Standard fields that should be in obs
    standard_fields = {
        'dataset': metadata.get('dataset_name', adata.uns.get('dataset', '')),
        'data_type': metadata.get('dataset_info', {}).get('data_type', 'RNA-seq'),
        'expression_unit': metadata.get('dataset_info', {}).get('expression_unit', 'TPM'),
        'species': metadata.get('dataset_info', {}).get('species', 'Homo sapiens'),
        'species_ontology': metadata.get('dataset_info', {}).get('species_ontology', 'NCBITaxon:9606'),
        'tissue': metadata.get('dataset_info', {}).get('tissue', 'blood')
    }
    
    # Add standard fields to obs
    for field, value in standard_fields.items():
        adata.obs[field] = value
        logger.info(f"Added obs[{field}] = {value}")
    
    return adata

def main():
    # Set paths
    metadata_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json'
    test_dir = '/tmp/anndata_metadata_test'
    os.makedirs(test_dir, exist_ok=True)
    
    # Test each dataset
    for dataset_name in ['adni', 'encode', 'gtex', 'mage']:
        logger.info(f"Testing metadata application for {dataset_name}")
        
        # Create test AnnData
        adata = create_test_anndata(dataset_name)
        logger.info(f"Created test AnnData for {dataset_name} with shape {adata.shape}")
        
        # Load metadata
        metadata = load_metadata_json(metadata_dir, dataset_name)
        if metadata:
            logger.info(f"Loaded metadata for {dataset_name}")
            
            # Check what's in the metadata
            if 'dataset_info' in metadata:
                logger.info(f"dataset_info keys: {list(metadata['dataset_info'].keys())}")
            else:
                logger.warning(f"No dataset_info found in {dataset_name} metadata")
            
            # Apply metadata
            adata = apply_metadata(adata, metadata)
            
            # Save test AnnData
            test_file = os.path.join(test_dir, f'{dataset_name}_test_metadata.h5ad')
            adata.write_h5ad(test_file)
            logger.info(f"Saved {dataset_name} test AnnData to {test_file}")
            
            # Check for required fields
            required_fields = ['dataset', 'sample_id', 'data_type', 'expression_unit', 'tissue', 'species']
            missing_fields = [field for field in required_fields if field not in adata.obs.columns]
            
            if missing_fields:
                logger.warning(f"Missing required fields in {dataset_name}: {missing_fields}")
            else:
                logger.info(f"All required fields present in {dataset_name}")
        else:
            logger.warning(f"No metadata found for {dataset_name}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
