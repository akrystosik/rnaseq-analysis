#!/usr/bin/env python3
"""
Test simplified metadata approach for AnnData objects.
Extract only the essential metadata we need and ensure it's all string-based.
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
logger = logging.getLogger('simplified_metadata_test')

def extract_essential_metadata(metadata):
    """
    Extract only the essential metadata fields we need for standardization.
    Convert everything to simple string values.
    """
    # Essential fields we want to preserve
    essential = {
        'data_type': str(metadata.get('data_type', '')),
        'platform': str(metadata.get('platform', '')),
        'reference_genome': str(metadata.get('reference_genome', '')),
        'harmonized_reference_genome': str(metadata.get('harmonized_reference_genome', 'hg38')),
        'gencode_version': str(metadata.get('gencode_version', '')),
        'harmonized_gencode_version': str(metadata.get('harmonized_gencode_version', 'v24')),
        'assay_ontology': str(metadata.get('assay_ontology', '')),
    }
    
    # Extract simplified protocol info
    if 'rna_seq_protocol' in metadata and isinstance(metadata['rna_seq_protocol'], dict):
        protocol = metadata['rna_seq_protocol']
        essential['protocol_type'] = str(protocol.get('protocol_type', ''))
        essential['selection_method'] = str(protocol.get('selection_method', ''))
    
    # Extract simplified tissue mapping for cell lines if present
    if 'cell_type_info' in metadata and isinstance(metadata['cell_type_info'], dict):
        if 'cell_line_mapping' in metadata['cell_type_info']:
            # Create a simplified format, not nested
            cell_line_mapping = {}
            for cell_line, info in metadata['cell_type_info']['cell_line_mapping'].items():
                cell_line_mapping[str(cell_line)] = {
                    'anatomical_entity_id': str(info.get('anatomical_entity_id', '')),
                    'cell_type_id': str(info.get('cell_type_id', ''))
                }
            essential['cell_line_mapping'] = cell_line_mapping
    
    # For ENCODE, extract the simplified cell line tissue information
    if ('dataset_info' in metadata and isinstance(metadata['dataset_info'], dict) and
            'cell_lines' in metadata['dataset_info']):
        cell_lines = {}
        for cell_line, info in metadata['dataset_info']['cell_lines'].items():
            cell_lines[str(cell_line)] = {
                'tissue': str(info.get('tissue', '')),
                'cell_type': str(info.get('cell_type', '')),
                'organism': str(info.get('organism', ''))
            }
        essential['cell_lines'] = cell_lines
    
    # Extract observation column mappings if present
    if 'obs_columns' in metadata and isinstance(metadata['obs_columns'], dict):
        obs_mappings = {}
        for key, value in metadata['obs_columns'].items():
            obs_mappings[str(key)] = str(value)
        essential['obs_mappings'] = obs_mappings
    
    return essential

def test_with_real_metadata():
    """Test with real metadata from the project"""
    metadata_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json"
    
    # Test each metadata file
    for dataset in ['adni', 'encode', 'gtex', 'mage']:
        metadata_file = os.path.join(metadata_dir, f"{dataset}_metadata.json")
        
        if not os.path.exists(metadata_file):
            logger.warning(f"Metadata file not found: {metadata_file}")
            continue
        
        logger.info(f"Testing {dataset} metadata file")
        
        try:
            # Load the metadata
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            
            # Extract essential metadata
            simplified = extract_essential_metadata(metadata)
            
            # Print what we've extracted
            logger.info(f"Extracted {len(simplified)} essential fields for {dataset}")
            logger.info(f"Fields: {', '.join(simplified.keys())}")
            
            # Verify we can serialize to JSON
            json_str = json.dumps(simplified)
            logger.info(f"Successfully serialized to JSON ({len(json_str)} characters)")
            
            # Create AnnData and save
            adata = sc.AnnData(X=np.zeros((2, 3)))
            adata.uns['dataset_metadata'] = simplified
            
            test_file = f"/tmp/{dataset}_simplified_test.h5ad"
            adata.write_h5ad(test_file)
            logger.info(f"Successfully saved {dataset} metadata to {test_file}")
            
            # Load it back
            loaded = sc.read_h5ad(test_file)
            logger.info(f"Successfully loaded {dataset} metadata from {test_file}")
            
            # Check if all fields were preserved
            for key in simplified.keys():
                if key in loaded.uns['dataset_metadata']:
                    logger.info(f"Field '{key}' preserved")
                else:
                    logger.warning(f"Field '{key}' lost")
            
        except Exception as e:
            logger.error(f"Error processing {dataset} metadata: {e}")
            import traceback
            logger.error(traceback.format_exc())

# Run the tests
if __name__ == "__main__":
    logger.info("=== Testing Simplified Metadata Approach ===")
    test_with_real_metadata()
