#!/usr/bin/env python3
"""
Test script to verify cell type ontology mapping functionality.

This script tests:
1. celltype_to_cl.json loading
2. CellTypeOntologyMapper functionality  
3. Integration with standardize_metadata.py
"""

import os
import sys
import json
import pandas as pd
import logging
from pathlib import Path

# Add the pipeline directory to path
sys.path.insert(0, '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2')

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger('test_celltype_mapping')

def test_celltype_mapping_file():
    """Test loading and parsing of celltype_to_cl.json"""
    logger.info("=== Testing celltype_to_cl.json file ===")
    
    mapping_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/celltype_to_cl.json'
    
    if not os.path.exists(mapping_file):
        logger.error(f"Mapping file not found: {mapping_file}")
        return False
    
    try:
        with open(mapping_file, 'r') as f:
            data = json.load(f)
        
        mappings = data.get('mappings', {})
        logger.info(f"Loaded {len(mappings)} cell type mappings")
        
        # Test some key mappings - validate against JSON structure
        test_cases = []
        key_cell_types = ['A549', 'lymphoblastoid cell line', 'K562', 'hepatocyte']
        
        for cell_type in key_cell_types:
            if cell_type in mappings:
                mapping_info = mappings[cell_type]
                if isinstance(mapping_info, dict):
                    expected_id = mapping_info.get('cl_term_id', '')
                    expected_name = mapping_info.get('cl_term_name', '')
                    test_cases.append((cell_type, expected_id, expected_name))
                else:
                    logger.warning(f"Cell type {cell_type} mapping is not a dict in JSON")
            else:
                logger.warning(f"Cell type {cell_type} not found in JSON mappings")
        
        for cell_type, expected_id, expected_name in test_cases:
            if cell_type in mappings:
                mapping_info = mappings[cell_type]
                if isinstance(mapping_info, dict):
                    actual_id = mapping_info.get('cl_term_id', '')
                    actual_name = mapping_info.get('cl_term_name', '')
                    if actual_id == expected_id and actual_name == expected_name:
                        logger.info(f"✓ {cell_type} -> {actual_id} ({actual_name})")
                    else:
                        logger.error(f"✗ {cell_type} mapping mismatch: expected {expected_id}, got {actual_id}")
                        return False
                else:
                    logger.error(f"✗ {cell_type} mapping is not a dict: {mapping_info}")
                    return False
            else:
                logger.error(f"✗ {cell_type} not found in mappings")
                return False
        
        logger.info("✓ celltype_to_cl.json file test passed")
        return True
        
    except Exception as e:
        logger.error(f"Error loading mapping file: {e}")
        return False

def test_celltype_mapper_class():
    """Test the CellTypeOntologyMapper class"""
    logger.info("=== Testing CellTypeOntologyMapper class ===")
    
    try:
        from standardize_metadata import CellTypeOntologyMapper
        
        # Initialize mapper
        metadata_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json'
        mapper = CellTypeOntologyMapper(default_mapping_dir=metadata_dir)
        
        # Test individual mappings
        test_cases = [
            ('A549', 'CL:0000183', 'epithelial cell of lung'),
            ('a549', 'CL:0000183', 'epithelial cell of lung'),  # Case insensitive
            ('lymphoblastoid cell line', 'CL:0000542', 'lymphoblast'),
            ('unknown_cell_type', '', '', 'none')  # Should not map
        ]
        
        for case in test_cases:
            if len(case) == 4:
                cell_type, expected_id, expected_name, expected_conf = case
            else:
                cell_type, expected_id, expected_name = case
                expected_conf = 'medium'  # Default
            
            cl_id, cl_name, confidence = mapper.map_celltype(cell_type)
            
            if cl_id == expected_id and cl_name == expected_name:
                logger.info(f"✓ {cell_type} -> {cl_id} ({cl_name}) [{confidence}]")
            else:
                logger.error(f"✗ {cell_type} mapping failed: expected ({expected_id}, {expected_name}), got ({cl_id}, {cl_name})")
                return False
        
        logger.info("✓ CellTypeOntologyMapper class test passed")
        return True
        
    except ImportError as e:
        logger.error(f"Could not import CellTypeOntologyMapper: {e}")
        return False
    except Exception as e:
        logger.error(f"Error testing CellTypeOntologyMapper: {e}")
        return False

def test_anndata_integration():
    """Test integration with AnnData objects"""
    logger.info("=== Testing AnnData integration ===")
    
    try:
        import scanpy as sc
        import numpy as np
        from standardize_metadata import CellTypeOntologyMapper
        
        # Create a mock AnnData object
        n_obs = 10
        n_vars = 100
        
        # Create mock data
        X = np.random.randn(n_obs, n_vars)
        obs = pd.DataFrame({
            'cell_type': ['A549', 'HepG2', 'lymphoblastoid cell line', 'K562', 'unknown'] * 2,
            'sample_id': [f'sample_{i}' for i in range(n_obs)]
        })
        obs.index = obs['sample_id']
        
        var = pd.DataFrame({
            'gene_id': [f'ENSG{i:011d}' for i in range(n_vars)]
        })
        var.index = var['gene_id']
        
        adata = sc.AnnData(X=X, obs=obs, var=var)
        
        logger.info(f"Created mock AnnData: {adata.n_obs} obs x {adata.n_vars} vars")
        logger.info(f"Cell types: {adata.obs['cell_type'].unique()}")
        
        # Apply cell type mapping
        metadata_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json'
        mapper = CellTypeOntologyMapper(default_mapping_dir=metadata_dir)
        adata_mapped = mapper.map_celltypes_in_adata(adata, celltype_field='cell_type')
        
        # Check results
        if 'cell_type_ontology_term_id' not in adata_mapped.obs.columns:
            logger.error("cell_type_ontology_term_id column not created")
            return False
        
        if 'cell_type_ontology_confidence' not in adata_mapped.obs.columns:
            logger.error("cell_type_ontology_confidence column not created")
            return False
        
        # Check specific mappings
        mapped_count = (adata_mapped.obs['cell_type_ontology_term_id'] != '').sum()
        logger.info(f"Mapped {mapped_count}/{adata_mapped.n_obs} samples to cell type ontology terms")
        
        # Verify specific mappings
        for idx, row in adata_mapped.obs.iterrows():
            cell_type = row['cell_type']
            cl_id = row['cell_type_ontology_term_id']
            confidence = row['cell_type_ontology_confidence']
            
            if cell_type == 'A549' and cl_id == 'CL:0000183':
                logger.info(f"✓ {cell_type} correctly mapped to {cl_id}")
            elif cell_type == 'unknown' and cl_id == '':
                logger.info(f"✓ {cell_type} correctly unmapped")
            elif cell_type in ['A549', 'HepG2', 'lymphoblastoid cell line', 'K562'] and cl_id.startswith('CL:'):
                logger.info(f"✓ {cell_type} mapped to {cl_id}")
        
        logger.info("✓ AnnData integration test passed")
        return True
        
    except Exception as e:
        logger.error(f"Error in AnnData integration test: {e}")
        return False

def main():
    """Run all tests"""
    logger.info("Starting cell type ontology mapping tests...")
    
    results = []
    
    # Test 1: Mapping file
    results.append(test_celltype_mapping_file())
    
    # Test 2: Mapper class
    results.append(test_celltype_mapper_class())
    
    # Test 3: AnnData integration
    results.append(test_anndata_integration())
    
    # Summary
    passed = sum(results)
    total = len(results)
    
    logger.info(f"\n=== Test Summary ===")
    logger.info(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        logger.info("✓ All tests passed! Cell type ontology mapping is working correctly.")
        return 0
    else:
        logger.error("✗ Some tests failed. Check the implementation.")
        return 1

if __name__ == '__main__':
    sys.exit(main())