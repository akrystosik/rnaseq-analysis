#!/usr/bin/env python3
"""
Test a more comprehensive deep conversion of metadata to handle complex nested types
"""
import scanpy as sc
import pandas as pd
import numpy as np
import logging
import sys
import os
import json
import collections

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('metadata_deep_test')

def deep_str_conversion(obj):
    """
    Recursively convert all objects to Python primitives, handling many possible types.
    This is a brute-force approach to ensure all objects are serializable.
    """
    if obj is None:
        return None
    elif isinstance(obj, (str, int, float, bool)):
        return obj
    elif isinstance(obj, (np.int64, np.int32, np.float64, np.float32, np.bool_)):
        return obj.item()  # Convert numpy scalars to Python types
    elif isinstance(obj, (list, tuple, set)):
        return [deep_str_conversion(item) for item in obj]
    elif isinstance(obj, dict):
        return {str(k): deep_str_conversion(v) for k, v in obj.items()}
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict('records')
    elif isinstance(obj, pd.Series):
        return obj.to_dict()
    elif isinstance(obj, np.ndarray):
        if np.issubdtype(obj.dtype, np.number):
            return obj.tolist()
        else:
            return [str(x) for x in obj]
    elif hasattr(obj, 'tolist'):
        return obj.tolist()
    else:
        return str(obj)

def create_test_dict():
    """Create a deliberately complex nested dictionary to test conversion"""
    test_dict = {
        'regular_string': 'test',
        'int': 42,
        'float': 3.14,
        'numpy_int': np.int64(100),
        'numpy_float': np.float32(2.718),
        'numpy_array': np.array([1, 2, 3]),
        'mixed_array': np.array([1, 'text', 3.5]),
        'nested_dict': {
            'pandas_series': pd.Series([1, 2, 3], index=['a', 'b', 'c']),
            'object_list': [object(), object()],
        },
        'list_of_dicts': [
            {'name': 'item1', 'value': np.float64(1.23)},
            {'name': 'item2', 'value': np.int32(42)},
        ],
        'complex_object': {'data_repository': [{'name': 'test', 'url': 'https://example.com'}]}
    }
    return test_dict

def save_as_anndata(metadata_dict):
    """Test saving metadata dictionary as part of an AnnData object"""
    # Create a minimal AnnData object
    adata = sc.AnnData(X=np.zeros((2, 3)))
    
    # Add the test dictionary to uns
    adata.uns['test_metadata'] = metadata_dict
    
    # Try to save
    test_file = '/tmp/test_metadata_conversion.h5ad'
    try:
        adata.write_h5ad(test_file)
        logger.info(f"Successfully saved AnnData with complex metadata to {test_file}")
        
        # Try to load it
        loaded_adata = sc.read_h5ad(test_file)
        logger.info("Successfully loaded AnnData with converted metadata")
        
        # Compare keys before and after
        original_keys = set(metadata_dict.keys())
        loaded_keys = set(loaded_adata.uns['test_metadata'].keys())
        logger.info(f"Original keys: {original_keys}")
        logger.info(f"Loaded keys: {loaded_keys}")
        
        return True
    except Exception as e:
        logger.error(f"Error saving/loading AnnData: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def test_metadata_conversion():
    """Test the metadata conversion process"""
    # Create test data
    test_dict = create_test_dict()
    logger.info("Created test dictionary with complex nested types")
    
    # Try to save as JSON (will fail with default conversion)
    try:
        json.dumps(test_dict)
        logger.info("Test dictionary can be serialized directly to JSON (unexpected)")
    except Exception as e:
        logger.info(f"As expected, test dictionary cannot be serialized directly to JSON: {e}")
        
        # Convert and try again
        converted_dict = deep_str_conversion(test_dict)
        logger.info("Converted test dictionary using deep_str_conversion")
        
        try:
            json_str = json.dumps(converted_dict)
            logger.info(f"Successfully serialized converted dictionary to JSON: {len(json_str)} characters")
            
            # Test keys preservation
            original_keys = set(test_dict.keys())
            converted_keys = set(converted_dict.keys())
            logger.info(f"Original keys: {original_keys}")
            logger.info(f"Converted keys: {converted_keys}")
            
            if original_keys == converted_keys:
                logger.info("All keys preserved after conversion")
            else:
                logger.warning("Some keys were lost or changed during conversion")
                missing = original_keys - converted_keys
                if missing:
                    logger.warning(f"Missing keys: {missing}")
                new_keys = converted_keys - original_keys
                if new_keys:
                    logger.warning(f"New keys: {new_keys}")
            
            # Now test saving as AnnData
            success = save_as_anndata(converted_dict)
            return success
            
        except Exception as e:
            logger.error(f"Error serializing converted dictionary to JSON: {e}")
            return False
    
    return True

def test_real_metadata_files():
    """Test conversion with real metadata files from the project"""
    metadata_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json'
    success = True
    
    for dataset in ['adni', 'encode', 'gtex', 'mage']:
        json_path = os.path.join(metadata_dir, f"{dataset}_metadata.json")
        logger.info(f"Testing metadata file: {json_path}")
        
        if not os.path.exists(json_path):
            logger.warning(f"Metadata file not found: {json_path}")
            continue
        
        try:
            with open(json_path, 'r') as f:
                metadata = json.load(f)
            
            logger.info(f"Loaded metadata for {dataset}")
            
            # Convert metadata
            converted_metadata = deep_str_conversion(metadata)
            logger.info(f"Converted metadata for {dataset}")
            
            # Test saving as AnnData
            adata = sc.AnnData(X=np.zeros((2, 3)))
            adata.uns['metadata'] = converted_metadata
            
            test_file = f'/tmp/{dataset}_test_metadata.h5ad'
            adata.write_h5ad(test_file)
            logger.info(f"Successfully saved {dataset} metadata to {test_file}")
            
            # Load it back
            loaded_adata = sc.read_h5ad(test_file)
            logger.info(f"Successfully loaded {dataset} metadata from {test_file}")
            
        except Exception as e:
            logger.error(f"Error processing {dataset} metadata: {e}")
            success = False
    
    return success

def main():
    logger.info("=== Testing Metadata Conversion for Complex Nested Types ===")
    
    test_success = test_metadata_conversion()
    logger.info(f"Test dictionary conversion {'successful' if test_success else 'failed'}")
    
    real_success = test_real_metadata_files()
    logger.info(f"Real metadata files conversion {'successful' if real_success else 'failed'}")
    
    return 0 if test_success and real_success else 1

if __name__ == "__main__":
    sys.exit(main())
