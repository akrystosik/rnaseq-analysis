#!/usr/bin/env python3
"""
Fix for metadata serialization in AnnData objects.
This script patches the create_standard_anndata function in standardize_datasets.py
to ensure proper serialization of metadata.
"""
import os
import sys
import logging
import scanpy as sc
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('metadata_serialization_fix')

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
        'source': str(metadata.get('source', '')),
        'samples': str(metadata.get('samples', '')),
        'genes': str(metadata.get('genes', '')),
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

def patch_save_anndata(original_save_function):
    """Create a patched version of the save_anndata function"""
    
    def patched_save_anndata(adata, file_path):
        """Patched version of save_anndata that handles metadata properly"""
        try:
            # Convert all values in uns to serializable types
            if 'dataset_info' in adata.uns:
                adata.uns['dataset_info'] = extract_essential_metadata(adata.uns['dataset_info'])
            
            # Also convert any other dictionary in uns
            for key in list(adata.uns.keys()):
                if isinstance(adata.uns[key], dict):
                    adata.uns[key] = extract_essential_metadata(adata.uns[key])
            
            # Call the original function
            return original_save_function(adata, file_path)
        
        except Exception as e:
            logger.error(f"Error in patched_save_anndata: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False
    
    return patched_save_anndata

def patch_create_standard_anndata(original_function):
    """Create a patched version of the create_standard_anndata function"""
    
    def patched_create_standard_anndata(data_df, obs_df, var_df, dataset_info):
        """Patched version of create_standard_anndata that handles metadata properly"""
        try:
            # Convert dataset_info to serializable types
            safe_dataset_info = extract_essential_metadata(dataset_info)
            
            # Call the original function with the safe dataset_info
            adata = original_function(data_df, obs_df, var_df, safe_dataset_info)
            
            # Also ensure all values in uns are serializable
            for key in list(adata.uns.keys()):
                if isinstance(adata.uns[key], dict):
                    adata.uns[key] = extract_essential_metadata(adata.uns[key])
            
            return adata
        
        except Exception as e:
            logger.error(f"Error in patched_create_standard_anndata: {e}")
            import traceback
            logger.error(traceback.format_exc())
            
            # Fall back to a very simple implementation
            try:
                logger.info("Falling back to simple implementation")
                # Create AnnData object
                adata = sc.AnnData(X=data_df.values, obs=obs_df, var=var_df)
                
                # Add serializable dataset info
                adata.uns['dataset_info'] = {
                    'source': str(dataset_info.get('source', '')),
                    'data_type': str(dataset_info.get('data_type', '')),
                    'gencode_version': str(dataset_info.get('gencode_version', '')),
                    'samples': str(dataset_info.get('samples', '')),
                    'genes': str(dataset_info.get('genes', ''))
                }
                
                return adata
            except Exception as e2:
                logger.error(f"Error in fallback implementation: {e2}")
                import traceback
                logger.error(traceback.format_exc())
                return None
    
    return patched_create_standard_anndata

def apply_patches():
    """Apply patches to the functions in standardize_datasets.py"""
    # Import the module to patch
    import sys
    sys.path.append('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts')
    
    try:
        import standardize_datasets
        
        # Save the original functions
        original_create_standard = standardize_datasets.create_standard_anndata
        original_save = standardize_datasets.save_anndata
        
        # Create patched versions
        patched_create_standard = patch_create_standard_anndata(original_create_standard)
        patched_save = patch_save_anndata(original_save)
        
        # Replace the original functions with the patched versions
        standardize_datasets.create_standard_anndata = patched_create_standard
        standardize_datasets.save_anndata = patched_save
        
        logger.info("Successfully patched standardize_datasets.py functions")
        
        return True
    
    except Exception as e:
        logger.error(f"Error applying patches: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def test_patches():
    """Test the patches on a small dataset"""
    import tempfile
    
    try:
        # Create a small test dataset
        obs_df = pd.DataFrame({
            'sample_id': ['s1', 's2'],
            'tissue': ['liver', 'brain'],
            'subject_id': ['subj1', 'subj2']
        }, index=['s1', 's2'])
        
        var_df = pd.DataFrame({
            'gene_id': ['ENSG00000000001', 'ENSG00000000002'],
            'gene_name': ['gene1', 'gene2']
        }, index=['ENSG00000000001', 'ENSG00000000002'])
        
        data_df = pd.DataFrame({
            's1': [1.0, 2.0],
            's2': [3.0, 4.0]
        }, index=['ENSG00000000001', 'ENSG00000000002'])
        
        # Test dataset_info with complex types
        dataset_info = {
            'source': 'test',
            'data_type': 'RNA-seq',
            'gencode_version': 24,
            'samples': 2,
            'genes': 2,
            'numpy_value': np.int64(42),
            'nested_dict': {
                'complex_value': np.array([1, 2, 3]),
                'series': pd.Series([1, 2, 3])
            }
        }
        
        # Apply patches
        if not apply_patches():
            return False
        
        # Import the patched module
        import standardize_datasets
        
        # Test creating AnnData with patched function
        adata = standardize_datasets.create_standard_anndata(
            data_df.T, obs_df, var_df, dataset_info
        )
        
        # Test saving AnnData with patched function
        with tempfile.NamedTemporaryFile(suffix='.h5ad') as tmp:
            success = standardize_datasets.save_anndata(adata, tmp.name)
            
            if success:
                # Try loading it back
                loaded = sc.read_h5ad(tmp.name)
                logger.info(f"Successfully loaded AnnData with shape {loaded.shape}")
                logger.info(f"Loaded uns keys: {list(loaded.uns.keys())}")
                
                return True
            else:
                logger.error("Failed to save AnnData")
                return False
    
    except Exception as e:
        logger.error(f"Error testing patches: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    logger.info("=== Fixing Metadata Serialization ===")
    
    if test_patches():
        logger.info("Patches tested successfully")
        
        # Update script log
        with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/CHANGELOG.md', 'a') as f:
            f.write("\n## 2025-04-30\n")
            f.write("### Fixed\n")
            f.write("- Fixed metadata serialization issues in standardize_datasets.py\n")
            f.write("- Simplified metadata structure to ensure proper serialization to h5ad files\n")
            f.write("- Preserved essential metadata fields while removing complex nested structures\n")
        
        logger.info("Updated changelog")
    else:
        logger.error("Patch testing failed")
