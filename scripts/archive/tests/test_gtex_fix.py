#!/usr/bin/env python3
"""
GTEx Data Standardization Fix Verification Script

This script tests the fixed standardization for GTEx data with a small subset of data
to verify that metadata is properly retained in the AnnData object.
"""

import os
import sys
import pandas as pd
import numpy as np
import anndata as ad
import gzip
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,  # Set to DEBUG to see all details
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gtex_fix_verification')

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
GTEX_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz")
GTEX_SAMPLE_ATTRS = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
GTEX_SUBJECT_ATTRS = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
OUTPUT_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/gtex_test_fixed.h5ad")

# Import required functions from the standardization script
from standardize_datasets import (
    standardize_ensembl_id,
    load_gencode_mapping, 
    add_gencode_annotations,
    standardize_metadata,
    prepare_for_anndata,
    create_standard_anndata,
    save_anndata
)

def load_gtex_metadata_sample():
    """Load a sample of GTEx metadata for testing."""
    logger.info("Loading GTEx metadata sample")
    
    try:
        # Check if files exist
        if not GTEX_SAMPLE_ATTRS.exists():
            logger.error(f"GTEx sample attributes file not found: {GTEX_SAMPLE_ATTRS}")
            return pd.DataFrame()
            
        if not GTEX_SUBJECT_ATTRS.exists():
            logger.error(f"GTEx subject attributes file not found: {GTEX_SUBJECT_ATTRS}")
            return pd.DataFrame()
        
        # Load sample attributes - just take first 100 for testing
        logger.info(f"Loading sample attributes from {GTEX_SAMPLE_ATTRS}")
        sample_attrs = pd.read_csv(GTEX_SAMPLE_ATTRS, sep='\t', nrows=100)
        
        # Load subject attributes
        logger.info(f"Loading subject attributes from {GTEX_SUBJECT_ATTRS}")
        subject_attrs = pd.read_csv(GTEX_SUBJECT_ATTRS, sep='\t')
        
        # Extract SUBJID from SAMPID
        sample_attrs['SUBJID'] = sample_attrs['SAMPID'].str.extract(r'(GTEX-[A-Z0-9]+)')
        
        # Merge sample and subject attributes
        merged_attrs = pd.merge(sample_attrs, subject_attrs, on='SUBJID', how='left')
        
        # Standardize column names
        column_mapping = {
            'SAMPID': 'sample_id',
            'SUBJID': 'subject_id',
            'SMTSD': 'tissue',  # Detailed tissue site
            'SMTS': 'broad_tissue',  # Broad tissue type
            'SEX': 'sex',
            'AGE': 'age',
            'SMGEBTCH': 'batch',
            'SMRIN': 'rna_integrity_number',
            'SMTSISCH': 'ischemic_time',
            'SMNABTCH': 'array_batch'
        }
        
        # Rename columns that exist in the DataFrame
        existing_cols = [col for col in column_mapping.keys() if col in merged_attrs.columns]
        merged_attrs = merged_attrs.rename(columns={col: column_mapping[col] for col in existing_cols})
        
        # Ensure we have a tissue column
        if 'tissue' not in merged_attrs.columns and 'broad_tissue' in merged_attrs.columns:
            merged_attrs['tissue'] = merged_attrs['broad_tissue']
            logger.info("Using broad_tissue as tissue column")
        
        # Set sample ID as index
        merged_attrs.set_index('sample_id', inplace=True)
        
        logger.info(f"Loaded metadata for {len(merged_attrs)} GTEx samples")
        return merged_attrs
    
    except Exception as e:
        logger.error(f"Error loading GTEx metadata: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return pd.DataFrame()

def load_gtex_expression_sample():
    """
    Load a small sample of GTEx expression data for testing.
    
    Returns:
    --------
    pandas.DataFrame
        Sample expression data with genes as rows and samples as columns
    """
    logger.info(f"Loading GTEx expression data sample from {GTEX_FILE}")
    
    try:
        # Read a limited number of genes and samples
        max_genes = 1000  # Number of genes to read
        
        with gzip.open(GTEX_FILE, 'rt') as f:
            # Read header lines
            version_line = next(f)
            dims_line = next(f)
            header_line = next(f).strip().split('\t')
            
            # Get sample IDs from header (first 50 samples)
            sample_ids = header_line[2:52]  # First 50 samples, skip Name and Description
            logger.info(f"Using {len(sample_ids)} samples for testing")
            
            # Read specified number of genes
            all_genes = []
            all_data = []
            
            for i in range(max_genes):
                try:
                    line = next(f)
                    fields = line.strip().split('\t')
                    gene_id = fields[0]  # First column is the gene ID
                    # Skip Description column (position 1)
                    expression = [float(x) for x in fields[2:52]]  # Only get first 50 samples
                    
                    all_genes.append(gene_id)
                    all_data.append(expression)
                except StopIteration:
                    break
            
            # Create the DataFrame
            data_df = pd.DataFrame(all_data, index=all_genes, columns=sample_ids)
            logger.info(f"Loaded sample expression data with {data_df.shape[0]} genes and {data_df.shape[1]} samples")
            
            return data_df
    
    except Exception as e:
        logger.error(f"Error loading GTEx expression data sample: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return pd.DataFrame()

def main():
    # Load a sample of GTEx data
    logger.info("Starting GTEx fix verification test")
    
    # Load expression data sample
    expr_df = load_gtex_expression_sample()
    if expr_df.empty:
        logger.error("Failed to load GTEx expression data sample")
        return
    
    # Load metadata sample
    metadata_df = load_gtex_metadata_sample()
    if metadata_df.empty:
        logger.warning("No GTEx metadata found, creating minimal metadata")
        metadata_df = pd.DataFrame(index=expr_df.columns)
        metadata_df['sample_id'] = metadata_df.index
    
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
    
    # Create a small standardized expression matrix
    sample_count = len(common_samples)
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    gene_count = len(unique_std_ids)
    
    logger.info(f"Creating standardized expression matrix with {gene_count} genes and {sample_count} samples")
    
    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    
    # Initialize expression matrix
    expr_matrix = np.zeros((gene_count, sample_count), dtype=np.float32)
    
    # Fill the matrix
    for sample_idx, sample_id in enumerate(common_samples):
        # Process genes for this sample
        for gene_id, expr_val in zip(expr_df.index, expr_df[sample_id]):
            std_id = gene_id_mapping[gene_id]
            if std_id in std_id_to_idx:  # Handle case where std_id might not be in our subset
                std_idx = std_id_to_idx[std_id]
                expr_matrix[std_idx, sample_idx] = max(expr_matrix[std_idx, sample_idx], expr_val)
    
    # Create DataFrame with the results
    std_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=common_samples)
    
    # Create variable (gene) metadata DataFrame
    var_df = pd.DataFrame(index=unique_std_ids)
    var_df['gene_id'] = var_df.index
    
    # Create a reverse mapping from std_id to original ids
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)
    
    # Add original IDs as a concatenated string
    var_df['original_ids'] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])) if x in reverse_mapping else "")
    
    # Add GENCODE annotations
    logger.info("Adding GENCODE annotations")
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)
    
    # Standardize observation metadata
    logger.info("Standardizing metadata")
    obs_df = standardize_metadata(metadata_df, 'GTEx')
    
    # Add missing standard fields
    if 'data_type' not in obs_df.columns:
        obs_df['data_type'] = 'RNA-seq'
    if 'expression_unit' not in obs_df.columns:
        obs_df['expression_unit'] = 'TPM'
    
    # Dataset info for uns
    dataset_info = {
        'source': 'GTEx',
        'version': 'v10',
        'gencode_version': 24,
        'data_type': 'RNA-seq',
        'expression_unit': 'TPM',
        'samples': len(obs_df),
        'genes': len(var_df),
        'tissue_count': len(obs_df['tissue'].unique()) if 'tissue' in obs_df.columns else 0,
        'subject_count': len(obs_df['subject_id'].unique()) if 'subject_id' in obs_df.columns else 0
    }
    
    # Create AnnData object
    logger.info("Creating AnnData object")
    
    # Debug info
    logger.info(f"Before AnnData creation:")
    logger.info(f"  - var_df shape: {var_df.shape}")
    logger.info(f"  - var_df columns: {list(var_df.columns)}")
    logger.info(f"  - obs_df shape: {obs_df.shape}")
    logger.info(f"  - obs_df columns: {list(obs_df.columns)}")
    
    # Creating AnnData with improved functions
    adata = create_standard_anndata(std_expr_df.T, obs_df, var_df, dataset_info)
    
    # Check if creation was successful
    logger.info(f"After AnnData creation:")
    logger.info(f"  - adata.var shape: {adata.var.shape}")
    logger.info(f"  - adata.var columns: {list(adata.var.columns)}")
    logger.info(f"  - adata.obs shape: {adata.obs.shape}")
    logger.info(f"  - adata.obs columns: {list(adata.obs.columns)}")
    
    # Save the standardized AnnData
    if save_anndata(adata, OUTPUT_FILE):
        logger.info(f"Successfully saved test GTEx data to {OUTPUT_FILE}")
        
        # Reload and verify
        logger.info("Reloading saved file to verify integrity")
        try:
            reload_adata = ad.read_h5ad(OUTPUT_FILE)
            logger.info(f"Verification - reloaded data:")
            logger.info(f"  - var shape: {reload_adata.var.shape}")
            logger.info(f"  - var columns: {list(reload_adata.var.columns)}")
            logger.info(f"  - obs shape: {reload_adata.obs.shape}")
            logger.info(f"  - obs columns: {list(reload_adata.obs.columns)}")
            
            # Success if metadata is preserved
            if reload_adata.var.shape[1] > 0 and reload_adata.obs.shape[1] > 0:
                logger.info("✅ TEST PASSED: Metadata was successfully preserved!")
            else:
                logger.error("❌ TEST FAILED: Metadata was lost after reload")
        except Exception as e:
            logger.error(f"Error verifying saved file: {e}")
    
    logger.info("GTEx fix verification test completed")

if __name__ == "__main__":
    main()