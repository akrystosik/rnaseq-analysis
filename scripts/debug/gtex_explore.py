#!/usr/bin/env python3
"""
GTEx AnnData Diagnostic Script

This script examines the problematic GTEx standardized AnnData file 
and compares it with the original GTEx source data to identify 
structural issues.
"""

import os
import pandas as pd
import numpy as np
import anndata as ad
import gzip
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gtex_diagnostic')

# Define paths
GTEX_STD_PATH = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/gtex_standardized.h5ad"
GTEX_RAW_PATH = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"
GENCODE_MAPPING_FILE = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gencode_v24_complete_mapping.csv"
GTEx_SAMPLE_ATTRS = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"

def load_standardized_gtex():
    """Load the standardized GTEx AnnData file and examine its structure."""
    logger.info(f"Loading standardized GTEx data from {GTEX_STD_PATH}")
    try:
        adata = ad.read_h5ad(GTEX_STD_PATH)
        
        # Basic structure check
        logger.info(f"AnnData dimensions: {adata.shape}")
        logger.info(f"var DataFrame shape: {adata.var.shape}")
        logger.info(f"obs DataFrame shape: {adata.obs.shape}")
        
        # Check var DataFrame contents
        logger.info(f"var DataFrame columns: {list(adata.var.columns)}")
        logger.info(f"var DataFrame index type: {type(adata.var.index)}")
        logger.info(f"First 5 var indices: {list(adata.var.index[:5])}")
        
        # Check obs DataFrame contents
        logger.info(f"obs DataFrame columns: {list(adata.obs.columns)}")
        logger.info(f"obs DataFrame index type: {type(adata.obs.index)}")
        logger.info(f"First 5 obs indices: {list(adata.obs.index[:5])}")
        
        # Check X matrix values
        logger.info(f"X matrix type: {type(adata.X)}")
        logger.info(f"X matrix shape: {adata.X.shape}")
        logger.info(f"X matrix data range: min={adata.X.min():.2f}, max={adata.X.max():.2f}, mean={adata.X.mean():.2f}")
        
        # Check uns metadata
        logger.info(f"uns keys: {list(adata.uns.keys())}")
        if 'dataset_info' in adata.uns:
            logger.info(f"dataset_info: {adata.uns['dataset_info']}")
        
        return adata
    except Exception as e:
        logger.error(f"Error loading standardized GTEx data: {e}")
        return None

def sample_gtex_raw_data():
    """Sample the raw GTEx data to compare with standardized version."""
    logger.info(f"Sampling raw GTEx data from {GTEX_RAW_PATH}")
    try:
        # Read the first few lines to get header and sample rows
        with gzip.open(GTEX_RAW_PATH, 'rt') as f:
            # Skip version line
            version = next(f).strip()
            # Read dimensions line
            dims_line = next(f).strip()
            dims = dims_line.split('\t')
            n_genes, n_samples = int(dims[0]), int(dims[1])
            logger.info(f"GTEx raw file dimensions: {n_genes} genes, {n_samples} samples")
            
            # Read column headers
            header_line = next(f).strip()
            headers = header_line.split('\t')
            logger.info(f"First few column headers: {headers[:5]}")
            
            # Sample a few gene rows
            sample_rows = []
            for i in range(5):
                sample_rows.append(next(f).strip().split('\t'))
            
            logger.info(f"First gene ID: {sample_rows[0][0]}")
            logger.info(f"Sample gene IDs: {[row[0] for row in sample_rows]}")
        
        # Also load GENCODE mapping to check ID formats
        if os.path.exists(GENCODE_MAPPING_FILE):
            gencode_df = pd.read_csv(GENCODE_MAPPING_FILE)
            logger.info(f"GENCODE mapping loaded with {len(gencode_df)} entries")
            logger.info(f"GENCODE columns: {list(gencode_df.columns)}")
            logger.info(f"First few GENCODE gene IDs: {list(gencode_df['gene_id'].head())}")
        
        # Check GTEx sample attributes
        if os.path.exists(GTEx_SAMPLE_ATTRS):
            sample_attrs = pd.read_csv(GTEx_SAMPLE_ATTRS, sep='\t')
            logger.info(f"GTEx sample attributes loaded with {len(sample_attrs)} entries")
            logger.info(f"Sample attribute columns: {list(sample_attrs.columns)}")
            logger.info(f"First few sample IDs: {list(sample_attrs['SAMPID'].head())}")
        
        return {
            "n_genes": n_genes,
            "n_samples": n_samples,
            "sample_gene_ids": [row[0] for row in sample_rows]
        }
    except Exception as e:
        logger.error(f"Error sampling raw GTEx data: {e}")
        return None

def diagnose_problems(adata, raw_data):
    """Diagnose problems in the standardized GTEx data."""
    if adata is None or raw_data is None:
        logger.error("Cannot diagnose problems: Missing data")
        return
    
    logger.info("=== DIAGNOSIS RESULTS ===")
    
    # Check if dimensions match
    logger.info(f"Standardized data has {adata.n_vars} genes and {adata.n_obs} samples")
    logger.info(f"Raw data has {raw_data['n_genes']} genes and {raw_data['n_samples']} samples")
    
    # Check for metadata presence
    if adata.var.shape[1] == 0:
        logger.error("PROBLEM: var DataFrame is empty - no gene metadata")
    
    if adata.obs.shape[1] == 0:
        logger.error("PROBLEM: obs DataFrame is empty - no sample metadata")
    
    # Check for proper gene identifiers
    has_gene_ids = False
    for sample_gene_id in raw_data['sample_gene_ids']:
        # Check if raw gene ID is preserved somewhere in the standardized data
        if sample_gene_id in adata.var.index:
            has_gene_ids = True
            logger.info(f"Found raw gene ID {sample_gene_id} in standardized data index")
            break
        elif 'gene_id' in adata.var.columns and sample_gene_id in adata.var['gene_id'].values:
            has_gene_ids = True
            logger.info(f"Found raw gene ID {sample_gene_id} in standardized data gene_id column")
            break
        elif 'original_ids' in adata.var.columns:
            for orig_ids in adata.var['original_ids']:
                if sample_gene_id in str(orig_ids).split(';'):
                    has_gene_ids = True
                    logger.info(f"Found raw gene ID {sample_gene_id} in standardized data original_ids column")
                    break
    
    if not has_gene_ids:
        logger.error("PROBLEM: Cannot find original gene IDs in standardized data")
    
    # Determine possible causes
    logger.info("=== POSSIBLE CAUSES ===")
    
    if adata.var.shape[1] == 0 and adata.obs.shape[1] == 0:
        logger.error(
            "1. AnnData creation issue: Both var and obs DataFrames are empty, suggesting "
            "the create_standard_anndata() function isn't properly transferring metadata."
        )
    
    logger.info(
        "2. Index conversion issue: The pipeline may be losing the mapping between "
        "standardized gene IDs and original gene IDs during processing."
    )
    
    logger.info(
        "3. Saving/loading issue: Metadata columns may be lost during the saving or "
        "loading of the AnnData object."
    )
    
    logger.info("=== RECOMMENDED FIXES ===")
    
    logger.info(
        "1. Inspect the create_standard_anndata() function to ensure metadata DataFrames "
        "are properly passed to AnnData constructor."
    )
    
    logger.info(
        "2. Add explicit checks for empty DataFrames before saving AnnData objects."
    )
    
    logger.info(
        "3. Ensure proper mapping between standardized IDs and original IDs is maintained "
        "throughout the processing pipeline."
    )

if __name__ == "__main__":
    logger.info("Running GTEx AnnData diagnostic")
    
    # Load standardized data
    adata = load_standardized_gtex()
    
    # Sample raw data
    raw_data = sample_gtex_raw_data()
    
    # Diagnose problems
    diagnose_problems(adata, raw_data)