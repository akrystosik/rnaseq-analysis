#!/usr/bin/env python3
"""
GTEx Single-Cell RNA-seq Processing with Pseudobulking

This script processes GTEx single-cell RNA-seq data to create pseudobulk samples
by aggregating cells by subject_id and cell_type. It applies the same
standardization pipeline as bulk RNA-seq data.

Usage:
  python process_gtex_single_cell.py --input /path/to/gtex_sc.h5ad \
                                    --output /path/to/gtex_sc_pseudobulk.h5ad \
                                    --metadata-dir /path/to/metadata/json
"""

import os
import sys
import logging
import argparse
import json
from pathlib import Path
import pandas as pd
import scanpy as sc
import numpy as np
from datetime import datetime

# Import shared utilities
from rnaseq_utils import (
    standardize_metadata,
    load_dataset_specific_metadata,
    apply_dataset_specific_metadata,
    add_gencode_annotations,
    load_gencode_mapping,
    create_standard_anndata,
    save_anndata
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gtex_single_cell_processor')

def load_celltype_to_cl_mapping(metadata_dir):
    """Load cell type to Cell Ontology mapping."""
    mapping_file = Path(metadata_dir) / "celltype_to_cl.json"
    if not mapping_file.exists():
        logger.warning(f"Cell type mapping file not found: {mapping_file}")
        return {}
    
    try:
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
            return mapping_data.get('mappings', {})
    except Exception as e:
        logger.error(f"Error loading cell type mappings: {e}")
        return {}

def standardize_cell_type_names(cell_type_series):
    """Standardize cell type names for consistency."""
    # Basic standardization rules
    standardized = cell_type_series.astype(str).str.strip()
    
    # Example mappings (expand as needed based on actual GTEx cell types)
    replacements = {
        'unknown': 'unknown',
        'nan': 'unknown',
        'None': 'unknown',
        '': 'unknown'
    }
    
    for old, new in replacements.items():
        standardized = standardized.replace(old, new, regex=False)
    
    return standardized

def create_pseudobulk_samples(adata, subject_id_col='subject_id', cell_type_col='cell_type'):
    """
    Create pseudobulk samples by aggregating cells by subject_id and cell_type.
    
    Args:
        adata: AnnData object with single-cell data
        subject_id_col: Column name for subject ID
        cell_type_col: Column name for cell type
        
    Returns:
        Dictionary with pseudobulk data and metadata
    """
    logger.info(f"Creating pseudobulk samples from {adata.n_obs} cells and {adata.n_vars} genes")
    
    # Ensure we have the required columns
    if subject_id_col not in adata.obs.columns:
        logger.error(f"Required column '{subject_id_col}' not found in adata.obs")
        return None
        
    if cell_type_col not in adata.obs.columns:
        logger.error(f"Required column '{cell_type_col}' not found in adata.obs")
        return None
    
    # Standardize cell type names
    adata.obs[cell_type_col] = standardize_cell_type_names(adata.obs[cell_type_col])
    
    # Create sample grouping key
    adata.obs['pseudobulk_sample_id'] = (adata.obs[subject_id_col].astype(str) + '_' + 
                                        adata.obs[cell_type_col].astype(str))
    
    # Get unique pseudobulk samples
    unique_samples = adata.obs['pseudobulk_sample_id'].unique()
    logger.info(f"Creating {len(unique_samples)} pseudobulk samples")
    
    # Initialize pseudobulk expression matrix
    pseudobulk_expr = np.zeros((len(unique_samples), adata.n_vars))
    pseudobulk_metadata = []
    
    for i, sample_id in enumerate(unique_samples):
        # Get cells for this pseudobulk sample
        sample_mask = adata.obs['pseudobulk_sample_id'] == sample_id
        sample_cells = adata[sample_mask]
        
        if sample_cells.n_obs == 0:
            continue
            
        # Sum expression across cells (pseudobulking)
        pseudobulk_expr[i, :] = np.array(sample_cells.X.sum(axis=0)).flatten()
        
        # Create metadata for this pseudobulk sample
        # Take the first cell's metadata as representative
        first_cell = sample_cells.obs.iloc[0]
        
        metadata = {
            'sample_id': sample_id,
            'subject_id': first_cell[subject_id_col],
            'cell_type': first_cell[cell_type_col],
            'n_cells_aggregated': sample_cells.n_obs,
            'dataset': 'GTEx_scRNA_pseudobulk',
            'data_type': 'RNA-seq',
            'expression_unit': 'pseudobulk_counts',
            'assay_type': 'single-cell RNA-seq (pseudobulked)',
        }
        
        # Copy other metadata fields if available
        for col in ['tissue', 'sex', 'age', 'species', 'race', 'ethnicity']:
            if col in sample_cells.obs.columns:
                metadata[col] = first_cell.get(col, 'unknown')
        
        pseudobulk_metadata.append(metadata)
    
    # Create pseudobulk AnnData object
    pseudobulk_obs = pd.DataFrame(pseudobulk_metadata).set_index('sample_id')
    pseudobulk_var = adata.var.copy()
    
    pseudobulk_adata = sc.AnnData(
        X=pseudobulk_expr,
        obs=pseudobulk_obs,
        var=pseudobulk_var
    )
    
    logger.info(f"Created pseudobulk dataset: {pseudobulk_adata.n_obs} samples x {pseudobulk_adata.n_vars} genes")
    
    return pseudobulk_adata

def process_gtex_single_cell(input_file, output_file, metadata_dir=None):
    """
    Process GTEx single-cell data to create standardized pseudobulk samples.
    
    Args:
        input_file: Path to GTEx single-cell H5AD file
        output_file: Path to save pseudobulk H5AD file
        metadata_dir: Directory containing metadata JSON files
    """
    logger.info(f"Processing GTEx single-cell data from {input_file}")
    
    # Load single-cell data
    try:
        adata_sc = sc.read_h5ad(input_file)
        logger.info(f"Loaded single-cell data: {adata_sc.n_obs} cells x {adata_sc.n_vars} genes")
    except Exception as e:
        logger.error(f"Error loading single-cell data from {input_file}: {e}")
        return None
    
    # Create pseudobulk samples
    adata_pb = create_pseudobulk_samples(adata_sc)
    if adata_pb is None:
        logger.error("Failed to create pseudobulk samples")
        return None
    
    # Load GENCODE mapping for gene annotation
    gencode_map = load_gencode_mapping()
    if gencode_map is not None:
        adata_pb.var = add_gencode_annotations(adata_pb.var, gencode_map)
    
    # Standardize metadata
    if metadata_dir:
        # Load GTEx-specific metadata configuration
        gtex_sc_config = load_dataset_specific_metadata(metadata_dir, "gtex_scrnaseq") or {}
        
        # Apply dataset-specific metadata
        if gtex_sc_config:
            adata_pb = apply_dataset_specific_metadata(adata_pb, gtex_sc_config)
        
        # Apply general metadata standardization
        from rnaseq_utils import load_mappings
        mappings = load_mappings()
        adata_pb.obs = standardize_metadata(adata_pb.obs, "GTEx_scRNA_pseudobulk", mappings)
        
        # Apply cell type ontology mapping
        celltype_mappings = load_celltype_to_cl_mapping(metadata_dir)
        if celltype_mappings and 'cell_type' in adata_pb.obs.columns:
            logger.info("Applying cell type ontology mapping")
            
            # Initialize ontology columns
            adata_pb.obs['cell_type_ontology_term_id'] = ''
            adata_pb.obs['cell_type_ontology_confidence'] = 'none'
            
            for idx, row in adata_pb.obs.iterrows():
                cell_type = str(row['cell_type']).strip().lower()
                if cell_type in celltype_mappings:
                    ontology_info = celltype_mappings[cell_type]
                    if isinstance(ontology_info, dict):
                        adata_pb.obs.loc[idx, 'cell_type_ontology_term_id'] = ontology_info.get('cl_term_id', '')
                        adata_pb.obs.loc[idx, 'cell_type_ontology_confidence'] = ontology_info.get('confidence', 'medium')
            
            mapped_count = (adata_pb.obs['cell_type_ontology_term_id'] != '').sum()
            logger.info(f"Mapped {mapped_count}/{adata_pb.n_obs} samples to cell type ontology terms")
    
    # Add harmonized metadata
    adata_pb.uns['harmonized_gencode_version'] = 'v24'
    adata_pb.uns['harmonized_reference_genome'] = 'hg38'
    adata_pb.uns['original_data_type'] = 'single-cell RNA-seq'
    adata_pb.uns['processing_method'] = 'pseudobulking by subject_id and cell_type'
    adata_pb.uns['processing_date'] = datetime.now().isoformat()
    
    # Save pseudobulk data
    if save_anndata(adata_pb, output_file):
        logger.info(f"Successfully saved GTEx pseudobulk data: {adata_pb.n_obs} samples x {adata_pb.n_vars} genes")
        logger.info(f"Output file: {output_file}")
        return adata_pb
    else:
        logger.error("Failed to save pseudobulk data")
        return None

def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(description='Process GTEx single-cell RNA-seq data')
    parser.add_argument('--input', required=True, help='Input GTEx single-cell H5AD file')
    parser.add_argument('--output', required=True, help='Output pseudobulk H5AD file')
    parser.add_argument('--metadata-dir', default='/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json',
                       help='Directory containing metadata JSON files')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Process data
    result = process_gtex_single_cell(args.input, args.output, args.metadata_dir)
    
    if result is None:
        logger.error("Processing failed")
        sys.exit(1)
    else:
        logger.info("Processing completed successfully")

if __name__ == '__main__':
    main()