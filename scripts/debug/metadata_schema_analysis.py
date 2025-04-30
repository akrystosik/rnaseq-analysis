#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/metadata_schema_analysis.py
#!/usr/bin/env python3
"""
Metadata Schema Analysis for Combined Dataset
"""
import sys
import pandas as pd
import numpy as np
import anndata as ad
import logging
from pathlib import Path
from collections import Counter

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('metadata_analysis')

# File paths
COMBINED_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_standardized.h5ad")
GTEX_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/gtex_standardized.h5ad")

def analyze_dataset(file_path, name="dataset"):
    """Analyze the metadata schema in an AnnData object."""
    logger.info(f"Loading {name} from {file_path}")
    
    if not file_path.exists():
        logger.error(f"File not found: {file_path}")
        return
    
    adata = ad.read_h5ad(file_path)
    logger.info(f"{name} dimensions: {adata.shape[0]} samples x {adata.shape[1]} genes")
    
    # Analyze observation (sample) metadata
    logger.info(f"\n=== {name} Sample Metadata Analysis ===")
    logger.info(f"Total columns in obs: {len(adata.obs.columns)}")
    logger.info(f"Columns: {', '.join(adata.obs.columns)}")
    
    # If source_dataset exists, analyze by dataset
    if 'source_dataset' in adata.obs.columns:
        dataset_counts = adata.obs['source_dataset'].value_counts().to_dict()
        logger.info("Sample counts by dataset:")
        for dataset, count in dataset_counts.items():
            logger.info(f"  {dataset}: {count} samples")
    
    # Check core metadata fields
    core_fields = [
        'subject_id', 'sex', 'age', 'tissue', 'dataset', 'data_type',
        'expression_unit', 'species', 'tissue_ontology', 'assay_ontology', 
        'developmental_stage_ontology', 'species_ontology'
    ]
    
    logger.info("\nCore metadata fields completeness:")
    for field in core_fields:
        if field in adata.obs.columns:
            non_empty = (~adata.obs[field].isna() & (adata.obs[field] != "")).sum()
            percentage = (non_empty / len(adata.obs)) * 100
            logger.info(f"  {field}: {non_empty}/{len(adata.obs)} ({percentage:.1f}%)")
            
            # For categorical fields, show category counts
            if pd.api.types.is_categorical_dtype(adata.obs[field]):
                top_values = adata.obs[field].value_counts().head(5).to_dict()
                logger.info(f"    Top values: {top_values}")
        else:
            logger.info(f"  {field}: Not present")
    
    # Analyze variable (gene) metadata
    logger.info(f"\n=== {name} Gene Metadata Analysis ===")
    logger.info(f"Total columns in var: {len(adata.var.columns)}")
    logger.info(f"Columns: {', '.join(adata.var.columns)}")
    
    # Check gene ID format
    if len(adata.var_names) > 0:
        ensembl_pattern = "^ENSG\\d+"
        ensembl_count = sum(adata.var_names.str.match(ensembl_pattern))
        logger.info(f"Ensembl ID format: {ensembl_count}/{len(adata.var_names)} ({ensembl_count/len(adata.var_names)*100:.1f}%)")
    
    # Check GENCODE mapping stats
    if 'mapping_source' in adata.var.columns:
        mapping_counts = adata.var['mapping_source'].value_counts().to_dict()
        logger.info("GENCODE mapping statistics:")
        for source, count in mapping_counts.items():
            logger.info(f"  {source}: {count} genes ({count/len(adata.var)*100:.1f}%)")
    
    # Check unstructured metadata (uns)
    logger.info(f"\n=== {name} Unstructured Metadata (uns) ===")
    if hasattr(adata, 'uns') and adata.uns:
        logger.info(f"uns keys: {list(adata.uns.keys())}")
        
        # If dataset_info exists, show it
        if 'dataset_info' in adata.uns:
            logger.info("dataset_info contents:")
            for key, value in adata.uns['dataset_info'].items():
                logger.info(f"  {key}: {value}")
        
        # If ontology_mappings exists, show it
        if 'ontology_mappings' in adata.uns:
            logger.info("Ontology mappings present:")
            for key in adata.uns['ontology_mappings'].keys():
                logger.info(f"  {key}")

    return adata

def main():
    """Analyze both the combined dataset and GTEx dataset."""
    # Analyze the combined dataset
    combined_adata = analyze_dataset(COMBINED_FILE, "Combined Dataset")
    
    # Analyze the GTEx dataset
    gtex_adata = analyze_dataset(GTEX_FILE, "GTEx Dataset")
    
    # Compare common genes if both datasets were loaded
    if combined_adata is not None and gtex_adata is not None:
        logger.info("\n=== Dataset Comparison ===")
        combined_genes = set(combined_adata.var_names)
        gtex_genes = set(gtex_adata.var_names)
        common_genes = combined_genes.intersection(gtex_genes)
        
        logger.info(f"Combined dataset: {len(combined_genes)} genes")
        logger.info(f"GTEx dataset: {len(gtex_genes)} genes")
        logger.info(f"Common genes: {len(common_genes)}")
        logger.info(f"Gene overlap: {len(common_genes)/len(combined_genes)*100:.1f}% of combined, {len(common_genes)/len(gtex_genes)*100:.1f}% of GTEx")

if __name__ == "__main__":
    main()