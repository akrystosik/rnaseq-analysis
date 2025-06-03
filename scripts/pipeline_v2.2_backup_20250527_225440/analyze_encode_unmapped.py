#!/usr/bin/env python3
"""
Analyze ENCODE Unmapped Gene IDs

This script investigates the 22,500 unmapped ENCODE gene IDs to understand:
1. What types of IDs they are (numeric Entrez, versioned Ensembl, etc.)
2. Whether they have low/zero expression values
3. Whether they could be mapped with additional processing

Usage:
    python analyze_encode_unmapped.py \
        --encode-h5ad /path/to/encode_standardized_preprocessed.h5ad \
        --encode-mapping /path/to/encode_id_to_ensembl_mapping.csv \
        --reference-mapping /path/to/gencode_v24_complete_mapping.csv
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('analyze_encode_unmapped')

def analyze_expression_levels(adata, unmapped_mask):
    """Analyze expression levels of unmapped genes."""
    logger.info("=== Analyzing Expression Levels ===")
    
    unmapped_genes = adata.var_names[unmapped_mask]
    mapped_genes = adata.var_names[~unmapped_mask]
    
    # Calculate statistics for unmapped genes
    unmapped_expr = adata[:, unmapped_mask].X
    mapped_expr = adata[:, ~unmapped_mask].X
    
    # Convert to dense if sparse
    if hasattr(unmapped_expr, 'toarray'):
        unmapped_expr = unmapped_expr.toarray()
    if hasattr(mapped_expr, 'toarray'):
        mapped_expr = mapped_expr.toarray()
    
    # Calculate mean expression per gene
    unmapped_means = np.mean(unmapped_expr, axis=0)
    mapped_means = np.mean(mapped_expr, axis=0)
    
    logger.info(f"Unmapped genes: {len(unmapped_genes)}")
    logger.info(f"Mapped genes: {len(mapped_genes)}")
    
    # Expression statistics
    logger.info(f"Unmapped gene expression - Mean: {np.mean(unmapped_means):.4f}, Median: {np.median(unmapped_means):.4f}")
    logger.info(f"Mapped gene expression - Mean: {np.mean(mapped_means):.4f}, Median: {np.median(mapped_means):.4f}")
    
    # Count zero/low expression genes
    zero_expr_unmapped = np.sum(unmapped_means == 0)
    low_expr_unmapped = np.sum(unmapped_means < 0.1)
    zero_expr_mapped = np.sum(mapped_means == 0)
    low_expr_mapped = np.sum(mapped_means < 0.1)
    
    logger.info(f"Zero expression - Unmapped: {zero_expr_unmapped}/{len(unmapped_genes)} ({100*zero_expr_unmapped/len(unmapped_genes):.1f}%)")
    logger.info(f"Zero expression - Mapped: {zero_expr_mapped}/{len(mapped_genes)} ({100*zero_expr_mapped/len(mapped_genes):.1f}%)")
    logger.info(f"Low expression (<0.1) - Unmapped: {low_expr_unmapped}/{len(unmapped_genes)} ({100*low_expr_unmapped/len(unmapped_genes):.1f}%)")
    logger.info(f"Low expression (<0.1) - Mapped: {low_expr_mapped}/{len(mapped_genes)} ({100*low_expr_mapped/len(mapped_genes):.1f}%)")
    
    return {
        'unmapped_zero_expr': zero_expr_unmapped,
        'unmapped_low_expr': low_expr_unmapped,
        'unmapped_total': len(unmapped_genes),
        'mapped_zero_expr': zero_expr_mapped,
        'mapped_low_expr': low_expr_mapped,
        'mapped_total': len(mapped_genes)
    }

def analyze_id_formats(adata):
    """Analyze the format and types of gene IDs."""
    logger.info("=== Analyzing Gene ID Formats ===")
    
    # Get mapping sources
    mapping_sources = adata.var['mapping_source'].value_counts()
    logger.info("Mapping source breakdown:")
    for source, count in mapping_sources.items():
        percentage = 100 * count / len(adata.var)
        logger.info(f"  {source}: {count} ({percentage:.1f}%)")
    
    # Analyze unmapped IDs specifically
    unmapped_mask = adata.var['mapping_source'] == 'unmapped'
    unmapped_original_ids = adata.var.loc[unmapped_mask, 'original_id_from_source']
    unmapped_input_ids = adata.var.loc[unmapped_mask, 'original_id_from_input']
    
    logger.info(f"\nAnalyzing {len(unmapped_original_ids)} unmapped gene IDs:")
    
    # Categorize original IDs
    entrez_like = 0
    ensembl_like = 0
    complex_like = 0
    other = 0
    
    for orig_id in unmapped_original_ids:
        orig_str = str(orig_id)
        if orig_str.isdigit():
            entrez_like += 1
        elif 'ENSG' in orig_str:
            ensembl_like += 1
        elif '|' in orig_str or 'ENST' in orig_str:
            complex_like += 1
        else:
            other += 1
    
    logger.info(f"  Numeric (Entrez-like): {entrez_like} ({100*entrez_like/len(unmapped_original_ids):.1f}%)")
    logger.info(f"  Contains ENSG: {ensembl_like} ({100*ensembl_like/len(unmapped_original_ids):.1f}%)")
    logger.info(f"  Complex (contains |): {complex_like} ({100*complex_like/len(unmapped_original_ids):.1f}%)")
    logger.info(f"  Other format: {other} ({100*other/len(unmapped_original_ids):.1f}%)")
    
    # Show examples
    logger.info(f"\nExamples of unmapped IDs:")
    logger.info(f"  First 10 original_id_from_source: {unmapped_original_ids.head(10).tolist()}")
    logger.info(f"  First 10 original_id_from_input: {unmapped_input_ids.head(10).tolist()}")
    
    return unmapped_mask, {
        'entrez_like': entrez_like,
        'ensembl_like': ensembl_like,
        'complex_like': complex_like,
        'other': other
    }

def check_mapping_potential(unmapped_original_ids, encode_mapping_df, reference_mapping_df):
    """Check if unmapped IDs could potentially be mapped with different strategies."""
    logger.info("=== Checking Mapping Potential ===")
    
    # Load the mapping files to understand their structure
    logger.info(f"ENCODE mapping file has {len(encode_mapping_df)} entries")
    logger.info(f"Reference mapping file has {len(reference_mapping_df)} entries")
    
    # Sample unmapped IDs for analysis
    sample_unmapped = unmapped_original_ids.head(20)
    logger.info(f"\nAnalyzing sample of unmapped IDs:")
    
    could_map_count = 0
    
    for orig_id in sample_unmapped:
        orig_str = str(orig_id)
        logger.info(f"  ID: {orig_str}")
        
        # Check if it's in ENCODE mapping
        encode_matches = encode_mapping_df[encode_mapping_df['original_id'] == orig_str]
        if not encode_matches.empty:
            logger.info(f"    -> Found in ENCODE mapping: {encode_matches['ensembl_id'].iloc[0]}")
            could_map_count += 1
            continue
        
        # Check if it's numeric and in reference mapping
        if orig_str.isdigit():
            ref_matches = reference_mapping_df[reference_mapping_df['numeric_id'] == orig_str]
            if not ref_matches.empty:
                logger.info(f"    -> Found in reference mapping: {ref_matches['gene_id'].iloc[0]}")
                could_map_count += 1
                continue
        
        # Check if it contains ENSG and could be processed
        if 'ENSG' in orig_str:
            # Try to extract ENSG ID
            parts = orig_str.split('|')
            ensg_parts = [p for p in parts if 'ENSG' in p]
            if ensg_parts:
                clean_ensg = ensg_parts[0].split('.')[0]  # Remove version
                ref_matches = reference_mapping_df[reference_mapping_df['gene_id'] == clean_ensg]
                if not ref_matches.empty:
                    logger.info(f"    -> Could extract ENSG: {clean_ensg}")
                    could_map_count += 1
                    continue
        
        logger.info(f"    -> No mapping found")
    
    logger.info(f"\nOut of {len(sample_unmapped)} sampled unmapped IDs, {could_map_count} could potentially be mapped")
    
    return could_map_count / len(sample_unmapped) if len(sample_unmapped) > 0 else 0

def main():
    """Main analysis function."""
    parser = argparse.ArgumentParser(description='Analyze ENCODE unmapped gene IDs')
    parser.add_argument('--encode-h5ad', required=True, help='ENCODE preprocessed H5AD file')
    parser.add_argument('--encode-mapping', required=True, help='ENCODE ID to Ensembl mapping CSV')
    parser.add_argument('--reference-mapping', required=True, help='Reference gene mapping CSV')
    
    args = parser.parse_args()
    
    logger.info("Starting ENCODE unmapped gene analysis...")
    
    # Load ENCODE dataset
    if not os.path.exists(args.encode_h5ad):
        logger.error(f"ENCODE H5AD file not found: {args.encode_h5ad}")
        return 1
    
    logger.info(f"Loading ENCODE dataset from {args.encode_h5ad}")
    adata = sc.read_h5ad(args.encode_h5ad)
    logger.info(f"Loaded dataset: {adata.n_obs} samples x {adata.n_vars} genes")
    
    # Load mapping files
    logger.info(f"Loading ENCODE mapping from {args.encode_mapping}")
    encode_mapping_df = pd.read_csv(args.encode_mapping)
    
    logger.info(f"Loading reference mapping from {args.reference_mapping}")
    reference_mapping_df = pd.read_csv(args.reference_mapping)
    
    # Analyze ID formats and get unmapped mask
    unmapped_mask, id_format_stats = analyze_id_formats(adata)
    
    # Analyze expression levels
    expr_stats = analyze_expression_levels(adata, unmapped_mask)
    
    # Check mapping potential for unmapped IDs
    unmapped_original_ids = adata.var.loc[unmapped_mask, 'original_id_from_source']
    mapping_potential = check_mapping_potential(unmapped_original_ids, encode_mapping_df, reference_mapping_df)
    
    # Generate summary report
    logger.info("\n" + "="*60)
    logger.info("SUMMARY REPORT")
    logger.info("="*60)
    logger.info(f"Total genes: {adata.n_vars}")
    logger.info(f"Unmapped genes: {expr_stats['unmapped_total']} ({100*expr_stats['unmapped_total']/adata.n_vars:.1f}%)")
    logger.info(f"Unmapped genes with zero expression: {expr_stats['unmapped_zero_expr']} ({100*expr_stats['unmapped_zero_expr']/expr_stats['unmapped_total']:.1f}%)")
    logger.info(f"Unmapped genes with low expression (<0.1): {expr_stats['unmapped_low_expr']} ({100*expr_stats['unmapped_low_expr']/expr_stats['unmapped_total']:.1f}%)")
    logger.info(f"Estimated mapping potential: {100*mapping_potential:.1f}%")
    
    # Recommendation
    zero_pct = 100*expr_stats['unmapped_zero_expr']/expr_stats['unmapped_total']
    low_pct = 100*expr_stats['unmapped_low_expr']/expr_stats['unmapped_total']
    
    logger.info("\nRECOMMENDATION:")
    if zero_pct > 80:
        logger.info("✓ Most unmapped genes have zero expression - likely safe to ignore")
    elif low_pct > 80:
        logger.info("✓ Most unmapped genes have very low expression - likely safe to ignore")
    elif mapping_potential > 0.5:
        logger.info("⚠ Significant mapping potential exists - consider improving mapping logic")
    else:
        logger.info("⚠ Mixed results - manual review of unmapped IDs recommended")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())