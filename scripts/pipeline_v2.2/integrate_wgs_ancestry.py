#!/usr/bin/env python3
"""
Integrate WGS ancestry inference data with RNA-seq metadata.

This script distinguishes between self-reported ethnicity (cultural/social identity) 
and genetically inferred ancestry (genomic evidence) by integrating results from 
the WGS ancestry inference pipeline.

Usage:
    python integrate_wgs_ancestry.py --preprocessed-dir <dir> --wgs-ancestry-file <file> [--force]
"""

import argparse
import logging
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
import sys
import json
from typing import Dict, Optional, Tuple

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def setup_ancestry_ontology_mapping() -> Dict[str, str]:
    """
    Create mapping from WGS ancestry labels to HANCESTRO ontology terms.
    
    Returns:
        Dictionary mapping ancestry codes to HANCESTRO terms
    """
    return {
        'EUR': 'HANCESTRO:0005',  # European ancestry
        'AFR': 'HANCESTRO:0008',  # African ancestry
        'EAS': 'HANCESTRO:0018',  # East Asian ancestry
        'SAS': 'HANCESTRO:0023',  # South Asian ancestry
        'AMR': 'HANCESTRO:0004',  # Admixed American ancestry
        'unknown': 'unknown'
    }

def setup_1000g_population_mapping() -> Dict[str, str]:
    """
    Create mapping from 1000 Genomes population codes to ancestry categories.
    
    Returns:
        Dictionary mapping 1000G population codes to EUR/AFR/EAS/SAS/AMR
    """
    return {
        # European ancestry populations
        'CEU': 'EUR',  # Utah residents (CEPH) with Northern and Western European ancestry
        'TSI': 'EUR',  # Toscani in Italia
        'FIN': 'EUR',  # Finnish in Finland
        'GBR': 'EUR',  # British in England and Scotland
        'IBS': 'EUR',  # Iberian population in Spain
        
        # African ancestry populations
        'YRI': 'AFR',  # Yoruba in Ibadan, Nigeria
        'LWK': 'AFR',  # Luhya in Webuye, Kenya
        'GWD': 'AFR',  # Gambian in Western Divisions in the Gambia
        'MSL': 'AFR',  # Mende in Sierra Leone
        'ESN': 'AFR',  # Esan in Nigeria
        'ASW': 'AFR',  # Americans of African Ancestry in SW USA
        'ACB': 'AFR',  # African Caribbeans in Barbados
        
        # East Asian ancestry populations
        'CHB': 'EAS',  # Han Chinese in Beijing, China
        'JPT': 'EAS',  # Japanese in Tokyo, Japan
        'CHS': 'EAS',  # Southern Han Chinese
        'CDX': 'EAS',  # Chinese Dai in Xishuangbanna, China
        'KHV': 'EAS',  # Kinh in Ho Chi Minh City, Vietnam
        
        # South Asian ancestry populations
        'GIH': 'SAS',  # Gujarati Indians from Houston, Texas
        'PJL': 'SAS',  # Punjabi from Lahore, Pakistan
        'BEB': 'SAS',  # Bengali from Bangladesh
        'STU': 'SAS',  # Sri Lankan Tamil from the UK
        'ITU': 'SAS',  # Indian Telugu from the UK
        
        # Admixed American ancestry populations
        'MXL': 'AMR',  # Mexican Ancestry from Los Angeles USA
        'PUR': 'AMR',  # Puerto Ricans from Puerto Rico
        'CLM': 'AMR',  # Colombians from Medellin, Colombia
        'PEL': 'AMR',  # Peruvians from Lima, Peru
    }

def load_wgs_ancestry_data(wgs_ancestry_file: Path) -> pd.DataFrame:
    """
    Load WGS ancestry inference results.
    
    Args:
        wgs_ancestry_file: Path to knn_ancestry_results.csv
        
    Returns:
        DataFrame with ancestry data indexed by sample ID
    """
    logger.info(f"Loading WGS ancestry data from {wgs_ancestry_file}")
    
    if not wgs_ancestry_file.exists():
        raise FileNotFoundError(f"WGS ancestry file not found: {wgs_ancestry_file}")
    
    try:
        ancestry_df = pd.read_csv(wgs_ancestry_file)
        logger.info(f"Loaded {len(ancestry_df)} samples with ancestry data")
        
        # Validate required columns
        required_cols = ['IID', 'final_ancestry', 'confidence', 'sample_type']
        missing_cols = [col for col in required_cols if col not in ancestry_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in ancestry file: {missing_cols}")
        
        # Create mapping dictionary for efficient lookup
        ancestry_mapping = ancestry_df.set_index('IID').to_dict('index')
        
        # Log ancestry distribution
        ancestry_counts = ancestry_df['final_ancestry'].value_counts()
        logger.info("Ancestry distribution:")
        for ancestry, count in ancestry_counts.items():
            percentage = (count / len(ancestry_df)) * 100
            logger.info(f"  {ancestry}: {count} samples ({percentage:.1f}%)")
        
        # Log confidence statistics
        conf_stats = ancestry_df['confidence'].describe()
        logger.info(f"Confidence statistics: mean={conf_stats['mean']:.3f}, "
                   f"median={conf_stats['50%']:.3f}, min={conf_stats['min']:.3f}, max={conf_stats['max']:.3f}")
        
        return ancestry_mapping
        
    except Exception as e:
        logger.error(f"Error loading WGS ancestry data: {e}")
        raise

def initialize_ancestry_fields(adata) -> None:
    """
    Initialize new ancestry fields in AnnData object.
    
    Args:
        adata: AnnData object to modify
    """
    # Initialize with default values
    n_samples = adata.n_obs
    
    # Define all possible categories for ancestry fields
    ancestry_categories = ['unknown', 'EUR', 'AFR', 'EAS', 'SAS', 'AMR']
    method_categories = ['unknown', 'KNN_PCA', '1000G_Classification']
    ontology_categories = ['unknown', 'HANCESTRO:0005', 'HANCESTRO:0008', 'HANCESTRO:0018', 'HANCESTRO:0023', 'HANCESTRO:0004']
    
    adata.obs['inferred_ancestry'] = pd.Categorical(['unknown'] * n_samples, categories=ancestry_categories)
    adata.obs['inferred_ancestry_confidence'] = np.full(n_samples, np.nan, dtype=float)
    adata.obs['inferred_ancestry_method'] = pd.Categorical(['unknown'] * n_samples, categories=method_categories)
    adata.obs['inferred_ancestry_ontology_term_id'] = pd.Categorical(['unknown'] * n_samples, categories=ontology_categories)
    adata.obs['has_genetic_ancestry_data'] = np.full(n_samples, False, dtype=bool)
    
    logger.info(f"Initialized ancestry fields for {n_samples} samples")

def integrate_ancestry_data(adata, ancestry_mapping: Dict, ontology_mapping: Dict) -> Tuple[int, int]:
    """
    Integrate ancestry data into AnnData object.
    
    Args:
        adata: AnnData object to update
        ancestry_mapping: Dictionary mapping sample IDs to ancestry info
        ontology_mapping: Dictionary mapping ancestry codes to HANCESTRO terms
        
    Returns:
        Tuple of (matched_samples, total_samples)
    """
    matched_samples = 0
    total_samples = adata.n_obs
    
    # Get sample ID column name (may vary by dataset)
    sample_id_col = None
    for col in ['sample_id', 'donor_id', 'subject_id']:
        if col in adata.obs.columns:
            sample_id_col = col
            break
    
    if sample_id_col is None:
        raise ValueError("No suitable sample ID column found in AnnData object")
    
    logger.info(f"Using '{sample_id_col}' column for sample ID matching")
    
    # Integrate ancestry data
    for idx, row in adata.obs.iterrows():
        sample_id = row[sample_id_col]
        
        if sample_id in ancestry_mapping:
            ancestry_info = ancestry_mapping[sample_id]
            
            # Update ancestry fields
            adata.obs.loc[idx, 'inferred_ancestry'] = ancestry_info['final_ancestry']
            adata.obs.loc[idx, 'inferred_ancestry_confidence'] = ancestry_info['confidence']
            adata.obs.loc[idx, 'inferred_ancestry_method'] = 'KNN_PCA'
            adata.obs.loc[idx, 'inferred_ancestry_ontology_term_id'] = ontology_mapping.get(
                ancestry_info['final_ancestry'], 'unknown'
            )
            adata.obs.loc[idx, 'has_genetic_ancestry_data'] = True
            
            matched_samples += 1
    
    return matched_samples, total_samples

def generate_integration_summary(matched_samples: int, total_samples: int, dataset_name: str) -> Dict:
    """
    Generate summary statistics for ancestry integration.
    
    Args:
        matched_samples: Number of samples with ancestry data
        total_samples: Total number of samples in dataset
        dataset_name: Name of the dataset
        
    Returns:
        Dictionary with summary statistics
    """
    integration_rate = (matched_samples / total_samples) * 100 if total_samples > 0 else 0
    
    summary = {
        'dataset': dataset_name,
        'total_samples': int(total_samples),
        'matched_samples': int(matched_samples),
        'integration_rate_percent': round(integration_rate, 2),
        'unmatched_samples': int(total_samples - matched_samples)
    }
    
    logger.info(f"{dataset_name}: {matched_samples}/{total_samples} samples matched "
               f"({integration_rate:.1f}% integration rate)")
    
    return summary

def analyze_concordance(adata, dataset_name: str = None) -> Dict:
    """
    Analyze concordance between self-reported ethnicity and inferred ancestry.
    
    Args:
        adata: AnnData object with both ethnicity and ancestry data
        dataset_name: Name of the dataset (used to exclude inappropriate concordance analysis)
        
    Returns:
        Dictionary with concordance statistics
    """
    # Skip concordance analysis for MAGE - ethnicity labels are derived from genetic data, not self-reported
    if dataset_name == 'mage':
        logger.info("Skipping concordance analysis for MAGE: ethnicity labels are derived from genetic classifications")
        return {'concordance_samples': 0, 'note': 'concordance_not_applicable_mage'}
    
    # Only analyze samples with both self-reported ethnicity and inferred ancestry
    has_both = (
        adata.obs['has_genetic_ancestry_data'] & 
        (adata.obs['self_reported_ethnicity'] != 'unknown') &
        (adata.obs['self_reported_ethnicity'].notna())
    )
    
    if not has_both.any():
        logger.warning("No samples have both self-reported ethnicity and inferred ancestry data")
        return {'concordance_samples': 0}
    
    analysis_subset = adata.obs[has_both]
    n_analysis = len(analysis_subset)
    
    # Create simplified mapping for concordance analysis
    ethnicity_to_ancestry = {
        'European': 'EUR', 'White': 'EUR', 'white': 'EUR',
        'African American': 'AFR', 'Black or African American': 'AFR', 'African': 'AFR',
        'Asian': 'EAS', 'East Asian': 'EAS',
        'South Asian': 'SAS', 'Indian': 'SAS',
        'Hispanic or Latino': 'AMR', 'American Indian or Alaska Native': 'AMR'
    }
    
    # Map self-reported ethnicity to expected ancestry
    analysis_subset = analysis_subset.copy()
    analysis_subset['expected_ancestry'] = analysis_subset['self_reported_ethnicity'].map(
        ethnicity_to_ancestry
    ).fillna('unknown')
    
    # Calculate concordance
    concordant = (
        analysis_subset['expected_ancestry'] == analysis_subset['inferred_ancestry']
    )
    
    concordance_rate = concordant.sum() / len(analysis_subset) * 100 if len(analysis_subset) > 0 else 0
    
    concordance_stats = {
        'concordance_samples': int(n_analysis),
        'concordant_samples': int(concordant.sum()),
        'discordant_samples': int((~concordant).sum()),
        'concordance_rate_percent': round(concordance_rate, 2)
    }
    
    logger.info(f"Concordance analysis: {concordant.sum()}/{n_analysis} samples concordant "
               f"({concordance_rate:.1f}%)")
    
    return concordance_stats

def integrate_mage_1000g_ancestry(adata, population_1000g_mapping: Dict, ontology_mapping: Dict) -> Tuple[int, int]:
    """
    Integrate MAGE 1000 Genomes population codes as genetic ancestry data.
    
    Args:
        adata: AnnData object to update
        population_1000g_mapping: Dictionary mapping 1000G codes to ancestry categories
        ontology_mapping: Dictionary mapping ancestry codes to HANCESTRO terms
        
    Returns:
        Tuple of (matched_samples, total_samples)
    """
    matched_samples = 0
    total_samples = adata.n_obs
    
    if 'population_code_1000g' not in adata.obs.columns:
        logger.warning("No population_code_1000g column found in MAGE dataset")
        return 0, total_samples
    
    logger.info("Integrating MAGE 1000 Genomes population codes as genetic ancestry")
    
    # Integrate 1000G population codes as ancestry data
    for idx, row in adata.obs.iterrows():
        pop_code = row['population_code_1000g']
        
        if pd.notna(pop_code) and pop_code in population_1000g_mapping:
            ancestry_code = population_1000g_mapping[pop_code]
            
            # Update ancestry fields
            adata.obs.loc[idx, 'inferred_ancestry'] = ancestry_code
            adata.obs.loc[idx, 'inferred_ancestry_confidence'] = 1.0  # High confidence for 1000G classifications
            adata.obs.loc[idx, 'inferred_ancestry_method'] = '1000G_Classification'
            adata.obs.loc[idx, 'inferred_ancestry_ontology_term_id'] = ontology_mapping.get(ancestry_code, 'unknown')
            adata.obs.loc[idx, 'has_genetic_ancestry_data'] = True
            
            matched_samples += 1
    
    return matched_samples, total_samples

def process_dataset(h5ad_file: Path, ancestry_mapping: Dict, ontology_mapping: Dict, 
                   force: bool = False) -> Dict:
    """
    Process a single H5AD dataset file.
    
    Args:
        h5ad_file: Path to H5AD file
        ancestry_mapping: Dictionary mapping sample IDs to ancestry info
        ontology_mapping: Dictionary mapping ancestry codes to HANCESTRO terms
        force: Whether to overwrite existing ancestry data
        
    Returns:
        Dictionary with processing summary
    """
    dataset_name = h5ad_file.stem.replace('_standardized_preprocessed', '')
    logger.info(f"Processing dataset: {dataset_name}")
    
    try:
        # Load AnnData object
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded {adata.n_obs} samples, {adata.n_vars} genes")
        
        # Check if ancestry data already exists
        has_ancestry_data = 'has_genetic_ancestry_data' in adata.obs.columns
        if has_ancestry_data and not force:
            existing_count = adata.obs['has_genetic_ancestry_data'].sum()
            logger.info(f"Ancestry data already exists for {existing_count} samples. "
                       f"Use --force to overwrite.")
            return {
                'dataset': dataset_name,
                'status': 'skipped',
                'reason': 'ancestry_data_exists',
                'existing_samples': int(existing_count)
            }
        
        # Initialize ancestry fields
        initialize_ancestry_fields(adata)
        
        # Handle MAGE dataset specially using 1000G population codes
        # Note: MAGE samples were used as reference populations for KNN ancestry inference,
        # so they should use their original 1000G classifications rather than KNN predictions
        if dataset_name == 'mage':
            population_1000g_mapping = setup_1000g_population_mapping()
            matched_samples, total_samples = integrate_mage_1000g_ancestry(
                adata, population_1000g_mapping, ontology_mapping
            )
        else:
            # Integrate WGS ancestry data for other datasets
            matched_samples, total_samples = integrate_ancestry_data(
                adata, ancestry_mapping, ontology_mapping
            )
        
        # Generate summary
        summary = generate_integration_summary(matched_samples, total_samples, dataset_name)
        
        # Analyze concordance if possible
        if 'self_reported_ethnicity' in adata.obs.columns:
            concordance_stats = analyze_concordance(adata, dataset_name)
            summary.update(concordance_stats)
        
        # Save updated H5AD file
        adata.write_h5ad(h5ad_file)
        logger.info(f"Updated H5AD file saved: {h5ad_file}")
        
        summary['status'] = 'completed'
        return summary
        
    except Exception as e:
        logger.error(f"Error processing {h5ad_file}: {e}")
        return {
            'dataset': dataset_name,
            'status': 'error',
            'error_message': str(e)
        }

def main():
    """Main function to integrate WGS ancestry data with RNA-seq metadata."""
    parser = argparse.ArgumentParser(
        description="Integrate WGS ancestry data with RNA-seq metadata"
    )
    parser.add_argument(
        '--preprocessed-dir', 
        type=Path, 
        required=True,
        help='Directory containing preprocessed H5AD files'
    )
    parser.add_argument(
        '--wgs-ancestry-file', 
        type=Path, 
        required=True,
        help='Path to WGS ancestry results CSV file'
    )
    parser.add_argument(
        '--force', 
        action='store_true',
        help='Overwrite existing ancestry data'
    )
    parser.add_argument(
        '--output-summary', 
        type=Path,
        help='Path to save integration summary JSON'
    )
    
    args = parser.parse_args()
    
    # Validate input paths
    if not args.preprocessed_dir.exists():
        logger.error(f"Preprocessed directory not found: {args.preprocessed_dir}")
        sys.exit(1)
    
    if not args.wgs_ancestry_file.exists():
        logger.error(f"WGS ancestry file not found: {args.wgs_ancestry_file}")
        sys.exit(1)
    
    try:
        # Load WGS ancestry data
        ancestry_mapping = load_wgs_ancestry_data(args.wgs_ancestry_file)
        ontology_mapping = setup_ancestry_ontology_mapping()
        
        # Find H5AD files to process
        h5ad_files = list(args.preprocessed_dir.glob("*_standardized_preprocessed.h5ad"))
        if not h5ad_files:
            logger.warning(f"No H5AD files found in {args.preprocessed_dir}")
            return
        
        logger.info(f"Found {len(h5ad_files)} H5AD files to process")
        
        # Process each dataset
        summaries = []
        for h5ad_file in h5ad_files:
            summary = process_dataset(h5ad_file, ancestry_mapping, ontology_mapping, args.force)
            summaries.append(summary)
        
        # Generate overall summary
        total_samples = sum(s.get('total_samples', 0) for s in summaries if s.get('status') == 'completed')
        total_matched = sum(s.get('matched_samples', 0) for s in summaries if s.get('status') == 'completed')
        overall_rate = (total_matched / total_samples * 100) if total_samples > 0 else 0
        
        logger.info("=" * 60)
        logger.info("ANCESTRY INTEGRATION SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Total samples processed: {total_samples:,}")
        logger.info(f"Total samples with ancestry: {total_matched:,}")
        logger.info(f"Overall integration rate: {overall_rate:.1%}")
        
        # Save summary if requested
        if args.output_summary:
            summary_data = {
                'total_samples': int(total_samples),
                'total_matched': int(total_matched),
                'overall_integration_rate_percent': round(overall_rate, 2),
                'dataset_summaries': summaries,
                'ancestry_mapping_source': str(args.wgs_ancestry_file),
                'processing_timestamp': pd.Timestamp.now().isoformat()
            }
            
            with open(args.output_summary, 'w') as f:
                json.dump(summary_data, f, indent=2)
            logger.info(f"Integration summary saved to: {args.output_summary}")
        
        logger.info("Ancestry integration completed successfully!")
        
    except Exception as e:
        logger.error(f"Fatal error during ancestry integration: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()