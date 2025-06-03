#!/usr/bin/env python3
"""
Dataset Gene ID Preprocessing Script

This script preprocesses datasets to ensure consistent gene identifiers.
It applies the reference mapping to standardize gene IDs across all datasets,
adds proper var columns, and includes comprehensive gene metadata.

Usage:
    python preprocess_dataset_gene_ids.py \
        --data-dir /path/to/standardized/data \
        --reference-mapping /path/to/gene_id_reference_mapping.csv \
        --output-dir /path/to/output/directory
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
import glob
import json
from typing import Dict, List, Set, Tuple, Optional


def ensure_string_columns(df):
    """Convert categorical columns to strings to avoid comparison issues."""
    for col in df.columns:
        if pd.api.types.is_categorical_dtype(df[col]):
            logger.debug(f"Converting categorical column {col} to string")
            df[col] = df[col].astype(str)
    return df

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('dataset_gene_id_preprocessor')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Preprocess dataset gene IDs')
    parser.add_argument('--data-dir', type=str, required=True,
                        help='Directory containing standardized datasets')
    parser.add_argument('--reference-mapping', type=str, required=True,
                        help='Path to gene ID reference mapping CSV')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Directory to save preprocessed datasets')
    parser.add_argument('--datasets', type=str, default='all',
                        help='Comma-separated list of datasets to process (or "all")')
    parser.add_argument('--force', action='store_true',
                        help='Force preprocessing even if output exists')
    
    return parser.parse_args()

def load_reference_mapping(mapping_file: str) -> pd.DataFrame:
    """Load gene ID reference mapping from CSV file."""
    logger.info(f"Loading gene ID reference mapping from {mapping_file}")
    try:
        mapping_df = pd.read_csv(mapping_file)
        logger.info(f"Loaded reference mapping with {len(mapping_df)} entries")
        
        # Create specialized mappings for quicker lookups
        numeric_to_ensembl = {}
        ensembl_to_info = {}
        
        for _, row in mapping_df.iterrows():
            gene_id = row['gene_id']
            numeric_id = row['numeric_id']
            
            # Map numeric IDs to Ensembl IDs
            if not pd.isna(numeric_id):
                numeric_to_ensembl[str(numeric_id)] = gene_id
            
            # Map Ensembl IDs to full gene info
            ensembl_to_info[gene_id] = {
                'gene_name': row['gene_name'],
                'gene_type': row['gene_type'],
                'chromosome': row['chromosome'],
                'mapping_confidence': row['mapping_confidence']
            }
        
        # Log some statistics about the mappings
        logger.info(f"Created numeric ID mapping with {len(numeric_to_ensembl)} entries")
        logger.info(f"Created Ensembl ID mapping with {len(ensembl_to_info)} entries")
        
        return mapping_df, numeric_to_ensembl, ensembl_to_info
    
    except Exception as e:
        logger.error(f"Error loading reference mapping: {e}")
        return None, {}, {}

def find_datasets(data_dir: str, datasets_list: str) -> Dict[str, str]:
    """Find datasets to process."""
    logger.info(f"Looking for datasets in {data_dir}")
    
    # Find all h5ad files
    h5ad_files = glob.glob(os.path.join(data_dir, '*_standardized*.h5ad'))
    
    # Extract dataset names
    datasets = {}
    for file_path in h5ad_files:
        dataset_name = os.path.basename(file_path).split('_')[0].lower()
        datasets[dataset_name] = file_path
    
    logger.info(f"Found {len(datasets)} datasets: {', '.join(datasets.keys())}")
    
    # Filter datasets if needed
    if datasets_list != 'all':
        selected_datasets = {}
        for dataset in datasets_list.split(','):
            dataset = dataset.strip().lower()
            if dataset in datasets:
                selected_datasets[dataset] = datasets[dataset]
            else:
                logger.warning(f"Requested dataset {dataset} not found")
        
        datasets = selected_datasets
        logger.info(f"Selected {len(datasets)} datasets for processing")
    
    return datasets

def preprocess_encode_dataset(adata: ad.AnnData, numeric_to_ensembl: Dict[str, str], ensembl_to_info: Dict[str, Dict], output_file: str) -> None:
    """Preprocess ENCODE dataset to add gene metadata and consistent IDs."""
    logger.info(f"Preprocessing ENCODE dataset with {adata.n_vars} genes")
    
    # Load the enhanced ENCODE ID mapping
    encode_mapping_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/encode_id_to_ensembl_mapping.csv'
    try:
        encode_mapping = pd.read_csv(encode_mapping_file)
        # Create a mapping dictionary from original IDs to Ensembl IDs
        encode_to_ensembl = {}
        for _, row in encode_mapping.iterrows():
            original_id = str(row['original_id'])
            ensembl_id = row['ensembl_id']
            
            # Skip ENTREZ: prefix for numeric IDs
            if ensembl_id.startswith('ENTREZ:'):
                entrez_id = ensembl_id[7:]  # Remove ENTREZ: prefix
                # Try to map using our enhanced Entrez mapping
                if entrez_id in numeric_to_ensembl:
                    ensembl_id = numeric_to_ensembl[entrez_id]
                else:
                    # Keep as is, will try to map later
                    ensembl_id = entrez_id
            
            encode_to_ensembl[original_id] = ensembl_id
        
        logger.info(f"Loaded ENCODE ID mapping with {len(encode_to_ensembl)} entries")
    except Exception as e:
        logger.error(f"Error loading ENCODE ID mapping: {e}")
        encode_to_ensembl = {}
    
    # Load the enhanced Entrez-to-Ensembl mapping as fallback
    enhanced_mapping_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/enhanced_entrez_to_ensembl_mapping.csv'
    try:
        enhanced_mapping = pd.read_csv(enhanced_mapping_file)
        # Create a direct mapping dictionary
        enhanced_entrez_to_ensembl = {}
        for _, row in enhanced_mapping.iterrows():
            entrez_id = str(row['entrez_id'])
            ensembl_id = row['ensembl_id']
            enhanced_entrez_to_ensembl[entrez_id] = ensembl_id
        
        logger.info(f"Loaded enhanced Entrez-to-Ensembl mapping with {len(enhanced_entrez_to_ensembl)} entries")
    except Exception as e:
        logger.error(f"Error loading enhanced Entrez-to-Ensembl mapping: {e}")
        enhanced_entrez_to_ensembl = {}
    
    # Create var DataFrame with mapping information
    var_columns = {
        'gene_id': [],             # Original ID 
        'original_gene_id': [],    # Original ID preserved
        'ensembl_id': [],          # Standard Ensembl ID (without version)
        'gene_name': [],           # Gene symbol
        'gene_type': [],           # Gene type (protein_coding, etc.)
        'chromosome': [],          # Chromosome
        'mapping_source': [],      # Source of mapping (GENCODE, etc.)
        'mapping_confidence': []   # Confidence of mapping (high, medium, low)
    }
    
    # Counters for tracking mapping results
    mapped_from_encode = 0
    mapped_from_enhanced = 0
    mapped_from_reference = 0
    unmapped_count = 0
    spikein_count = 0
    
    # Attempt to retrieve original ENCODE IDs when possible
    # This is challenging because the standardized dataset uses sequential indices
    # We'll try to use the original ID from ENCODE files if it matches in length
    
    # Get original ENCODE IDs
    original_ids = list(encode_to_ensembl.keys())
    seq_indices = [str(i) for i in range(adata.n_vars)]
    
    # Create a mapping table if lengths match
    original_id_map = {}
    if len(original_ids) >= adata.n_vars:
        # For each sequential index, map to a potential original ID
        for i, seq_idx in enumerate(seq_indices):
            if i < len(original_ids):
                original_id_map[seq_idx] = original_ids[i]
    
    # Process each gene in the dataset
    for gene_idx in adata.var_names:
        # Store sequential index as gene_id
        var_columns['gene_id'].append(gene_idx)
        
        # Try to map to original ID if possible
        original_id = original_id_map.get(str(gene_idx), gene_idx)
        var_columns['original_gene_id'].append(original_id)
        
        # Map to Ensembl ID using ENCODE mapping
        ensembl_id = ''
        mapping_source = ''
        
        # Try ENCODE mapping first
        if str(original_id) in encode_to_ensembl:
            ensembl_id = encode_to_ensembl[str(original_id)]
            if ensembl_id.startswith('gSpikein'):
                mapping_source = 'spike_in'
                spikein_count += 1
            else:
                mapping_source = 'encode_mapping'
                mapped_from_encode += 1
        
        # If not found or it's a numeric ID, try enhanced mapping
        elif str(original_id).isdigit() and str(original_id) in enhanced_entrez_to_ensembl:
            ensembl_id = enhanced_entrez_to_ensembl[str(original_id)]
            mapping_source = 'enhanced_mapping'
            mapped_from_enhanced += 1
        
        # Try reference mapping as last resort
        elif str(original_id).isdigit() and str(original_id) in numeric_to_ensembl:
            ensembl_id = numeric_to_ensembl[str(original_id)]
            # Skip placeholder IDs and keep original numeric ID
            if ensembl_id.startswith('PLACEHOLDER_'):
                # Use original ID as Ensembl ID to avoid empty values
                ensembl_id = f"ENTREZ:{str(original_id)}"
                mapping_source = 'entrez_id'
                mapped_from_reference += 1
            else:
                mapping_source = 'reference_mapping'
                mapped_from_reference += 1
        else:
            mapping_source = 'unmapped'
            unmapped_count += 1
        
        var_columns['ensembl_id'].append(ensembl_id)
        var_columns['mapping_source'].append(mapping_source)
        
        # Add gene information if Ensembl ID is available
        if ensembl_id and ensembl_id in ensembl_to_info:
            info = ensembl_to_info[ensembl_id]
            var_columns['gene_name'].append(info['gene_name'])
            var_columns['gene_type'].append(info['gene_type'])
            var_columns['chromosome'].append(info['chromosome'])
            var_columns['mapping_confidence'].append(info['mapping_confidence'])
        elif ensembl_id.startswith('gSpikein'):
            # Handle spike-in controls
            var_columns['gene_name'].append(ensembl_id)
            var_columns['gene_type'].append('spike_in_control')
            var_columns['chromosome'].append('spike_in')
            var_columns['mapping_confidence'].append('high')
        else:
            # Unmapped genes
            var_columns['gene_name'].append('')
            var_columns['gene_type'].append('unknown')
            var_columns['chromosome'].append('')
            var_columns['mapping_confidence'].append('none')
    
    # Create new var DataFrame
    new_var = pd.DataFrame(var_columns)
    
    # Fix gene_id column to use ensembl_id (or ENTREZ: id) instead of sequential index
    for i, row in new_var.iterrows():
        if row['ensembl_id']:
            # Use the ensembl_id as the gene_id
            new_var.at[i, 'gene_id'] = row['ensembl_id']
        elif row['original_gene_id'].startswith('gSpikein'):
            # For spike-in controls, use the original ID
            new_var.at[i, 'gene_id'] = row['original_gene_id']
    
    # Calculate mapping statistics
    mapped_count = mapped_from_encode + mapped_from_enhanced + mapped_from_reference + spikein_count
    mapping_percentage = mapped_count / adata.n_vars * 100
    
    logger.info(f"ENCODE mapping details:")
    logger.info(f"  - Mapped from ENCODE mapping: {mapped_from_encode}")
    logger.info(f"  - Mapped from enhanced mapping: {mapped_from_enhanced}")
    logger.info(f"  - Mapped from reference mapping: {mapped_from_reference}")
    logger.info(f"  - Spike-in controls: {spikein_count}")
    logger.info(f"  - Unmapped: {unmapped_count}")
    logger.info(f"Total mapped: {mapped_count}/{adata.n_vars} genes ({mapping_percentage:.2f}%)")
    
    # Replace var DataFrame
    adata.var = new_var
    
    # Store additional metadata in uns
    adata.uns['gene_mapping_stats'] = {
        'mapped_genes': mapped_count,
        'mapped_from_encode': mapped_from_encode,
        'mapped_from_enhanced': mapped_from_enhanced,
        'mapped_from_reference': mapped_from_reference,
        'spikein_count': spikein_count,
        'unmapped_count': unmapped_count,
        'total_genes': adata.n_vars,
        'mapping_percentage': mapping_percentage
    }
    
    # Save preprocessed dataset
    logger.info(f"Saving preprocessed ENCODE dataset to {output_file}")
    adata.write(output_file)
    
def preprocess_entex_dataset(adata: ad.AnnData, numeric_to_ensembl: Dict[str, str], ensembl_to_info: Dict[str, Dict], output_file: str) -> None:
    """Preprocess ENTEx dataset to handle mixed ID formats."""
    logger.info(f"Preprocessing ENTEx dataset with {adata.n_vars} genes")
    
    # Create var DataFrame with mapping information
    var_columns = {
        'gene_id': [],             # Original ID 
        'original_gene_id': [],    # Original Ensembl ID with version
        'ensembl_id': [],          # Standard Ensembl ID
        'gene_name': [],           # Gene symbol
        'gene_type': [],           # Gene type (protein_coding, etc.)
        'chromosome': [],          # Chromosome
        'mapping_source': [],      # Source of mapping (GENCODE, etc.)
        'mapping_confidence': []   # Confidence of mapping (high, medium, low)
    }
    
    for gene_id in adata.var_names:
        # Store original ID
        var_columns['gene_id'].append(gene_id)
        var_columns['original_gene_id'].append(gene_id)  # Always preserve original
        
        # Handle different ID formats
        if str(gene_id).startswith('ENSG'):
            # Already an Ensembl ID
            ensembl_id = str(gene_id).split('.')[0]  # Remove version for mapping
            var_columns['ensembl_id'].append(ensembl_id)
            
            # Look up gene info
            if ensembl_id in ensembl_to_info:
                info = ensembl_to_info[ensembl_id]
                var_columns['gene_name'].append(info['gene_name'])
                var_columns['gene_type'].append(info['gene_type'])
                var_columns['chromosome'].append(info['chromosome'])
                var_columns['mapping_confidence'].append(info['mapping_confidence'])
                var_columns['mapping_source'].append('reference_mapping')
            else:
                # Ensembl ID not in reference mapping
                var_columns['gene_name'].append('')
                var_columns['gene_type'].append('unknown')
                var_columns['chromosome'].append('')
                var_columns['mapping_confidence'].append('medium')
                var_columns['mapping_source'].append('original_ensembl')
        
        elif str(gene_id).isdigit():
            # Numeric ID - look up Ensembl ID
            ensembl_id = numeric_to_ensembl.get(str(gene_id), '')
            var_columns['ensembl_id'].append(ensembl_id)
            
            # Look up gene info if we have an Ensembl ID
            if ensembl_id and ensembl_id in ensembl_to_info:
                info = ensembl_to_info[ensembl_id]
                var_columns['gene_name'].append(info['gene_name'])
                var_columns['gene_type'].append(info['gene_type'])
                var_columns['chromosome'].append(info['chromosome'])
                var_columns['mapping_confidence'].append(info['mapping_confidence'])
                var_columns['mapping_source'].append('reference_mapping')
            elif ensembl_id.startswith('PLACEHOLDER_'):
                # Handle placeholders
                var_columns['gene_name'].append('')
                var_columns['gene_type'].append('unknown')
                var_columns['chromosome'].append('')
                var_columns['mapping_confidence'].append('low')
                var_columns['mapping_source'].append('placeholder')
            else:
                # Unmapped genes
                var_columns['gene_name'].append('')
                var_columns['gene_type'].append('unknown')
                var_columns['chromosome'].append('')
                var_columns['mapping_confidence'].append('none')
                var_columns['mapping_source'].append('unmapped')
        
        elif str(gene_id).startswith('gSpikein'):
            # Spike-in control
            var_columns['ensembl_id'].append(str(gene_id))
            var_columns['gene_name'].append(str(gene_id))
            var_columns['gene_type'].append('spike_in_control')
            var_columns['chromosome'].append('spike_in')
            var_columns['mapping_confidence'].append('high')
            var_columns['mapping_source'].append('spike_in')
        
        else:
            # Other ID formats
            var_columns['ensembl_id'].append('')
            var_columns['gene_name'].append('')
            var_columns['gene_type'].append('unknown')
            var_columns['chromosome'].append('')
            var_columns['mapping_confidence'].append('none')
            var_columns['mapping_source'].append('unmapped')
    
    # Create new var DataFrame
    new_var = pd.DataFrame(var_columns)
    
    # Fix gene_id column to use ensembl_id (or ENTREZ: id) instead of sequential index
    for i, row in new_var.iterrows():
        if row['ensembl_id']:
            # Use the ensembl_id as the gene_id
            new_var.at[i, 'gene_id'] = row['ensembl_id']
        elif row['original_gene_id'].startswith('gSpikein'):
            # For spike-in controls, use the original ID
            new_var.at[i, 'gene_id'] = row['original_gene_id']
    
    # Calculate mapping statistics
    mapped_count = sum(1 for x in var_columns['ensembl_id'] if x and not x.startswith('PLACEHOLDER_'))
    mapping_percentage = mapped_count / adata.n_vars * 100
    
    logger.info(f"Mapped {mapped_count}/{adata.n_vars} genes ({mapping_percentage:.2f}%)")
    
    # Replace var DataFrame
    adata.var = new_var
    
    # Store additional metadata in uns
    adata.uns['gene_mapping_stats'] = {
        'mapped_genes': mapped_count,
        'total_genes': adata.n_vars,
        'mapping_percentage': mapping_percentage,
        'mapping_source': 'reference_mapping'
    }
    
    # Save preprocessed dataset
    logger.info(f"Saving preprocessed ENTEx dataset to {output_file}")
    adata.write(output_file)

def preprocess_other_dataset(adata: ad.AnnData, dataset_name: str, ensembl_to_info: Dict[str, Dict], output_file: str) -> None:
    """Preprocess other datasets (ADNI, GTEx, MAGE) with standard Ensembl IDs."""
    logger.info(f"Preprocessing {dataset_name} dataset with {adata.n_vars} genes")
    
    # Check if var already has necessary columns
    required_columns = ['gene_id', 'original_gene_id', 'ensembl_id', 'gene_name', 'gene_type', 'chromosome']
    missing_columns = [col for col in required_columns if col not in adata.var.columns]
    
    # If the existing var DataFrame already has all columns
    if not missing_columns:
        logger.info(f"{dataset_name} dataset already has all required columns")
        
        # No need to add columns, but ensure they're populated correctly
        pass
    else:
        logger.info(f"Adding missing columns to {dataset_name} dataset: {missing_columns}")
        
        # Create a new DataFrame with all required columns
        var_data = {}
        
        # Start with existing columns
        for col in adata.var.columns:
            var_data[col] = adata.var[col].tolist()
        
        # Add missing columns
        if 'gene_id' not in var_data:
            var_data['gene_id'] = adata.var_names.tolist()
        
        if 'original_gene_id' not in var_data:
            var_data['original_gene_id'] = adata.var_names.tolist()
        
        if 'ensembl_id' not in var_data:
            # For datasets like ADNI, GTEx, MAGE that already use Ensembl IDs
            # Use the var_names but strip version numbers
            var_data['ensembl_id'] = [str(id).split('.')[0] for id in adata.var_names]
        
        # Add other missing columns as needed
        for col in ['gene_name', 'gene_type', 'chromosome', 'mapping_source', 'mapping_confidence']:
            if col not in var_data:
                var_data[col] = [''] * adata.n_vars
        
        # Create the new var DataFrame
        new_var = pd.DataFrame(var_data, index=adata.var_names)
        
        # Update gene information where available
        for i, gene_id in enumerate(adata.var_names):
            ensembl_id = str(gene_id).split('.')[0]  # Remove version if present
            
            if ensembl_id in ensembl_to_info:
                info = ensembl_to_info[ensembl_id]
                new_var.loc[gene_id, 'gene_name'] = info['gene_name']
                new_var.loc[gene_id, 'gene_type'] = info['gene_type']
                new_var.loc[gene_id, 'chromosome'] = info['chromosome']
                new_var.loc[gene_id, 'mapping_confidence'] = info['mapping_confidence']
                new_var.loc[gene_id, 'mapping_source'] = 'reference_mapping'
            else:
                new_var.loc[gene_id, 'mapping_source'] = 'unmapped'
                new_var.loc[gene_id, 'mapping_confidence'] = 'none'
        
        # Replace the existing var DataFrame
        adata.var = new_var
    
    # Calculate mapping statistics
    mapped_count = sum(1 for gene_id in adata.var_names if str(gene_id).split('.')[0] in ensembl_to_info)
    mapping_percentage = mapped_count / adata.n_vars * 100
    
    logger.info(f"Mapped {mapped_count}/{adata.n_vars} genes ({mapping_percentage:.2f}%)")
    
    # Store additional metadata in uns
    adata.uns['gene_mapping_stats'] = {
        'mapped_genes': mapped_count,
        'total_genes': adata.n_vars,
        'mapping_percentage': mapping_percentage,
        'mapping_source': 'reference_mapping'
    }
    
    # Save preprocessed dataset
    logger.info(f"Saving preprocessed {dataset_name} dataset to {output_file}")
    adata.write(output_file)
    
def main():
    """Main function to preprocess dataset gene IDs."""
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load reference mapping
    mapping_df, numeric_to_ensembl, ensembl_to_info = load_reference_mapping(args.reference_mapping)
    if mapping_df is None:
        logger.error("Failed to load reference mapping")
        return
    
    # Find datasets to process
    datasets = find_datasets(args.data_dir, args.datasets)
    if not datasets:
        logger.error("No datasets found to process")
        return
    
    # Process each dataset
    for dataset_name, file_path in datasets.items():
        output_file = os.path.join(args.output_dir, f"{dataset_name}_standardized_preprocessed.h5ad")
        
        # Check if output already exists
        if os.path.exists(output_file) and not args.force:
            logger.info(f"Output file {output_file} already exists. Use --force to regenerate.")
            continue
        
        # Load dataset
        logger.info(f"Loading {dataset_name} dataset from {file_path}")
        try:
            adata = sc.read_h5ad(file_path)
            logger.info(f"Loaded {dataset_name} dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        except Exception as e:
            logger.error(f"Error loading {file_path}: {e}")
            continue
        
        # Preprocess dataset based on type
        if dataset_name.lower() == 'encode':
            preprocess_encode_dataset(adata, numeric_to_ensembl, ensembl_to_info, output_file)
        elif dataset_name.lower() == 'entex':
            preprocess_entex_dataset(adata, numeric_to_ensembl, ensembl_to_info, output_file)
        else:
            # ADNI, GTEx, MAGE
            preprocess_other_dataset(adata, dataset_name, ensembl_to_info, output_file)
    
    logger.info("Preprocessing completed successfully")

if __name__ == "__main__":
    main()