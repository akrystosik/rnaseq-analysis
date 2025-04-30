#!/usr/bin/env python3
"""
Create Combined Dataset with All Genes (Sparse Matrix Implementation)

This script creates a combined dataset that includes all genes from all component datasets,
representing missing genes as null/NaN values rather than zeros.
It uses a sparse matrix representation for memory efficiency and ensures
consistent gene IDs and metadata across all datasets.

Usage:
    python create_combined_dataset_all_genes_sparse.py \
        --input-dir /path/to/preprocessed/datasets \
        --output-file /path/to/output/combined_all_genes_standardized.h5ad \
        --reference-mapping /path/to/gene_id_reference_mapping.csv
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
from pathlib import Path
import glob
from typing import Dict, List, Set, Tuple, Optional, Union

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('combined_dataset_creator')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Create combined dataset with all genes')
    parser.add_argument('--input-dir', type=str, required=True,
                        help='Directory containing preprocessed datasets')
    parser.add_argument('--output-file', type=str, required=True,
                        help='Output file for combined dataset')
    parser.add_argument('--reference-mapping', type=str, required=True,
                        help='Path to gene ID reference mapping CSV')
    parser.add_argument('--dataset-pattern', type=str, default='*_standardized_preprocessed.h5ad',
                        help='Pattern to match dataset files')
    parser.add_argument('--include-datasets', type=str, default='encode,gtex,mage,adni,entex',
                        help='Comma-separated list of datasets to include (e.g., "encode,gtex,mage,adni,entex")')
    parser.add_argument('--sparse', action='store_true', default=True,
                        help='Use sparse matrix representation (default: True)')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration of combined dataset even if output exists')
    
    return parser.parse_args()

def find_preprocessed_datasets(input_dir: str, pattern: str, include_datasets: List[str]) -> Dict[str, str]:
    """Find preprocessed datasets to combine based on include list."""
    logger.info(f"Looking for preprocessed datasets in {input_dir} matching pattern {pattern}")
    logger.info(f"Will include only these datasets: {', '.join(include_datasets)}")
    
    # Find matching files
    dataset_files = glob.glob(os.path.join(input_dir, pattern))
    
    # Extract dataset names and filter based on include list
    datasets = {}
    for file_path in dataset_files:
        filename = os.path.basename(file_path)
        # Skip any files with "combined" in the name
        if "combined" in filename.lower():
            logger.info(f"Skipping combined dataset: {filename}")
            continue
            
        # Extract dataset name from filename
        dataset_name = filename.split('_')[0].lower()
        
        # Only include datasets in the include list
        if dataset_name in include_datasets:
            datasets[dataset_name] = file_path
            logger.info(f"Including dataset: {dataset_name} from {filename}")
    
    logger.info(f"Found {len(datasets)} datasets to include: {', '.join(datasets.keys())}")
    return datasets

def load_reference_mapping(mapping_file: str) -> pd.DataFrame:
    """Load gene ID reference mapping from CSV file."""
    logger.info(f"Loading gene ID reference mapping from {mapping_file}")
    try:
        mapping_df = pd.read_csv(mapping_file)
        logger.info(f"Loaded reference mapping with {len(mapping_df)} entries")
        return mapping_df
    except Exception as e:
        logger.error(f"Error loading reference mapping: {e}")
        return None

def load_preprocessed_datasets(datasets: Dict[str, str]) -> Dict[str, ad.AnnData]:
    """Load preprocessed datasets."""
    loaded_datasets = {}
    
    for dataset_name, file_path in datasets.items():
        logger.info(f"Loading {dataset_name} dataset from {file_path}")
        try:
            adata = sc.read_h5ad(file_path)
            
            # Check if any samples have 'combined' as their dataset
            if 'dataset' in adata.obs.columns:
                combined_samples = sum(adata.obs['dataset'] == 'combined')
                if combined_samples > 0:
                    logger.warning(f"Found {combined_samples} samples with 'combined' as their dataset in {dataset_name}")
                    logger.warning(f"These will be renamed to their source dataset: {dataset_name}")
                    # Replace 'combined' with the actual dataset name
                    adata.obs.loc[adata.obs['dataset'] == 'combined', 'dataset'] = dataset_name
            
            loaded_datasets[dataset_name] = adata
            logger.info(f"Loaded {dataset_name} dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        except Exception as e:
            logger.error(f"Error loading {file_path}: {e}")
    
    return loaded_datasets

def create_gene_id_index(datasets: Dict[str, ad.AnnData], reference_mapping: pd.DataFrame) -> Tuple[List[str], Dict[str, int]]:
    """Create a unified gene index for the combined dataset."""
    logger.info("Creating unified gene index")
    
    # Collect all gene IDs
    all_gene_ids = set()
    
    for dataset_name, adata in datasets.items():
        # Use ensembl_id column if available, otherwise use var_names
        if 'ensembl_id' in adata.var.columns:
            for idx, row in adata.var.iterrows():
                if pd.notna(row['ensembl_id']) and row['ensembl_id'] != '':
                    all_gene_ids.add(row['ensembl_id'])
                else:
                    # If no ensembl_id, use original gene_id if it's an Ensembl ID
                    gene_id = str(idx)
                    if gene_id.startswith('ENSG') or gene_id.startswith('gSpikein'):
                        all_gene_ids.add(gene_id)
        else:
            # Use var_names directly if they're Ensembl IDs
            for gene_id in adata.var_names:
                if str(gene_id).startswith('ENSG') or str(gene_id).startswith('gSpikein'):
                    all_gene_ids.add(str(gene_id))
    
    logger.info(f"Collected {len(all_gene_ids)} unique gene IDs across all datasets")
    
    # Create sorted list of gene IDs
    # - Ensembl IDs first, sorted
    # - Spike-in controls next, sorted
    # - Other IDs last, sorted
    ensembl_ids = sorted([gene_id for gene_id in all_gene_ids if str(gene_id).startswith('ENSG')])
    spike_in_ids = sorted([gene_id for gene_id in all_gene_ids if str(gene_id).startswith('gSpikein')])
    other_ids = sorted([gene_id for gene_id in all_gene_ids 
                        if not str(gene_id).startswith('ENSG') and not str(gene_id).startswith('gSpikein')])
    
    sorted_gene_ids = ensembl_ids + spike_in_ids + other_ids
    
    # Create mapping from gene ID to index
    gene_to_idx = {gene_id: idx for idx, gene_id in enumerate(sorted_gene_ids)}
    
    logger.info(f"Created gene index with {len(sorted_gene_ids)} entries")
    logger.info(f"- Ensembl IDs: {len(ensembl_ids)}")
    logger.info(f"- Spike-in controls: {len(spike_in_ids)}")
    logger.info(f"- Other IDs: {len(other_ids)}")
    
    return sorted_gene_ids, gene_to_idx

def create_combined_sparse_matrix(datasets: Dict[str, ad.AnnData], gene_to_idx: Dict[str, int]) -> Tuple[sp.csr_matrix, List[str], List[str], pd.DataFrame]:
    """Create combined sparse matrix of gene expression data."""
    logger.info("Creating combined sparse matrix")
    
    # Calculate dimensions
    total_genes = len(gene_to_idx)
    total_samples = sum(adata.n_obs for adata in datasets.values())
    
    logger.info(f"Combined matrix dimensions: {total_samples} samples x {total_genes} genes")
    
    # Create empty sparse matrix in LIL format (good for incremental construction)
    combined_sparse = sp.lil_matrix((total_samples, total_genes), dtype=np.float32)
    
    # Track sample metadata
    sample_ids = []
    dataset_labels = []
    
    # Sample metadata as DataFrame
    obs_data = {
        'sample_id': [],
        'subject_id': [],
        'dataset': [],
        'data_type': [],
        'expression_unit': [],
        'tissue': [],
        'sex': [],
        'age': []
    }
    
    # Track which samples have which genes
    gene_presence = np.zeros((total_samples, total_genes), dtype=bool)
    
    # Track dataset of each gene
    gene_datasets = {}
    for gene_id in gene_to_idx:
        gene_datasets[gene_id] = []

    
    # Fill matrix with data from each dataset
    sample_offset = 0
    
    for dataset_name, adata in datasets.items():
        logger.info(f"Processing {dataset_name} dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Map genes from this dataset to combined matrix indices
        gene_indices = {}
        for gene_id in adata.var_names:
            # Check if this gene has an ensembl_id in var
            if 'ensembl_id' in adata.var.columns:
                ensembl_id = adata.var.loc[gene_id, 'ensembl_id']
                if pd.notna(ensembl_id) and ensembl_id != '' and ensembl_id in gene_to_idx:
                    gene_indices[gene_id] = gene_to_idx[ensembl_id]
                    gene_datasets[ensembl_id].append(dataset_name)
                elif str(gene_id).startswith('ENSG') or str(gene_id).startswith('gSpikein'):
                    # Use original ID if it's an Ensembl ID or spike-in
                    if str(gene_id) in gene_to_idx:
                        gene_indices[gene_id] = gene_to_idx[str(gene_id)]
                        gene_datasets[str(gene_id)].append(dataset_name)
            else:
                # Use gene_id directly if ensembl_id column doesn't exist
                if str(gene_id) in gene_to_idx:
                    gene_indices[gene_id] = gene_to_idx[str(gene_id)]
                    gene_datasets[str(gene_id)].append(dataset_name)
        
        # Map samples to combined matrix indices
        sample_indices = range(sample_offset, sample_offset + adata.n_obs)
        
        # Add sample metadata
        for i, (idx, row) in enumerate(adata.obs.iterrows()):
            sample_ids.append(str(idx))
            dataset_labels.append(dataset_name)
            
            # Add to obs data
            obs_data['sample_id'].append(str(idx))
            
            # Add other metadata if available
            for col in ['subject_id', 'data_type', 'expression_unit', 'tissue', 'sex', 'age']:
                if col in adata.obs.columns:
                    obs_data[col].append(row[col])
                else:
                    obs_data[col].append('')
            
            # Always add dataset
            obs_data['dataset'].append(dataset_name)
        
        # Fill expression values
        for gene_id, combined_idx in gene_indices.items():
            gene_idx = adata.var_names.get_loc(gene_id)
            
            # Get expression values for this gene
            if isinstance(adata.X, np.ndarray):
                # Dense matrix
                expr_values = adata.X[:, gene_idx]
            else:
                # Sparse matrix
                expr_values = adata.X[:, gene_idx].toarray().flatten()
            
            # Fill combined matrix
            for i, sample_idx in enumerate(sample_indices):
                combined_sparse[sample_idx, combined_idx] = expr_values[i]
                gene_presence[sample_idx, combined_idx] = True
        
        # Increment sample offset
        sample_offset += adata.n_obs
    
    # Create obs DataFrame
    obs_df = pd.DataFrame(obs_data)
    obs_df.index = sample_ids
    
    # Convert to CSR format for efficient storage and operations
    combined_sparse = combined_sparse.tocsr()
    
    # Calculate sparsity
    non_zero = combined_sparse.nnz
    total_elements = combined_sparse.shape[0] * combined_sparse.shape[1]
    sparsity = 1.0 - (non_zero / total_elements)
    
    logger.info(f"Combined matrix sparsity: {sparsity:.4f} ({non_zero} non-zero elements)")
    
    # Calculate memory savings
    dense_memory = total_elements * 4  # 4 bytes per float32
    sparse_memory = (non_zero * 4) + (non_zero * 4) + (combined_sparse.shape[0] + 1) * 4  # data + indices + indptr
    memory_ratio = dense_memory / sparse_memory
    
    logger.info(f"Memory efficiency: {memory_ratio:.2f}x (sparse: {sparse_memory/(1024*1024):.2f} MB, dense: {dense_memory/(1024*1024):.2f} MB)")
    
    return combined_sparse, sample_ids, dataset_labels, obs_df, gene_presence, gene_datasets

def create_var_dataframe(gene_ids: List[str], gene_datasets: Dict[str, List[str]], reference_mapping: pd.DataFrame) -> pd.DataFrame:
    """Create var DataFrame with gene metadata."""
    logger.info("Creating var DataFrame with gene metadata")
    
    # Create var DataFrame
    var_data = {
        'gene_id': gene_ids,
        'present_in_datasets': [],
        'chromosome': [],
        'mapping_source': [],
        'gene_name': [],
        'original_ids': [],
        'gene_type': []
    }
    
    # Create mapping from gene ID to reference row
    gene_id_to_ref = {}
    for _, row in reference_mapping.iterrows():
        gene_id = row['gene_id']
        gene_id_to_ref[gene_id] = row



    # Fill var data
    for gene_id in gene_ids:
        # Dataset presence - ensure no 'combined' values
        datasets = gene_datasets.get(gene_id, [])
        # Filter out any 'combined' values to ensure clean dataset tracking
        datasets = [ds for ds in datasets if ds != 'combined']
        var_data['present_in_datasets'].append(','.join(sorted(set(datasets))))
            
        
        # Get gene metadata from reference mapping
        if gene_id in gene_id_to_ref:
            ref_row = gene_id_to_ref[gene_id]
            var_data['chromosome'].append(ref_row['chromosome'])
            var_data['gene_name'].append(ref_row['gene_name'])
            var_data['gene_type'].append(ref_row['gene_type'])
            var_data['mapping_source'].append('reference_mapping')
            var_data['original_ids'].append(gene_id)
        else:
            # Gene not in reference mapping
            var_data['chromosome'].append('')
            var_data['gene_name'].append('')
            var_data['gene_type'].append('')
            var_data['mapping_source'].append('unmapped')
            var_data['original_ids'].append(gene_id)
    
    # Create DataFrame
    var_df = pd.DataFrame(var_data)
    var_df.index = gene_ids
    
    return var_df

def create_combined_dataset(datasets: Dict[str, ad.AnnData], reference_mapping: pd.DataFrame, output_file: str, use_sparse: bool = True) -> None:
    """Create combined dataset with all genes."""
    logger.info("Creating combined dataset with all genes")
    
    # Step 1: Create gene ID index
    gene_ids, gene_to_idx = create_gene_id_index(datasets, reference_mapping)
    
    # Step 2: Create combined expression matrix
    combined_matrix, sample_ids, dataset_labels, obs_df, gene_presence, gene_datasets = create_combined_sparse_matrix(
        datasets, gene_to_idx
    )
    
    # Step 3: Create var DataFrame
    var_df = create_var_dataframe(gene_ids, gene_datasets, reference_mapping)
    
    # Step 4: Create AnnData object
    if use_sparse:
        adata = ad.AnnData(X=combined_matrix, obs=obs_df, var=var_df)
    else:
        # Convert to dense matrix if requested
        logger.info("Converting to dense matrix")
        dense_matrix = combined_matrix.toarray()
        adata = ad.AnnData(X=dense_matrix, obs=obs_df, var=var_df)
    
    # Step 5: Add additional metadata
    adata.uns['dataset_info'] = {
        'datasets_combined': list(datasets.keys()),
        'total_samples': adata.n_obs,
        'total_genes': adata.n_vars,
        'creation_date': pd.Timestamp.now().strftime('%Y-%m-%d'),
        'approach': 'all_genes_sparse',
        'description': 'Combined dataset containing all genes from all datasets using sparse matrix representation'
    }
    
    # Add dataset overlap information
    adata.uns['dataset_overlap'] = {
        dataset_name: {
            'samples': (adata.obs['dataset'] == dataset_name).sum(),
            'genes': sum(1 for gene_id in gene_ids if dataset_name in gene_datasets[gene_id])
        } for dataset_name in datasets.keys()
    }
    
    # Add gene counts by dataset
    adata.uns['gene_counts'] = {
        dataset_name: sum(1 for gene_id in gene_ids if dataset_name in gene_datasets[gene_id])
        for dataset_name in datasets.keys()
    }
    
    # Add sparsity statistics
    non_zero = combined_matrix.nnz
    total_elements = combined_matrix.shape[0] * combined_matrix.shape[1]
    sparsity = 1.0 - (non_zero / total_elements)
    
    adata.uns['sparsity_stats'] = {
        'sparsity': sparsity,
        'non_zero_elements': int(non_zero),
        'total_elements': int(total_elements),
        'memory_saving_factor': float(total_elements / non_zero) if non_zero > 0 else 0.0
    }
    
    # Add standardization info
    adata.uns['harmonized_gencode_version'] = 'v24'
    adata.uns['harmonized_reference_genome'] = 'hg38'
    
    # Step 6: Save combined dataset
    logger.info(f"Saving combined dataset to {output_file}")
    adata.write(output_file)
    
    logger.info(f"Combined dataset saved with {adata.n_obs} samples and {adata.n_vars} genes")


def main():
    """Main function to create combined dataset with all genes."""
    args = parse_arguments()
    
    # Parse include_datasets argument
    include_datasets = [ds.strip().lower() for ds in args.include_datasets.split(',')]
    
    # Check if output already exists
    if os.path.exists(args.output_file) and not args.force:
        logger.info(f"Output file {args.output_file} already exists. Use --force to regenerate.")
        return
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
    
    # Load reference mapping
    reference_mapping = load_reference_mapping(args.reference_mapping)
    if reference_mapping is None:
        logger.error("Failed to load reference mapping")
        return
    
    # Find preprocessed datasets with filtering
    dataset_files = find_preprocessed_datasets(args.input_dir, args.dataset_pattern, include_datasets)
    if not dataset_files:
        logger.error("No preprocessed datasets found that match the inclusion criteria")
        return
    
    # Load preprocessed datasets
    datasets = load_preprocessed_datasets(dataset_files)
    if not datasets:
        logger.error("Failed to load any preprocessed datasets")
        return
    
    # Create combined dataset
    create_combined_dataset(
        datasets=datasets,
        reference_mapping=reference_mapping,
        output_file=args.output_file,
        use_sparse=args.sparse
    )
    
    logger.info("Combined dataset creation completed successfully")

if __name__ == "__main__":
    main()