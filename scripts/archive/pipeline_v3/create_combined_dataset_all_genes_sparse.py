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
import scipy.sparse as sp # Ensure scipy.sparse is imported
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
            
            # Ensure all samples have the correct dataset label
            if 'dataset' in adata.obs.columns:
                # Set all samples to the correct dataset name
                if not (adata.obs['dataset'] == dataset_name).all():
                    logger.warning(f"Setting all samples in {dataset_name} to have consistent dataset label")
                    adata.obs['dataset'] = dataset_name
            
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

def create_combined_sparse_matrix(datasets: Dict[str, ad.AnnData], gene_to_idx: Dict[str, int]) -> Tuple[sp.csr_matrix, List[str], List[str], pd.DataFrame, np.ndarray, Dict[str, List[str]]]:
    """Create combined sparse matrix of gene expression data using COO format."""
    logger.info("Creating combined sparse matrix using COO intermediate")

    total_genes = len(gene_to_idx)
    total_samples = sum(adata.n_obs for adata in datasets.values())
    logger.info(f"Combined matrix dimensions: {total_samples} samples x {total_genes} genes")

    # --- Prepare for COO construction ---
    rows = []
    cols = []
    data = []
    # --- End COO preparation ---

    # Track sample metadata
    sample_ids = []
    dataset_labels = []
    obs_data = { # Use lists for efficiency
        'sample_id': [], 'subject_id': [], 'dataset': [], 'data_type': [],
        'expression_unit': [], 'tissue': [], 'sex': [], 'age': []
    }
    gene_presence = np.zeros((total_samples, total_genes), dtype=bool) # Keep this for now
    gene_datasets = {gene_id: [] for gene_id in gene_to_idx}

    sample_offset = 0
    processed_samples_count = 0

    for dataset_name, adata in datasets.items():
        n_samples_ds = adata.n_obs
        n_genes_ds = adata.n_vars
        logger.info(f"Processing {dataset_name} dataset with {n_samples_ds} samples and {n_genes_ds} genes")

        # --- Efficiently create gene mapping for this dataset ---
        local_gene_idx_to_global_col_idx = {}
        valid_local_gene_indices = [] # Store local indices that map globally

        # Determine the gene identifier to use (prioritize ensembl_id)
        has_ensembl_id_col = 'ensembl_id' in adata.var.columns

        for local_idx, gene_var_name in enumerate(adata.var_names):
            target_gene_id = None
            if has_ensembl_id_col:
                # Use 'ensembl_id' column if available and valid in gene_to_idx
                ensembl_id = adata.var.iloc[local_idx]['ensembl_id']
                if pd.notna(ensembl_id) and ensembl_id != '' and ensembl_id in gene_to_idx:
                    target_gene_id = ensembl_id
                # Fallback to var_name if it's a standard ID and ensembl_id is missing/invalid
                elif str(gene_var_name).startswith(('ENSG', 'gSpikein', 'ENTREZ:')):
                     if str(gene_var_name) in gene_to_idx:
                        target_gene_id = str(gene_var_name)
            # If no ensembl_id column, use var_name if it's a standard ID
            elif str(gene_var_name).startswith(('ENSG', 'gSpikein', 'ENTREZ:')):
                 if str(gene_var_name) in gene_to_idx:
                    target_gene_id = str(gene_var_name)

            if target_gene_id and target_gene_id in gene_to_idx:
                global_col_idx = gene_to_idx[target_gene_id]
                local_gene_idx_to_global_col_idx[local_idx] = global_col_idx
                gene_datasets[target_gene_id].append(dataset_name)
                valid_local_gene_indices.append(local_idx)
        logger.info(f"  Mapped {len(local_gene_idx_to_global_col_idx)} local genes to global index for {dataset_name}")
        # --- End gene mapping ---

        # --- Process sample metadata ---
        current_sample_ids = adata.obs.index.astype(str).tolist()
        sample_ids.extend(current_sample_ids)
        dataset_labels.extend([dataset_name] * n_samples_ds)

        obs_data['sample_id'].extend(current_sample_ids)
        obs_data['dataset'].extend([dataset_name] * n_samples_ds)
        for col in ['subject_id', 'data_type', 'expression_unit', 'tissue', 'sex', 'age']:
            if col in adata.obs.columns:
                obs_data[col].extend(adata.obs[col].tolist())
            else:
                obs_data[col].extend([''] * n_samples_ds)
        # --- End sample metadata ---

        # --- Efficiently extract and map non-zero expression data ---
        logger.info(f"  Extracting non-zero values for {dataset_name}...")
        # Ensure adata.X is in COO format for efficient iteration
        if isinstance(adata.X, np.ndarray): # Handle dense input
            logger.warning(f"  Input matrix for {dataset_name} is dense. Converting to COO for processing.")
            source_matrix_coo = sp.coo_matrix(adata.X)
        elif not isinstance(adata.X, sp.coo_matrix):
            source_matrix_coo = adata.X.tocoo()
        else:
            source_matrix_coo = adata.X

        # Get data pointers - this avoids repeated lookups
        source_rows = source_matrix_coo.row
        source_cols = source_matrix_coo.col
        source_data = source_matrix_coo.data

        # Iterate through non-zero elements of the source matrix
        num_vals = len(source_data)
        count_added = 0
        for i in range(num_vals):
            local_sample_idx = source_rows[i]
            local_gene_idx = source_cols[i]
            value = source_data[i]

            # Check if this gene is one we want to keep and map
            if local_gene_idx in local_gene_idx_to_global_col_idx:
                global_col_idx = local_gene_idx_to_global_col_idx[local_gene_idx]
                global_row_idx = local_sample_idx + sample_offset

                # Append to COO lists
                rows.append(global_row_idx)
                cols.append(global_col_idx)
                data.append(value)
                count_added += 1

                # Update gene presence matrix (if still needed, consider alternatives)
                gene_presence[global_row_idx, global_col_idx] = True

            # Log progress occasionally within large datasets
            if num_vals > 1000000 and i > 0 and i % (num_vals // 20) == 0 : # Log more frequently
                 progress = (i / num_vals) * 100
                 logger.info(f"    {dataset_name} data extraction: {progress:.0f}% complete ({count_added} values added)")

        logger.info(f"  Finished extracting {count_added} values for {dataset_name}.")
        # --- End data extraction ---

        processed_samples_count += n_samples_ds
        logger.info(f"  Processed {processed_samples_count}/{total_samples} total samples.")
        sample_offset += n_samples_ds
        # --- End dataset loop ---

    # Create obs DataFrame
    logger.info("Creating final obs DataFrame...")
    obs_df = pd.DataFrame(obs_data)
    # Use sample_id as the index
    if len(obs_df) == len(sample_ids):
        obs_df.index = sample_ids
    else:
        logger.error(f"Length mismatch between obs_data ({len(obs_df)}) and sample_ids ({len(sample_ids)}). Cannot set index correctly.")
        # Fallback or raise error
        obs_df = obs_df.set_index('sample_id', drop=False) # Keep sample_id col

    # --- Build the final sparse matrix ---
    logger.info(f"Constructing final sparse matrix from {len(data)} non-zero values...")
    combined_sparse_coo = sp.coo_matrix((data, (rows, cols)), shape=(total_samples, total_genes), dtype=np.float32) # Specify dtype

    # Convert to CSR for efficiency (duplicates are summed automatically by CSR constructor)
    logger.info("Converting to CSR format (summing duplicates)...")
    combined_sparse_csr = combined_sparse_coo.tocsr()
    logger.info("CSR conversion complete.")
    # --- End matrix build ---

    # Calculate sparsity
    non_zero = combined_sparse_csr.nnz
    total_elements = combined_sparse_csr.shape[0] * combined_sparse_csr.shape[1]
    sparsity = 1.0 - (non_zero / total_elements) if total_elements > 0 else 0.0
    logger.info(f"Combined matrix sparsity: {sparsity:.4f} ({non_zero} non-zero elements)")

    # Memory usage estimate (less precise than before but gives an idea)
    csr_memory_mb = (combined_sparse_csr.data.nbytes + combined_sparse_csr.indices.nbytes + combined_sparse_csr.indptr.nbytes) / (1024*1024)
    dense_memory_mb = (total_samples * total_genes * np.dtype(np.float32).itemsize) / (1024*1024)
    logger.info(f"Memory usage: CSR ~{csr_memory_mb:.2f} MB, Dense ~{dense_memory_mb:.2f} MB")

    # Ensure gene_datasets contains unique dataset names per gene
    for gene_id in gene_datasets:
         gene_datasets[gene_id] = sorted(list(set(gene_datasets[gene_id])))


    return combined_sparse_csr, sample_ids, dataset_labels, obs_df, gene_presence, gene_datasets


def create_var_dataframe(gene_ids: List[str], gene_datasets: Dict[str, List[str]], reference_mapping: pd.DataFrame) -> pd.DataFrame:
    """Create var DataFrame with gene metadata using optimized mapping."""
    logger.info("Creating var DataFrame with gene metadata (optimized)")

    # Create var DataFrame indexed by gene_ids
    var_df = pd.DataFrame(index=gene_ids)
    # Important: Set index name to avoid potential conflicts later
    var_df.index.name = 'feature_id'
    var_df['gene_id'] = gene_ids # Keep the gene_id column as well

    # --- Optimized Mapping ---
    logger.info("Preparing reference mapping for optimized lookup...")
    # Prepare reference mapping for faster lookup using map, ensure index is string
    reference_mapping['gene_id'] = reference_mapping['gene_id'].astype(str)
    ref_map = reference_mapping.set_index('gene_id')
    logger.info(f"Reference map created with {len(ref_map)} entries.")

    # Map basic info using .map()
    logger.info("Mapping basic gene info (chromosome, name, type)...")
    var_df['chromosome'] = var_df.index.map(ref_map.get('chromosome', pd.Series(dtype='str'))).fillna('')
    var_df['gene_name'] = var_df.index.map(ref_map.get('gene_name', pd.Series(dtype='str'))).fillna('')
    var_df['gene_type'] = var_df.index.map(ref_map.get('gene_type', pd.Series(dtype='str'))).fillna('')

    # map 'mapping_source' based on whether the index was found in ref_map
    logger.info("Determining mapping source...")
    var_df['mapping_source'] = np.where(var_df.index.isin(ref_map.index), 'reference_mapping', 'unmapped')

    # Keep original_ids as the index (which is our target gene_id)
    # If a different original ID is needed, it should be fetched from ref_map
    var_df['original_ids'] = var_df['gene_id']
    # Example: fetch original ID *if* it exists in reference mapping
    # var_df['original_from_ref'] = var_df.index.map(ref_map.get('original_gene_id', pd.Series(dtype='str'))).fillna('')


    # Map dataset presence (loop is still reasonable here due to set logic)
    logger.info("Mapping dataset presence...")
    present_col = []
    for gene_id in var_df.index:
        datasets = gene_datasets.get(str(gene_id), []) # Ensure lookup key is string
        datasets = [ds for ds in datasets if ds != 'combined'] # Ensure no combined
        present_col.append(','.join(sorted(set(datasets))))
    var_df['present_in_datasets'] = present_col
    # --- End Optimized Mapping ---

    # Handle spike-ins specifically if needed (example)
    logger.info("Handling spike-in metadata...")
    spike_mask = var_df.index.str.startswith('gSpikein', na=False) # Added na=False
    if spike_mask.any():
         num_spikes = spike_mask.sum()
         logger.info(f"Assigning metadata for {num_spikes} spike-ins")
         var_df.loc[spike_mask, 'gene_name'] = var_df.index[spike_mask]
         var_df.loc[spike_mask, 'gene_type'] = 'spike_in_control'
         var_df.loc[spike_mask, 'chromosome'] = 'spike_in'
         var_df.loc[spike_mask, 'mapping_source'] = 'spike_in'

    # --- Final Cleanup ---
    # Reset index if needed, or ensure it's named correctly
    # var_df = var_df.reset_index().rename(columns={'index': 'gene_id'}) # Option 1
    # Option 2: Keep index as feature_id, ensure gene_id column is correct
    var_df['gene_id'] = var_df.index.astype(str) # Ensure gene_id col matches index

    # Ensure correct order of columns (optional but good practice)
    desired_cols = ['gene_id', 'present_in_datasets', 'chromosome', 'mapping_source', 'gene_name', 'original_ids', 'gene_type']
    # Add any other columns that might exist
    all_cols = desired_cols + [col for col in var_df.columns if col not in desired_cols]
    var_df = var_df[all_cols]

    logger.info(f"Created var DataFrame with {len(var_df)} genes.")
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