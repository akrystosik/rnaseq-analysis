# File List with Contents

The following files were found, along with their contents:

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/anndata_save_wrapper.py`

````
#!/usr/bin/env python3
"""
Wrapper script to save AnnData objects with proper string conversion.
Usage: python anndata_save_wrapper.py input_file output_file
"""

import sys
import os
import scanpy as sc
import pandas as pd
import numpy as np
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('anndata_save_wrapper')

def convert_to_serializable(obj):
    """Convert dict values to serializable types."""
    if isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif pd.isna(obj):
        return ""
    elif isinstance(obj, (pd.Series, pd.DataFrame, np.ndarray)):
        if hasattr(obj, "tolist"):
            return obj.tolist()
        else:
            return str(obj)
    else:
        return str(obj)

def save_adata_safely(input_file, output_file):
    """Load AnnData and save it with proper string conversion."""
    try:
        logger.info(f"Loading AnnData from {input_file}")
        adata = sc.read_h5ad(input_file)
        
        logger.info(f"Loaded AnnData with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Make a copy of the original uns
        original_uns = adata.uns.copy()
        
        # Convert all uns values to serializable types
        logger.info("Converting uns dictionary to serializable types")
        adata.uns = convert_to_serializable(original_uns)
        
        # Save the AnnData object
        logger.info(f"Saving AnnData to {output_file}")
        adata.write(output_file)
        
        # Verify the save worked
        logger.info("Verifying saved file")
        test_adata = sc.read_h5ad(output_file)
        logger.info(f"Verification successful: {test_adata.n_obs} samples, {test_adata.n_vars} genes")
        
        return True
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python anndata_save_wrapper.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    success = save_adata_safely(input_file, output_file)
    
    if success:
        print(f"Successfully saved AnnData to {output_file}")
        sys.exit(0)
    else:
        print("Failed to save AnnData")
        sys.exit(1)

````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/create_combined_dataset_all_genes_sparse.py`

````
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
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/entrez-to-ensembl-mapping.py`

````
#!/usr/bin/env python3
"""
Entrez to Ensembl Mapping Creator

This script creates a mapping between Entrez Gene IDs and Ensembl IDs
using NCBI's gene2ensembl file, which is the definitive source for this mapping.

Usage:
  python entrez-to-ensembl-mapping.py --output /path/to/output.csv [--species human]
"""

import os
import sys
import argparse
import pandas as pd
import urllib.request
import gzip
import io
import logging
import time
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('entrez_to_ensembl')

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
DEFAULT_OUTPUT = BASE_DIR / "metadata/json/entrez_to_ensembl_mapping.csv"

def get_taxon_id(species):
    """Get NCBI taxonomy ID for a species."""
    species_tax_map = {
        'human': '9606',
        'mouse': '10090',
        'rat': '10116',
        'fly': '7227',
        'worm': '6239'
    }
    return species_tax_map.get(species.lower())

def create_gene_mapping(output_file=DEFAULT_OUTPUT, species="human"):
    """
    Create a mapping file between Entrez Gene IDs and Ensembl IDs.
    
    Parameters:
    -----------
    output_file : str
        Path to save the mapping file
    species : str
        Species to generate mapping for ('human', 'mouse', etc.)
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Get taxon ID for species
        taxon_id = get_taxon_id(species)
        if not taxon_id:
            logger.error(f"Invalid species: {species}")
            return False
            
        logger.info(f"Creating mapping for {species} (taxon ID {taxon_id})")
        
        # Download NCBI gene2ensembl file directly to memory
        logger.info("Downloading gene2ensembl from NCBI...")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
        
        response = urllib.request.urlopen(url)
        compressed_file = io.BytesIO(response.read())
        
        # Prepare to store the mapping data
        mapping_data = []
        total_lines = 0
        species_count = 0
        
        # Parse the file
        logger.info("Parsing gene2ensembl file...")
        with gzip.open(compressed_file, 'rt') as f:
            # Skip header
            next(f)
            
            # Process each line
            for line in f:
                total_lines += 1
                fields = line.strip().split('\t')
                
                # Check if this is the species we want
                if fields[0] == taxon_id:
                    entrez_id = fields[1]
                    ensembl_gene = fields[2]
                    
                    # Store the mapping
                    mapping_data.append({
                        'entrez_id': entrez_id,
                        'ensembl_id': ensembl_gene
                    })
                    
                    species_count += 1
                
                # Show progress
                if total_lines % 1000000 == 0:
                    logger.info(f"Processed {total_lines:,} lines, found {species_count:,} {species} gene mappings")
        
        # Create a DataFrame directly from the mapping data
        mapping_df = pd.DataFrame(mapping_data)
        mapping_df.to_csv(output_file, index=False)
        
        logger.info(f"Created mapping file with {len(mapping_df):,} entries")
        logger.info(f"Saved mapping to {output_file}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error creating gene mapping: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Entrez to Ensembl Mapping Creator')
    parser.add_argument('--output', default=DEFAULT_OUTPUT, help='Path to save the mapping file')
    parser.add_argument('--species', default='human', choices=['human', 'mouse', 'rat', 'fly', 'worm'],
                       help='Species to generate mapping for')
    parser.add_argument('--force', action='store_true', help='Force regeneration even if output file exists')
    
    args = parser.parse_args()
    
    # Check if output already exists and we're not forcing regeneration
    if os.path.exists(args.output) and not args.force:
        logger.info(f"Output file {args.output} already exists. Use --force to regenerate.")
        return
    
    start_time = time.time()
    success = create_gene_mapping(args.output, args.species)
    
    if success:
        logger.info(f"Successfully created mapping in {time.time() - start_time:.2f} seconds")
    else:
        logger.error("Failed to create mapping")
        sys.exit(1)

if __name__ == '__main__':
    main()
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/fix_placeholder_ids.py`

````
#!/usr/bin/env python3
"""
Fix placeholder IDs in preprocessed datasets by converting them to proper Entrez IDs.
"""
import os
import sys
import re
import logging
import scanpy as sc
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('placeholder_id_fixer')

def fix_placeholder_ids(h5ad_file, output_file):
    """
    Fix placeholder IDs in the specified h5ad file.
    
    Parameters:
    -----------
    h5ad_file : str
        Path to the h5ad file to fix
    output_file : str
        Path to save the fixed h5ad file
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Load the dataset
        logger.info(f"Loading dataset from {h5ad_file}")
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Check if there are placeholder IDs
        placeholder_pattern = re.compile(r'^PLACEHOLDER_(.+)$')
        
        # Find placeholder IDs in gene_id and ensembl_id columns
        placeholders = []
        for col in ['gene_id', 'ensembl_id']:
            if col in adata.var.columns:
                # Count placeholders
                placeholder_count = sum(1 for g in adata.var[col] if isinstance(g, str) and g.startswith('PLACEHOLDER_'))
                if placeholder_count > 0:
                    logger.info(f"Found {placeholder_count} placeholder IDs in {col} column")
                    placeholders.append(col)
        
        if not placeholders:
            logger.info("No placeholder IDs found, no fix needed")
            return True
        
        # Make a copy of var DataFrame
        var_df = adata.var.copy()
        
        # Fix categorical columns to prevent comparison issues
        for col in var_df.columns:
            if pd.api.types.is_categorical_dtype(var_df[col]):
                logger.info(f"Converting categorical column {col} to string")
                var_df[col] = var_df[col].astype(str)
        
        # Replace placeholder IDs with Entrez IDs
        entrez_count = 0
        other_count = 0
        
        for col in placeholders:
            for idx in var_df.index:
                value = var_df.loc[idx, col]
                if isinstance(value, str) and value.startswith('PLACEHOLDER_'):
                    # Extract the ID from the placeholder
                    match = placeholder_pattern.match(value)
                    if match:
                        id_part = match.group(1)
                        
                        # If it's a numeric ID, use Entrez prefix
                        if id_part.isdigit():
                            var_df.loc[idx, col] = f"ENTREZ:{id_part}"
                            entrez_count += 1
                        else:
                            # For non-numeric IDs, keep as is but remove PLACEHOLDER_
                            var_df.loc[idx, col] = id_part
                            other_count += 1
        
        # Replace the var DataFrame
        adata.var = var_df
        
        # Log results
        logger.info(f"Fix results: {entrez_count + other_count} placeholder IDs converted")
        logger.info(f"Total Entrez IDs: {entrez_count}")
        logger.info(f"Total other IDs: {other_count}")
        
        # Save the fixed dataset
        logger.info(f"Saving fixed dataset to {output_file}")
        adata.write_h5ad(output_file)
        
        logger.info("Fix completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error fixing placeholder IDs: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fix_placeholder_ids.py <input_h5ad_file> [output_h5ad_file]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        # Use same file with .fixed suffix
        output_file = input_file + ".fixed"
    
    success = fix_placeholder_ids(input_file, output_file)
    
    if success:
        print(f"Successfully fixed placeholder IDs and saved to {output_file}")
        sys.exit(0)
    else:
        print(f"Failed to fix placeholder IDs")
        sys.exit(1)

````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/gene_id_mapping_reference.py`

````
#!/usr/bin/env python3
"""
Gene ID Reference Mapping Generator

This script creates a comprehensive gene ID reference mapping database
that serves as the single source of truth for all gene identifiers.
It combines GENCODE v24 annotations with ENCODE/ENTEx numeric ID mappings
to create a unified mapping between different gene ID formats.

Usage:
    python gene_id_mapping_reference.py \
        --encode-dir /path/to/encode/data \
        --entex-dir /path/to/entex/data \
        --gencode-gtf /path/to/gencode.v24.gtf \
        --entrez-mapping /path/to/entrez_to_ensembl_mapping.csv \
        --output /path/to/gene_id_reference_mapping.csv
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import gzip
import re
from pathlib import Path
import glob
import urllib.request
import shutil
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gene_id_mapping_generator')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate comprehensive gene ID reference mapping')
    parser.add_argument('--encode-dir', type=str, required=True,
                        help='Directory containing ENCODE data')
    parser.add_argument('--entex-dir', type=str, required=True,
                        help='Directory containing ENTEx data')
    parser.add_argument('--gencode-gtf', type=str,
                        help='Path to GENCODE v24 GTF file (will download if not provided)')
    parser.add_argument('--entrez-mapping', type=str, required=True,
                        help='Path to Entrez to Ensembl mapping CSV')
    parser.add_argument('--output', type=str, required=True,
                        help='Output file for gene ID reference mapping')
    parser.add_argument('--temp-dir', type=str, default='/tmp',
                        help='Temporary directory for downloaded files')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration of mapping even if output exists')
    
    return parser.parse_args()

def download_gencode_gtf(output_path, temp_dir='/tmp'):
    """Download GENCODE v24 GTF file if not already present."""
    logger.info("Downloading GENCODE v24 GTF file...")
    
    # GENCODE v24 URL (hg38)
    gencode_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz"
    
    # Create temp directory if it doesn't exist
    os.makedirs(temp_dir, exist_ok=True)
    
    # Download gzipped file
    temp_file = os.path.join(temp_dir, "gencode.v24.annotation.gtf.gz")
    try:
        urllib.request.urlretrieve(gencode_url, temp_file)
        
        # Decompress and save to output path
        with gzip.open(temp_file, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        logger.info(f"GENCODE v24 GTF file downloaded and saved to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Error downloading GENCODE v24 GTF file: {e}")
        return False

def parse_gencode_gtf(gtf_file):
    """Parse GENCODE GTF file to extract gene information."""
    logger.info(f"Parsing GENCODE GTF file: {gtf_file}")
    
    genes = {}
    gene_count = 0
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                
                # Extract gene information
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                
                # Parse attributes
                attr_dict = {}
                attr_string = fields[8]
                for attr in attr_string.split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    
                    key_value = attr.split(' ', 1)
                    if len(key_value) != 2:
                        continue
                    
                    key, value = key_value
                    value = value.strip('"')
                    attr_dict[key] = value
                
                # Required attributes for a gene
                if 'gene_id' not in attr_dict:
                    continue
                
                original_gene_id = attr_dict['gene_id']  # Keep full ID with version
                gene_id = original_gene_id.split('.')[0]  # Base ID for mapping                
                gene_name = attr_dict.get('gene_name', '')
                gene_type = attr_dict.get('gene_type', '')

                genes[gene_id] = {
                    'gene_id': gene_id,
                    'original_gene_id': original_gene_id,  # Add this field
                    'gene_name': gene_name,                
                    'gene_type': gene_type,
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                }
                
                gene_count += 1
                if gene_count % 10000 == 0:
                    logger.info(f"Processed {gene_count} genes...")
        
        logger.info(f"Parsed {gene_count} genes from GENCODE GTF file")
        return genes
    except Exception as e:
        logger.error(f"Error parsing GENCODE GTF file: {e}")
        return {}

def load_encode_numeric_ids(encode_dir):
    """Load numeric IDs from ENCODE dataset."""
    logger.info(f"Loading numeric IDs from ENCODE dataset: {encode_dir}")
    
    numeric_ids = set()
    
    # Look for TSV files in ENCODE directory
    tsv_files = glob.glob(os.path.join(encode_dir, '**/*.tsv'), recursive=True)
    
    if not tsv_files:
        logger.warning(f"No TSV files found in {encode_dir}")
        return numeric_ids
    
    # Process a sample of files to identify gene IDs
    for file_path in tsv_files[:10]:  # Process first 10 files
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#')
            if 'gene_id' in df.columns:
                # Check if these are numeric IDs
                sample_ids = df['gene_id'].iloc[:5].astype(str).tolist()
                if any(id.isdigit() for id in sample_ids):
                    numeric_ids.update(df['gene_id'].astype(str))
                    logger.info(f"Found numeric IDs in {file_path}")
        except Exception as e:
            logger.warning(f"Error reading {file_path}: {e}")
    
    logger.info(f"Found {len(numeric_ids)} numeric IDs in ENCODE dataset")
    return numeric_ids

def load_entex_gene_ids(entex_dir):
    """Load gene IDs from ENTEx dataset."""
    logger.info(f"Loading gene IDs from ENTEx dataset: {entex_dir}")
    
    gene_ids = defaultdict(set)  # Dictionary to store different ID types
    
    # Look for TSV files in ENTEx directory
    tsv_files = glob.glob(os.path.join(entex_dir, '**/*.tsv'), recursive=True)
    
    if not tsv_files:
        logger.warning(f"No TSV files found in {entex_dir}")
        return gene_ids
    
    # Process a sample of files to identify gene IDs
    for file_path in tsv_files[:10]:  # Process first 10 files
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#')
            if 'gene_id' in df.columns:
                # Categorize IDs by type
                for id_str in df['gene_id'].astype(str):
                    if id_str.startswith('ENSG'):
                        gene_ids['ensembl'].add(id_str.split('.')[0])  # Remove version
                    elif id_str.isdigit():
                        gene_ids['numeric'].add(id_str)
                    elif id_str.startswith('gSpikein'):
                        gene_ids['spike_in'].add(id_str)
                    else:
                        gene_ids['other'].add(id_str)
        except Exception as e:
            logger.warning(f"Error reading {file_path}: {e}")
    
    # Log counts by ID type
    for id_type, ids in gene_ids.items():
        logger.info(f"Found {len(ids)} {id_type} IDs in ENTEx dataset")
    
    return gene_ids

def load_entrez_to_ensembl_mapping(mapping_file):
    """Load Entrez to Ensembl ID mapping."""
    logger.info(f"Loading Entrez to Ensembl mapping from {mapping_file}")
    
    try:
        # Load mapping file
        mapping_df = pd.read_csv(mapping_file)
        logger.info(f"Loaded mapping with {len(mapping_df)} entries")
        
        # Create mapping dictionary
        entrez_to_ensembl = {}
        for _, row in mapping_df.iterrows():
            entrez_id = str(row['entrez_id'])
            ensembl_id = row['ensembl_id']
            
            entrez_to_ensembl[entrez_id] = ensembl_id
        
        logger.info(f"Created Entrez to Ensembl mapping with {len(entrez_to_ensembl)} entries")
        return entrez_to_ensembl
    
    except Exception as e:
        logger.error(f"Error loading Entrez to Ensembl mapping: {e}")
        return {}

def load_standardized_datasets(data_dir):
    """Load standardized datasets to extract gene ID information."""
    logger.info(f"Loading standardized datasets from {data_dir}")
    
    datasets = {}
    
    # Find standardized h5ad files
    h5ad_files = glob.glob(os.path.join(data_dir, '*_standardized_v*.h5ad'))
    
    for file_path in h5ad_files:
        dataset_name = os.path.basename(file_path).split('_')[0]
        logger.info(f"Loading {dataset_name} dataset from {file_path}")
        
        try:
            adata = sc.read_h5ad(file_path)
            datasets[dataset_name] = adata
            
            # Log dataset information
            logger.info(f"Loaded {dataset_name} with {adata.n_obs} samples and {adata.n_vars} genes")
            
            # Check if var contains gene_id column
            if 'gene_id' in adata.var.columns:
                logger.info(f"Found gene_id column in {dataset_name} var DataFrame")
            
            # Check var_names format
            var_names_sample = list(adata.var_names[:5])
            logger.info(f"Sample var_names in {dataset_name}: {var_names_sample}")
            
        except Exception as e:
            logger.error(f"Error loading {file_path}: {e}")
    
    return datasets

def create_gene_id_mapping(gencode_genes, encode_numeric_ids, entex_gene_ids, entrez_to_ensembl, standardized_datasets):
    """Create comprehensive gene ID mapping."""
    logger.info("Creating comprehensive gene ID mapping...")
    
    # Initialize mapping DataFrame with GENCODE genes
    mapping_data = []
    for gene_id, gene_info in gencode_genes.items():
        mapping_data.append({
            'gene_id': gene_id,
            'gene_name': gene_info['gene_name'],
            'gene_type': gene_info['gene_type'],
            'chromosome': gene_info['chromosome'],
            'start': gene_info['start'],
            'end': gene_info['end'],
            'strand': gene_info['strand'],
            'numeric_id': None,
            'source': 'GENCODE',
            'mapping_confidence': 'high'
        })
    
    # Add Entrez to Ensembl mappings
    entrez_added = set()
    
    for entrez_id, ensembl_id in entrez_to_ensembl.items():
        # Skip if Ensembl ID is not in GENCODE (which would be unusual)
        if ensembl_id not in gencode_genes:
            continue
        
        # Find the GENCODE entry for this Ensembl ID
        for item in mapping_data:
            if item['gene_id'] == ensembl_id:
                item['numeric_id'] = entrez_id
                entrez_added.add(entrez_id)
                break
    
    logger.info(f"Added {len(entrez_added)} Entrez to Ensembl mappings")
    
    # Add remaining Entrez IDs (ones not in GENCODE)
    missing_entrez = encode_numeric_ids - entrez_added
    logger.info(f"Found {len(missing_entrez)} Entrez IDs not in GENCODE mapping")
    
    for entrez_id in missing_entrez:
        # Add as placeholder
        mapping_data.append({
            'gene_id': f"PLACEHOLDER_{entrez_id}",
            'gene_name': f"Unknown_{entrez_id}",
            'gene_type': 'unknown',
            'chromosome': '',
            'start': 0,
            'end': 0,
            'strand': '',
            'numeric_id': entrez_id,
            'source': 'ENCODE',
            'mapping_confidence': 'low'
        })
    
    # Handle ENTEx gene IDs
    # 1. Add Ensembl IDs that aren't in GENCODE
    for ensembl_id in entex_gene_ids['ensembl']:
        if ensembl_id not in gencode_genes:
            mapping_data.append({
                'gene_id': ensembl_id,
                'gene_name': '',
                'gene_type': 'unknown',
                'chromosome': '',
                'start': 0,
                'end': 0,
                'strand': '',
                'numeric_id': None,
                'source': 'ENTEx',
                'mapping_confidence': 'medium'
            })
    
    # 2. Add spike-in controls with special handling
    for spike_in_id in entex_gene_ids['spike_in']:
        mapping_data.append({
            'gene_id': spike_in_id,
            'gene_name': spike_in_id,
            'gene_type': 'spike_in_control',
            'chromosome': 'spike_in',
            'start': 0,
            'end': 0,
            'strand': '',
            'numeric_id': None,
            'source': 'ENTEx',
            'mapping_confidence': 'high'
        })
    
    # Create DataFrame from mapping data
    mapping_df = pd.DataFrame(mapping_data)
    logger.info(f"Created mapping with {len(mapping_df)} entries")
    
    # Log some statistics
    logger.info(f"Mapped numeric IDs: {sum(mapping_df['numeric_id'].notna())}")
    logger.info(f"Mapping confidence counts: {mapping_df['mapping_confidence'].value_counts().to_dict()}")
    
    return mapping_df

def main():
    """Main function to create gene ID reference mapping."""
    args = parse_arguments()
    
    # Check if output already exists and we're not forcing regeneration
    if os.path.exists(args.output) and not args.force:
        logger.info(f"Output file {args.output} already exists. Use --force to regenerate.")
        return
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Step 1: Get GENCODE GTF file
    gencode_gtf_file = args.gencode_gtf
    if not gencode_gtf_file or not os.path.exists(gencode_gtf_file):
        temp_gtf_file = os.path.join(args.temp_dir, "gencode.v24.annotation.gtf")
        if not download_gencode_gtf(temp_gtf_file, args.temp_dir):
            logger.error("Failed to download GENCODE GTF file. Please provide it manually.")
            return
        gencode_gtf_file = temp_gtf_file
    
    # Step 2: Parse GENCODE GTF to extract gene information
    gencode_genes = parse_gencode_gtf(gencode_gtf_file)
    if not gencode_genes:
        logger.error("Failed to parse GENCODE GTF file.")
        return
    
    # Step 3: Load NCBI Entrez to Ensembl mapping
    entrez_to_ensembl = load_entrez_to_ensembl_mapping(args.entrez_mapping)
    if not entrez_to_ensembl:
        logger.error("Failed to load Entrez to Ensembl mapping.")
        return
    
    # Step 4: Load numeric IDs from ENCODE dataset
    encode_numeric_ids = load_encode_numeric_ids(args.encode_dir)
    
    # Step 5: Load gene IDs from ENTEx dataset
    entex_gene_ids = load_entex_gene_ids(args.entex_dir)
    
    # Step 6: Load standardized datasets for additional information
    standardized_dir = os.path.dirname(args.output)
    standardized_datasets = load_standardized_datasets(standardized_dir)
    
    # Step 7: Create comprehensive gene ID mapping
    mapping_df = create_gene_id_mapping(
        gencode_genes, encode_numeric_ids, entex_gene_ids, entrez_to_ensembl, standardized_datasets
    )
    
    # Step 8: Save mapping to output file
    mapping_df.to_csv(args.output, index=False)
    logger.info(f"Gene ID reference mapping saved to {args.output}")
    
    # Step 9: Create a JSON version for easier loading
    json_output = args.output.replace('.csv', '.json')
    mapping_df.to_json(json_output, orient='records')
    logger.info(f"Gene ID reference mapping also saved as JSON to {json_output}")

if __name__ == "__main__":
    main()
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/generate_encode_mapping.py`

````
#!/usr/bin/env python3
"""
ENCODE Gene ID Mapping Generator

This script processes ENCODE gene IDs from raw TSV files and creates a mapping
between original IDs and standard Ensembl IDs without version numbers.
"""

import os
import glob
import pandas as pd
import numpy as np
import logging
import argparse
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('encode_id_mapper')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate ENCODE gene ID mapping')
    parser.add_argument('--encode-dir', type=str, required=True,
                        help='Directory containing raw ENCODE data')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Directory to save mapping files')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration of mapping even if exists')
    
    return parser.parse_args()

def collect_gene_ids_from_file(file_path):
    """Extract gene IDs from a single TSV file."""
    try:
        df = pd.read_csv(file_path, sep='\t')
        
        if 'gene_id' not in df.columns:
            logger.warning(f"No gene_id column in {file_path}")
            return []
        
        # Get unique gene IDs
        gene_ids = df['gene_id'].unique().tolist()
        logger.info(f"Collected {len(gene_ids)} unique gene IDs from {os.path.basename(file_path)}")
        return gene_ids
    
    except Exception as e:
        logger.error(f"Error processing {file_path}: {e}")
        return []

def collect_all_gene_ids(encode_dir):
    """Collect gene IDs from all ENCODE TSV files."""
    # Find all TSV files
    tsv_files = []
    for root, dirs, files in os.walk(encode_dir):
        for file in files:
            if file.endswith('.tsv'):
                tsv_files.append(os.path.join(root, file))
    
    logger.info(f"Found {len(tsv_files)} TSV files in ENCODE directory")
    
    # Process each file and collect unique gene IDs
    all_gene_ids = set()
    for file in tsv_files:
        gene_ids = collect_gene_ids_from_file(file)
        all_gene_ids.update(gene_ids)
    
    # Convert to list
    gene_id_list = list(all_gene_ids)
    logger.info(f"Collected {len(gene_id_list)} unique gene IDs from all files")
    
    return gene_id_list

def extract_ensembl_id(id_str):
    """Extract Ensembl gene ID from various formats."""
    id_str = str(id_str)
    
    # Direct Ensembl ID
    if id_str.startswith('ENSG'):
        return id_str.split('.')[0]  # Remove version number
    
    # Pipe-delimited composite entry
    if '|' in id_str:
        fields = id_str.split('|')
        if len(fields) > 1:
            # Look for ENSG in the second field
            ensembl_field = fields[1]
            if ensembl_field.startswith('ENSG'):
                return ensembl_field.split('.')[0]  # Remove version number
    
    # Numeric ID (potential Entrez)
    if id_str.isdigit():
        return f"ENTREZ:{id_str}"
    
    # Special case: spike-in controls
    if id_str.startswith('gSpikein'):
        return id_str
    
    # Other format - can't extract Ensembl ID
    return None

def process_gene_ids(gene_ids):
    """Process collected gene IDs and extract Ensembl IDs."""
    # Create DataFrame with gene IDs
    gene_ids_df = pd.DataFrame({'gene_id': gene_ids})
    
    # Apply the extraction function
    gene_ids_df['extracted_ensembl_id'] = gene_ids_df['gene_id'].apply(extract_ensembl_id)
    
    # Calculate statistics
    total_ids = len(gene_ids_df)
    extracted_count = gene_ids_df['extracted_ensembl_id'].notna().sum()
    extraction_rate = extracted_count / total_ids * 100
    
    direct_ensembl = sum(1 for id in gene_ids_df['gene_id'] if str(id).startswith('ENSG'))
    composite_with_ensembl = sum(1 for i, row in gene_ids_df.iterrows() 
                                if not str(row['gene_id']).startswith('ENSG') and 
                                row['extracted_ensembl_id'] is not None and
                                str(row['extracted_ensembl_id']).startswith('ENSG'))
    entrez_ids = sum(1 for id in gene_ids_df['extracted_ensembl_id'] if id is not None and str(id).startswith('ENTREZ:'))
    spikein_ids = sum(1 for id in gene_ids_df['extracted_ensembl_id'] if id is not None and str(id).startswith('gSpikein'))
    unmapped_ids = total_ids - extracted_count
    
    # Report statistics
    logger.info(f"Extraction statistics:")
    logger.info(f"  - Total IDs processed: {total_ids}")
    logger.info(f"  - Successfully extracted: {extracted_count} ({extraction_rate:.2f}%)")
    logger.info(f"  - Direct Ensembl IDs: {direct_ensembl} ({direct_ensembl/total_ids*100:.2f}%)")
    logger.info(f"  - Composite entries with Ensembl: {composite_with_ensembl} ({composite_with_ensembl/total_ids*100:.2f}%)")
    logger.info(f"  - Entrez IDs: {entrez_ids} ({entrez_ids/total_ids*100:.2f}%)")
    logger.info(f"  - Spike-in controls: {spikein_ids} ({spikein_ids/total_ids*100:.2f}%)")
    logger.info(f"  - Unmapped IDs: {unmapped_ids} ({unmapped_ids/total_ids*100:.2f}%)")
    
    return gene_ids_df

def create_mapping_file(processed_df, output_dir):
    """Create mapping files from processed gene IDs."""
    # Save processed IDs with details
    processed_file = os.path.join(output_dir, 'encode_processed_gene_ids.csv')
    processed_df.to_csv(processed_file, index=False)
    logger.info(f"Saved processed gene IDs to {processed_file}")
    
    # Create a mapping file for integration with pipeline
    mapping_df = processed_df[processed_df['extracted_ensembl_id'].notna()].copy()
    mapping_df = mapping_df[['gene_id', 'extracted_ensembl_id']].rename(
        columns={'gene_id': 'original_id', 'extracted_ensembl_id': 'ensembl_id'}
    )
    
    # Save mapping file
    mapping_file = os.path.join(output_dir, 'encode_id_to_ensembl_mapping.csv')
    mapping_df.to_csv(mapping_file, index=False)
    logger.info(f"Saved ID mapping to {mapping_file}")
    
    return mapping_file

def main():
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    mapping_file = os.path.join(args.output_dir, 'encode_id_to_ensembl_mapping.csv')
    if os.path.exists(mapping_file) and not args.force:
        logger.info(f"Mapping file {mapping_file} already exists. Use --force to regenerate.")
        return mapping_file
    
    # Collect gene IDs from ENCODE files
    gene_ids = collect_all_gene_ids(args.encode_dir)
    
    # Process gene IDs
    processed_df = process_gene_ids(gene_ids)
    
    # Create mapping file
    mapping_file = create_mapping_file(processed_df, args.output_dir)
    
    return mapping_file

if __name__ == "__main__":
    main()
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/preprocess_dataset_gene_ids.py`

````
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

# In scripts/pipeline/preprocess_dataset_gene_ids.py

def preprocess_encode_dataset(adata: ad.AnnData, numeric_to_ensembl: Dict[str, str], ensembl_to_info: Dict[str, Dict], output_file: str) -> None:
    """Preprocess ENCODE dataset to add gene metadata and consistent IDs."""
    logger.info(f"Preprocessing ENCODE dataset with {adata.n_vars} genes")
    
    # --- Load Mappings ---
    # Load the enhanced ENCODE ID mapping (original_id -> target_id)
    encode_mapping_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_mapping/encode_id_to_ensembl_mapping.csv'
    try:
        encode_mapping_df = pd.read_csv(encode_mapping_file)
        # Create a simple dict for lookup: original_id -> target_id (ensembl or ENTREZ:...)
        encode_orig_to_target = dict(zip(encode_mapping_df['original_id'].astype(str), encode_mapping_df['ensembl_id'].astype(str)))
        logger.info(f"Loaded ENCODE ID mapping with {len(encode_orig_to_target)} entries")
    except Exception as e:
        logger.error(f"Error loading ENCODE ID mapping: {e}")
        encode_orig_to_target = {}

    # Enhanced Entrez map (used as fallback if ENCODE mapping gave Entrez ID)
    # This numeric_to_ensembl comes from the reference mapping (gene_id_mapping_reference.py)
    # Let's rename it for clarity within this function
    ref_numeric_to_ensembl = numeric_to_ensembl 
    logger.info(f"Using reference numeric->ensembl map with {len(ref_numeric_to_ensembl)} entries")
    
    # --- Prepare New Var Columns ---
    var_columns = {
        'original_id_from_input': [], # Store the var_name from the input AnnData (numeric string)
        'original_id_from_source': [],# Store the complex ID string from V1's 'original_ids' column
        'ensembl_id': [],          # Final Standard Ensembl ID (or ENTREZ: if unmapped)
        'gene_name': [],           # Gene symbol
        'gene_type': [],           # Gene type (protein_coding, etc.)
        'chromosome': [],          # Chromosome
        'mapping_source': [],      # Source of mapping
        'mapping_confidence': []   # Confidence of mapping
    }
    
    mapped_count = 0
    unmapped_count = 0
    # --- Iterate through input adata's var ---
    # Use .iterrows() to access both index (numeric string) and row data easily
    for idx, row in adata.var.iterrows():
        input_var_name = str(idx) # The numeric string index, e.g., '10904'
        original_id_source = str(row.get('original_ids', '')) # The complex ID from stage 1, e.g., '10904' or 'ENST...|ENSG...|...'
        
        var_columns['original_id_from_input'].append(input_var_name)
        var_columns['original_id_from_source'].append(original_id_source)

        target_id = '' # This will hold the result of the mapping
        mapping_source = 'unmapped' # Default
        
        # --- Mapping Logic ---
        # 1. Try mapping using the complex original_id_source via encode_orig_to_target
        if original_id_source in encode_orig_to_target:
            target_id = encode_orig_to_target[original_id_source]
            mapping_source = 'encode_mapping'
            
        # 2. If encode mapping resulted in an ENTREZ ID, try to map it further using reference map
        if target_id.startswith('ENTREZ:'):
            entrez_numeric = target_id[7:]
            if entrez_numeric in ref_numeric_to_ensembl:
                potential_ensembl = ref_numeric_to_ensembl[entrez_numeric]
                # Check if the reference map gave a valid Ensembl ID (not placeholder)
                if potential_ensembl.startswith('ENSG'):
                    target_id = potential_ensembl # Update target_id to Ensembl
                    mapping_source = 'encode_mapping -> ref_numeric'
                else:
                    # Keep the ENTREZ: ID, mapping source remains encode_mapping
                    pass 
            # else: Keep the ENTREZ: ID from encode_mapping

        # 3. If still no valid target_id AND the input var name is numeric, try mapping it via reference numeric map
        elif not target_id.startswith('ENSG') and input_var_name.isdigit():
             if input_var_name in ref_numeric_to_ensembl:
                  potential_ensembl = ref_numeric_to_ensembl[input_var_name]
                  if potential_ensembl.startswith('ENSG'):
                       target_id = potential_ensembl
                       mapping_source = 'ref_numeric_fallback'
                  else: # Reference map gave placeholder or something else
                       target_id = f"ENTREZ:{input_var_name}" # Fallback to ENTREZ ID
                       mapping_source = 'entrez_id_fallback'
             else:
                  # Numeric input ID not in reference map
                  target_id = f"ENTREZ:{input_var_name}" # Fallback to ENTREZ ID
                  mapping_source = 'entrez_id_fallback'

        # 4. Handle spike-ins explicitly if they appear in original ID
        elif original_id_source.startswith('gSpikein'):
              target_id = original_id_source
              mapping_source = 'spike_in'

        # 5. If target_id is still empty or not standard, mark as unmapped
        if not target_id.startswith(('ENSG', 'ENTREZ:', 'gSpikein')):
             target_id = input_var_name # Use the input numeric string as fallback ID
             mapping_source = 'unmapped'
             unmapped_count += 1
        else:
             mapped_count += 1

        var_columns['ensembl_id'].append(target_id) # Store the final ID (Ensembl, Entrez, Spikein, or fallback)
        var_columns['mapping_source'].append(mapping_source)
        
        # --- Add Gene Info based on final target_id ---
        final_lookup_id = target_id.split('.')[0] # Use base ID for info lookup
        if final_lookup_id in ensembl_to_info:
            info = ensembl_to_info[final_lookup_id]
            var_columns['gene_name'].append(info['gene_name'])
            var_columns['gene_type'].append(info['gene_type'])
            var_columns['chromosome'].append(info['chromosome'])
            var_columns['mapping_confidence'].append(info['mapping_confidence'])
        elif target_id.startswith('gSpikein'):
            var_columns['gene_name'].append(target_id)
            var_columns['gene_type'].append('spike_in_control')
            var_columns['chromosome'].append('spike_in')
            var_columns['mapping_confidence'].append('high')
        else: # Entrez ID or Unmapped
            var_columns['gene_name'].append('')
            var_columns['gene_type'].append('unknown')
            var_columns['chromosome'].append('')
            var_columns['mapping_confidence'].append('low' if target_id.startswith('ENTREZ:') else 'none')

    # Create new var DataFrame
    new_var = pd.DataFrame(var_columns)
    
    # --- Set the final index ---
    # Use the 'ensembl_id' column (which now contains ENSG, ENTREZ:, gSpikein, or fallback) as the new index
    # Handle potential duplicates in the target IDs by adding a suffix using ad.utils.make_index_unique
    new_var.index = ad.utils.make_index_unique(pd.Index(new_var['ensembl_id'].astype(str))) # <<< CORRECTED LINE
    new_var.index.name = 'feature_id' # Set index name

    # Add a final 'gene_id' column that matches the new index
    new_var['gene_id'] = new_var.index.astype(str) # Ensure the column is also string
        
    # Calculate mapping statistics
    mapping_percentage = mapped_count / adata.n_vars * 100 if adata.n_vars > 0 else 0
    
    logger.info(f"ENCODE mapping details (Stage 2.5):")
    logger.info(f"  - Mapped source counts: {pd.Series(var_columns['mapping_source']).value_counts().to_dict()}")
    logger.info(f"  Total mapped: {mapped_count}/{adata.n_vars} genes ({mapping_percentage:.2f}%)")
    
    # Replace var DataFrame
    adata.var = new_var
    
    # Store additional metadata in uns
    adata.uns['gene_mapping_stats'] = {
        'mapped_genes': mapped_count,
        'unmapped_count': unmapped_count,
        'total_genes': adata.n_vars,
        'mapping_percentage': mapping_percentage,
        'mapping_sources': pd.Series(var_columns['mapping_source']).value_counts().to_dict()
    }
    
    # Save preprocessed dataset
    logger.info(f"Saving preprocessed ENCODE dataset to {output_file}")
    
    # Use the safe save wrapper script
    import subprocess
    save_script = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline/anndata_save_wrapper.py'
    temp_save_path = output_file + ".tmp.h5ad"
    try:
        adata.write_h5ad(temp_save_path) # Write temporarily first
        logger.info(f"Temporarily wrote AnnData to {temp_save_path}")
        
        # Call the wrapper to perform the final safe save
        result = subprocess.run(
            ['python', save_script, temp_save_path, output_file],
            capture_output=True, text=True, check=True
        )
        logger.info("anndata_save_wrapper script output:")
        logger.info(result.stdout)
        if result.stderr:
            logger.error("anndata_save_wrapper script error output:")
            logger.error(result.stderr)
        
        logger.info(f"Successfully saved preprocessed ENCODE dataset via wrapper to {output_file}")
        os.remove(temp_save_path) # Clean up temp file

    except subprocess.CalledProcessError as e:
         logger.error(f"Error running anndata_save_wrapper script: {e}")
         logger.error(f"Stderr: {e.stderr}")
         logger.error(f"Stdout: {e.stdout}")
         # Attempt direct write as fallback, though it might fail
         try:
              logger.warning(f"Attempting direct write to {output_file} as fallback.")
              adata.write(output_file)
         except Exception as direct_write_e:
              logger.error(f"Direct write fallback also failed: {direct_write_e}")
    except Exception as e:
        logger.error(f"Error during final save process for ENCODE: {e}")
        if os.path.exists(temp_save_path):
             os.remove(temp_save_path)
    
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
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/rnaseq_utils.py`

````
#!/usr/bin/env python3
"""
RNA-seq Data Processing Utilities

This module contains utility functions for RNA-seq data processing,
standardization, and validation used across the RNA-seq standardization pipeline.
"""

import os
import json
import logging
import pandas as pd
import numpy as np
from pathlib import Path

# Set up logging
logger = logging.getLogger('rnaseq_utils')

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
GENCODE_MAPPING_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gencode_v24_complete_mapping.csv")
METADATA_JSON_DIR = BASE_DIR / "metadata/json"

# Constants
CORE_METADATA_FIELDS = [
    'sample_id',       # Unique identifier for each sample
    'subject_id',      # Identifier for the subject/donor
    'sex',             # Biological sex (standardized as 'male', 'female', 'unknown')
    'age',             # Age at sample collection
    'tissue',          # Tissue or brain region
    'dataset',         # Source dataset identifier
    'data_type',       # Type of data (RNA-seq, microarray, etc.)
    'expression_unit'  # Unit of measurement (TPM, FPKM, counts, etc.)
]

# ===================================================================================
# File Loading Functions
# ===================================================================================

def load_json_mapping(filename, default=None):
    """
    Load a JSON mapping file from the metadata directory.
    
    Args:
        filename: Name of the JSON file
        default: Default value to return if file not found
        
    Returns:
        Dictionary with mapping data
    """
    filepath = METADATA_JSON_DIR / filename
    
    try:
        if filepath.exists():
            with open(filepath, 'r') as f:
                return json.load(f)
        else:
            logger.warning(f"Mapping file not found: {filepath}")
            return default if default is not None else {}
    except Exception as e:
        logger.error(f"Error loading mapping file {filepath}: {e}")
        return default if default is not None else {}

def load_mappings():
    """
    Load all standard mapping files into a single dictionary.
    
    Returns:
        Dictionary containing all mappings
    """
    mappings = {
        'TISSUE_TO_UBERON': load_json_mapping('tissue_to_uberon.json', {}),
        'ASSAY_TO_EFO': load_json_mapping('assay_to_efo.json', {}),
        'AGE_TO_HSAPDV': load_json_mapping('age_to_hsapdv.json', {}),
        'SEX_STANDARDIZATION': load_json_mapping('sex_standardization.json', {}),
        'SPECIES_TO_NCBI_TAXON': load_json_mapping('species_to_taxon.json', {})
    }
    
    # Load dataset-specific metadata
    for dataset in ['encode', 'gtex', 'mage', 'adni', 'entex']:
        key = f"{dataset.upper()}_METADATA"
        mappings[key] = load_json_mapping(f"{dataset}_metadata.json", {})
    
    return mappings

def load_gencode_mapping():
    """
    Load GENCODE mapping from CSV file.
    
    Returns:
        Dictionary mapping Ensembl IDs to gene information
    """
    try:
        if not os.path.exists(GENCODE_MAPPING_FILE):
            logger.warning(f"GENCODE mapping file not found: {GENCODE_MAPPING_FILE}")
            return {}
        
        logger.info(f"Loading GENCODE mapping from {GENCODE_MAPPING_FILE}")
        gencode_df = pd.read_csv(GENCODE_MAPPING_FILE)
        
        # Create dictionary for faster lookups
        gencode_dict = {}
        
        for _, row in gencode_df.iterrows():
            if 'gene_id' in row:
                # Get the gene ID (with or without version)
                gene_id = row['gene_id']
                
                # Also extract the base ID without version
                base_id = standardize_ensembl_id(gene_id)
                
                info = {
                    'gene_name': row['gene_name'] if 'gene_name' in row else '',
                    'gene_type': row['gene_type'] if 'gene_type' in row else '',
                    'chromosome': row['chromosome'] if 'chromosome' in row else ''
                }
                
                # Map both the full ID and base ID
                gencode_dict[gene_id] = info
                gencode_dict[base_id] = info
        
        logger.info(f"Loaded GENCODE mapping for {len(gencode_dict)} gene IDs")
        return gencode_dict
    
    except Exception as e:
        logger.warning(f"Error loading GENCODE mapping: {e}")
        return {}

# ===================================================================================
# Gene ID Processing Functions
# ===================================================================================

def standardize_ensembl_id(gene_id):
    """
    Standardize Ensembl gene ID by removing version number.
    
    Args:
        gene_id: Gene ID to standardize
        
    Returns:
        Standardized gene ID
    """
    if pd.isna(gene_id):
        return None
    
    # Convert to string if needed
    gene_id = str(gene_id)
    
    # Remove version if present
    if '.' in gene_id and gene_id.startswith('ENSG'):
        return gene_id.split('.')[0]
    return gene_id

def add_gencode_annotations(var_df, gencode_mapping):
    """
    Add GENCODE annotations to variable DataFrame.
    
    Args:
        var_df: DataFrame with gene information
        gencode_mapping: Dictionary mapping gene IDs to gene information
        
    Returns:
        DataFrame with added GENCODE annotations
    """
    # Pre-allocate arrays for efficiency
    gene_count = len(var_df)
    gene_names = np.empty(gene_count, dtype=object)
    gene_types = np.empty(gene_count, dtype=object)
    chromosomes = np.empty(gene_count, dtype=object)
    mapping_sources = np.empty(gene_count, dtype=object)
    
    for i, (index, _) in enumerate(var_df.iterrows()):
        ensembl_id = index
        
        # Try exact match first
        if ensembl_id in gencode_mapping:
            info = gencode_mapping[ensembl_id]
            gene_names[i] = info['gene_name']
            gene_types[i] = info['gene_type']
            chromosomes[i] = info['chromosome']
            mapping_sources[i] = 'exact_match'
        else:
            # Try base ID match (without version)
            base_id = standardize_ensembl_id(ensembl_id)
            if base_id in gencode_mapping:
                info = gencode_mapping[base_id]
                gene_names[i] = info['gene_name']
                gene_types[i] = info['gene_type']
                chromosomes[i] = info['chromosome']
                mapping_sources[i] = 'base_id_match'
            else:
                # No mapping found
                gene_names[i] = ""
                gene_types[i] = ""
                chromosomes[i] = ""
                mapping_sources[i] = 'unmapped'
    
    # Add annotations to DataFrame efficiently
    var_df['gene_name'] = gene_names
    var_df['gene_type'] = gene_types
    var_df['chromosome'] = chromosomes
    var_df['mapping_source'] = mapping_sources
    
    return var_df

# ===================================================================================
# Metadata Processing Functions
# ===================================================================================

def map_age_to_hsapdv(age_value, age_mapping=None):
    """
    Map an age value to the appropriate HsapDv ontology term.
    
    Args:
        age_value: Age value as a string
        age_mapping: Dictionary mapping age to HsapDv terms
        
    Returns:
        HsapDv ontology ID or empty string if mapping not found
    """
    if not age_mapping:
        age_mapping = load_json_mapping('age_to_hsapdv.json', {})
        
    if not age_value or age_value == '' or pd.isna(age_value):
        return ""
    
    # Check if age is already a developmental stage
    if age_value.lower() in age_mapping:
        return age_mapping[age_value.lower()]
    
    # Try to convert to integer age
    try:
        age_years = int(age_value)
        
        # Map to age ranges
        if age_years < 2:
            return age_mapping.get("0-1", "")  # infant
        elif age_years < 13:
            return age_mapping.get("2-12", "")  # child
        elif age_years < 20:
            return age_mapping.get("13-19", "")  # adolescent
        elif age_years < 40:
            return age_mapping.get("20-39", "")  # young adult
        elif age_years < 60:
            return age_mapping.get("40-59", "")  # middle aged adult
        else:
            return age_mapping.get("60+", "")  # elderly
    except ValueError:
        # Not a simple integer - check for ranges or other formats
        # This could be expanded with more pattern matching
        return ""


def map_tissue_to_ontology(tissue_name, mappings):
    """
    Map a tissue name to an ontology ID.
    
    Args:
        tissue_name: Name of the tissue
        mappings: Dictionary containing mapping information
        
    Returns:
        Tuple of (ontology_id, confidence)
    """
    if not tissue_name or pd.isna(tissue_name):
        return "", "none"
    
    # Normalize tissue name
    tissue_lower = str(tissue_name).lower().strip()
    
    # Get tissue mappings from the mappings dictionary
    tissue_to_uberon = mappings.get('TISSUE_TO_UBERON', {})
    
    # Check for direct match
    if tissue_lower in tissue_to_uberon:
        return tissue_to_uberon[tissue_lower], "high"
    
    # Check for case-insensitive match
    for known_tissue, ontology_id in tissue_to_uberon.items():
        if known_tissue.lower() == tissue_lower:
            return ontology_id, "high"
    
    # Try substring matching for approximate matches
    for known_tissue, ontology_id in tissue_to_uberon.items():
        if known_tissue.lower() in tissue_lower or tissue_lower in known_tissue.lower():
            return ontology_id, "medium"
    
    # No match found
    return "", "none"

def standardize_metadata(metadata_df, dataset_name, mappings=None):
    """
    Standardize metadata across datasets to ensure consistent fields and values.
    Maps fields to standard ontology terms where possible.
    
    Args:
        metadata_df: Metadata DataFrame
        dataset_name: Source dataset name
        mappings: Dictionary with mapping information
        
    Returns:
        Standardized metadata DataFrame with ontology mappings
    """
    # Load mappings if not provided
    if mappings is None:
        mappings = load_mappings()
    
    # Make a copy to avoid modifying the original
    df = metadata_df.copy()
    
    # Get specific mappings
    tissue_to_uberon = mappings.get('TISSUE_TO_UBERON', {})
    assay_to_efo = mappings.get('ASSAY_TO_EFO', {})
    sex_standardization = mappings.get('SEX_STANDARDIZATION', {})
    species_to_taxon = mappings.get('SPECIES_TO_NCBI_TAXON', {})
    
    # Add dataset identifier
    df['dataset'] = dataset_name
    
    # Ensure all core fields exist
    for field in CORE_METADATA_FIELDS:
        if field not in df.columns:
            df[field] = ""
    
    # Add species field if not present - default to human
    if 'species' not in df.columns:
        df['species'] = "human"
        df['species_ontology'] = species_to_taxon.get("human", "")
    else:
        # Map species to NCBI Taxon
        df['species_ontology'] = df['species'].map(
            lambda x: species_to_taxon.get(x, species_to_taxon.get("unknown", ""))
        )
    
    # Standardize sex values
    if 'sex' in df.columns:
        df['sex'] = df['sex'].map(lambda x: sex_standardization.get(x, "unknown"))
    
    # Standardize tissue to Uberon ontology
    if 'tissue' in df.columns:
        # Keep original tissue names
        df['tissue_original'] = df['tissue']
        # Normalize tissue names (lowercase, strip whitespace)
        normalized_tissue = df['tissue'].astype(str).str.lower().str.strip()
        
        # Map to standardized tissue names
        df['tissue'] = normalized_tissue.map(
            lambda x: x if pd.isna(x) else x
        )
        
        # Map to Uberon IDs
        df['tissue_ontology'] = normalized_tissue.map(
            lambda x: tissue_to_uberon.get(x, "")
        )
        
        # For missing mappings, try a more flexible approach
        mask = df['tissue_ontology'] == ""
        if mask.any():
            for tissue_name, uberon_id in tissue_to_uberon.items():
                # Find tissues that contain this tissue name as a substring
                substr_mask = normalized_tissue.str.contains(tissue_name.lower(), na=False)
                # Update only those that don't already have a mapping
                update_mask = mask & substr_mask
                if update_mask.any():
                    df.loc[update_mask, 'tissue_ontology'] = uberon_id
    
    # Standardize data_type to EFO ontology
    if 'data_type' in df.columns:
        # Keep original data type
        df['data_type_original'] = df['data_type']
        
        # Map to EFO IDs
        df['assay_ontology'] = df['data_type'].map(
            lambda x: assay_to_efo.get(x, "")
        )
        
        # Also try to map extraction_method if available
        if 'extraction_method' in df.columns:
            # For rows without assay_ontology, try mapping from extraction_method
            mask = df['assay_ontology'] == ""
            df.loc[mask, 'assay_ontology'] = df.loc[mask, 'extraction_method'].map(
                lambda x: assay_to_efo.get(x, "")
            )

    # Handle missing subject_id by using donor_id as fallback
    if 'subject_id' not in df.columns and 'donor_id' in df.columns:
        logger.info(f"Using donor_id as subject_id for {dataset_name}")
        df['subject_id'] = df['donor_id']
    elif 'subject_id' not in df.columns and 'donor' in df.columns:
        logger.info(f"Using donor as subject_id for {dataset_name}")
        df['subject_id'] = df['donor']    
    
    
    # For unmapped tissues, add them to a pending file for curator review
    if 'tissue' in df.columns:
        unmapped_tissues = []
        for tissue in df['tissue'].dropna().unique():
            tissue_ontology, confidence = map_tissue_to_ontology(tissue, mappings)
            if not tissue_ontology and tissue not in unmapped_tissues:
                unmapped_tissues.append(tissue)
        

        if unmapped_tissues:
            logger.warning(f"{dataset_name}: Unmapped tissues: {', '.join(unmapped_tissues)}")
            # Only log the unmapped tissues, don't try to save to file as mappings_dir might not be defined
            logger.info(f"Unmapped tissues for potential mapping: {', '.join(unmapped_tissues)}")        
            # Add to pending file for curator review
            # Just log the unmapped tissues instead of trying to write to a file
            logger.info(f"Unmapped tissues for {dataset_name}: {', '.join([str(t) for t in unmapped_tissues])}")
    

    
    # Standardize age - try to map age ranges to developmental stages
    if 'age' in df.columns:
        # Keep original age values
        df['age_original'] = df['age']
        
        # Convert age to string and clean
        df['age'] = df['age'].astype(str).replace('nan', '').replace('None', '')
        
        # Map age ranges to developmental stages
        df['developmental_stage_ontology'] = df['age'].apply(
            lambda x: map_age_to_hsapdv(x, mappings.get('AGE_TO_HSAPDV', {}))
        )
    
    # Convert categorical fields to categorical data type
    categorical_fields = ['sex', 'tissue', 'dataset', 'cell_type', 'disease', 'ethnicity', 'data_type']
    for field in categorical_fields:
        if field in df.columns:
            df[field] = pd.Categorical(df[field].astype(str))
    
    return df

# ===================================================================================
# AnnData Processing Functions
# ===================================================================================

def prepare_for_anndata(df, var_df, data_df):
    """
    Prepare dataframes for AnnData creation with consistent data types.
    
    Args:
        df: Observation metadata DataFrame (samples)
        var_df: Variable metadata DataFrame (genes)
        data_df: Expression data DataFrame
        
    Returns:
        Processed observation, variable, and data DataFrames
    """
    # Make copies to avoid modifying the original DataFrames
    df = df.copy()
    var_df = var_df.copy()
    
    # Handle NA values in observation metadata
    for col in df.columns:
        if pd.api.types.is_numeric_dtype(df[col]):
            # For numeric columns, fill NaNs with 0
            df[col] = df[col].fillna(0)
        elif isinstance(df[col].dtype, pd.CategoricalDtype):
            # For categorical columns, convert to string first
            # Then back to categorical with empty string as a category
            str_series = df[col].astype(str).fillna("")
            categories = list(df[col].cat.categories)
            if "" not in categories:
                categories.append("")
            df[col] = pd.Categorical(str_series, categories=categories)
        else:
            # For other columns (object/string), fill NaNs with empty string
            df[col] = df[col].fillna("")
    
    # Handle NA values in variable metadata
    for col in var_df.columns:
        if col != 'gene_id':  # Don't convert index/gene_id
            if isinstance(var_df[col].dtype, pd.CategoricalDtype):
                # Same approach as above for categorical columns
                str_series = var_df[col].astype(str).fillna("")
                categories = list(var_df[col].cat.categories)
                if "" not in categories:
                    categories.append("")
                var_df[col] = pd.Categorical(str_series, categories=categories)
            elif pd.api.types.is_numeric_dtype(var_df[col]):
                var_df[col] = var_df[col].fillna(0)
            else:
                var_df[col] = var_df[col].fillna("")
    
    # Ensure data matrix has no NaN values
    data_df = data_df.fillna(0)
    
    return df, var_df, data_df

# ===================================================================================
# Validation Functions
# ===================================================================================

def validate_metadata(metadata_df, dataset_name=None):
    """
    Validate standardized metadata and log warnings for missing or invalid values.
    
    Args:
        metadata_df: Standardized metadata DataFrame
        dataset_name: Dataset name for logging purposes
        
    Returns:
        (validated_df, validation_report) tuple
    """
    # Make a copy to avoid modifying the original
    df = metadata_df.copy()
    
    # Initialize validation report
    report = {
        'dataset': dataset_name or 'unknown',
        'sample_count': len(df),
        'missing_fields': [],
        'missing_values': {},
        'unmapped_tissues': [],
        'unmapped_assays': [],
        'unmapped_ages': [],
        'validation_status': 'passed'
    }
    
    # Check required fields
    required_fields = ['sample_id', 'subject_id', 'sex', 'tissue', 'dataset', 'data_type']
    missing_fields = [field for field in required_fields if field not in df.columns]
    
    if missing_fields:
        logger.warning(f"{dataset_name}: Missing required fields in metadata: {', '.join(missing_fields)}")
        report['missing_fields'] = missing_fields
        report['validation_status'] = 'warning'
    
    # Check for missing values in key fields
    for field in df.columns:
        missing_count = df[field].isna().sum()
        if missing_count > 0:
            logger.warning(f"{dataset_name}: Field '{field}' has {missing_count} missing values out of {len(df)} samples")
            report['missing_values'][field] = int(missing_count)
            report['validation_status'] = 'warning'
    
    # Check ontology mappings
    if 'tissue_ontology' in df.columns:
        missing_ontology = (df['tissue_ontology'] == "").sum()
        if missing_ontology > 0:
            logger.warning(f"{dataset_name}: Missing tissue ontology mappings for {missing_ontology} samples out of {len(df)}")
            # Log the unmapped tissues
            unmapped_tissues = df.loc[df['tissue_ontology'] == "", 'tissue'].unique()
            unmapped_tissues = [t for t in unmapped_tissues if t and not pd.isna(t)]
            logger.warning(f"{dataset_name}: Unmapped tissues: {', '.join(unmapped_tissues)}")
            report['unmapped_tissues'] = list(unmapped_tissues)
            report['validation_status'] = 'warning'
    
    if 'assay_ontology' in df.columns:
        missing_ontology = (df['assay_ontology'] == "").sum()
        if missing_ontology > 0:
            logger.warning(f"{dataset_name}: Missing assay ontology mappings for {missing_ontology} samples out of {len(df)}")
            # Log the unmapped assays
            unmapped_assays = df.loc[df['assay_ontology'] == "", 'data_type'].unique()
            unmapped_assays = [a for a in unmapped_assays if a and not pd.isna(a)]
            logger.warning(f"{dataset_name}: Unmapped assays: {', '.join(unmapped_assays)}")
            report['unmapped_assays'] = list(unmapped_assays)
            report['validation_status'] = 'warning'
    
    if 'developmental_stage_ontology' in df.columns:
        missing_ontology = (df['developmental_stage_ontology'] == "").sum()
        if missing_ontology > 0 and 'age' in df.columns:
            logger.warning(f"{dataset_name}: Missing developmental stage mappings for {missing_ontology} samples out of {len(df)}")
            # Log the unmapped ages
            unmapped_ages = df.loc[df['developmental_stage_ontology'] == "", 'age'].unique()
            unmapped_ages = [a for a in unmapped_ages if a and not pd.isna(a)]
            logger.warning(f"{dataset_name}: Unmapped ages: {', '.join(unmapped_ages)}")
            report['unmapped_ages'] = list(unmapped_ages)
            report['validation_status'] = 'warning'
    
    # Add validation report to the DataFrame as metadata
    df.attrs['validation_report'] = report
    
    return df, report

def validate_dataset_metadata(metadata, dataset_name):
    """
    Validate that dataset metadata meets requirements.
    
    Args:
        metadata: Dictionary with dataset metadata
        dataset_name: Name of the dataset
        
    Returns:
        True if valid, raises ValueError otherwise
    """
    # Check required fields
    required_fields = ['harmonized_gencode_version', 'harmonized_reference_genome']
    for field in required_fields:
        if field not in metadata:
            raise ValueError(f"Dataset {dataset_name} metadata is missing required field '{field}'")
    
    # Validate GENCODE version
    harmonized_gencode = metadata['harmonized_gencode_version']
    if harmonized_gencode not in ["24", "v24"]:
        raise ValueError(f"Dataset {dataset_name} must have harmonized_gencode_version set to 'v24', got '{harmonized_gencode}'")
    
    # Validate reference genome
    harmonized_genome = metadata['harmonized_reference_genome']
    if harmonized_genome not in ["hg38", "GRCh38"]:
        raise ValueError(f"Dataset {dataset_name} must have harmonized_reference_genome set to 'hg38' or 'GRCh38', got '{harmonized_genome}'")
    
    return True

def load_dataset_specific_metadata(metadata_dir, dataset_name):
    """
    Load dataset-specific metadata from JSON file.
    
    Args:
        metadata_dir: Directory containing dataset metadata JSON files
        dataset_name: Name of the dataset to load metadata for
        
    Returns:
        Dictionary with dataset metadata or None if not found
    """
    if not metadata_dir:
        return None
    
    # Try both lowercase and original case filenames
    potential_filenames = [
        f"{dataset_name.lower()}_metadata.json",
        f"{dataset_name}_metadata.json"
    ]
    
    for filename in potential_filenames:
        metadata_path = os.path.join(metadata_dir, filename)
        if os.path.exists(metadata_path):
            try:
                logger.info(f"Loading dataset-specific metadata from {metadata_path}")
                with open(metadata_path, 'r') as f:
                    metadata = json.load(f)
                
                # Validate metadata meets requirements
                validate_dataset_metadata(metadata, dataset_name)
                
                return metadata
            except ValueError as ve:
                # Re-raise validation errors
                raise ve
            except Exception as e:
                logger.warning(f"Error loading dataset metadata from {metadata_path}: {e}")
                return None
    
    logger.info(f"No dataset-specific metadata found for {dataset_name}")
    return None


def apply_dataset_specific_metadata(adata, metadata):
    """
    Apply dataset-specific metadata from JSON to AnnData object.

    Args:
        adata: AnnData object to update
        metadata: Dictionary with metadata to apply

    Returns:
        Updated AnnData object
    """
    if not metadata:
        return adata

    logger.info("Applying dataset-specific metadata from JSON file")

    # Apply top-level attributes to uns
    for key, value in metadata.items():
        if key not in ['obs_columns', 'metadata_source']:  # Handle these separately
            # Store directly in uns - ensure_serializable will handle it later
            adata.uns[key] = value
            # Note: Complex nested structures are handled by ensure_serializable before saving

    # Apply observation column updates if present
    if 'obs_columns' in metadata:
        for col_name, value in metadata['obs_columns'].items():
            logger.info(f"  Updating obs column {col_name}")
            # Check if value should be applied to all rows or is a map
            if isinstance(value, dict):
                 logger.warning(f"Applying dict values to obs column '{col_name}' is not directly supported here. Check metadata structure.")
                 # Potentially map using index if value is a dict keyed by sample ID
            elif isinstance(value, list) and len(value) == adata.n_obs:
                 adata.obs[col_name] = value # Assign list if length matches
            elif not isinstance(value, (dict, list)):
                 adata.obs[col_name] = value # Assign scalar value to all rows
            else:
                 logger.warning(f"Could not apply value of type {type(value)} to obs column '{col_name}'")


    # Add metadata source information
    if 'metadata_source' in metadata:
        if 'metadata_sources' not in adata.uns:
            adata.uns['metadata_sources'] = []
        elif not isinstance(adata.uns['metadata_sources'], list):
            # If it exists but isn't a list, make it a list containing the original item
            adata.uns['metadata_sources'] = [str(adata.uns['metadata_sources'])] # Convert existing item to string

        # --- Modification: ALWAYS Append as JSON string ---
        source_dict = metadata['metadata_source']
        try:
            # Serialize the source_dict itself first (using the util from standardize_datasets)
            # This requires importing ensure_serializable or duplicating it here
            # Simpler: directly dump to JSON, assuming source_dict is simple enough
            source_str = json.dumps(source_dict, default=str) # Use default=str for safety
            adata.uns['metadata_sources'].append(source_str)
            logger.debug(f"Appended metadata_source as JSON string: {source_str[:100]}...")
        except Exception as e:
            logger.error(f"Could not serialize metadata_source to JSON: {e}. Appending as raw string.")
            adata.uns['metadata_sources'].append(str(source_dict))

    logger.info("Dataset-specific metadata applied successfully")
    return adata

````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/run_rnaseq_pipeline.sh`

````
#!/bin/bash
# Run the complete RNA-seq standardization pipeline with improved gene ID mapping
# Modified to add date suffix to output directory for version control

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR="${BASE_DIR}/scripts/pipeline"
METADATA_DIR="${BASE_DIR}/metadata/json"

# Create a timestamp for output directory and log file
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="${BASE_DIR}/standardized_data/run_${TIMESTAMP}"
PREPROCESSED_DIR="${BASE_DIR}/preprocessed_data/run_${TIMESTAMP}"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/pipeline_${TIMESTAMP}.log"

# Create directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$PREPROCESSED_DIR"
mkdir -p "$METADATA_DIR"

# Function to log messages
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# Function to run a command with logging

# Function to safely save AnnData files
save_anndata_safely() {
    local input_file="$1"
    local output_file="$2"
    
    if [ -f "$input_file" ]; then
        log_message "Using save wrapper to safely save AnnData: $output_file"
        python "${SCRIPTS_DIR}/anndata_save_wrapper.py" "$input_file" "$output_file"
        return $?
    else
        log_message "Error: Input file not found: $input_file"
        return 1
    fi
}

run_command() {
    log_message "Running: $1"
    eval "$1" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

# Log the start of the pipeline with versioned directories
log_message "=== Starting RNA-seq Pipeline with Versioned Output ==="
log_message "Output directory: $OUTPUT_DIR"
log_message "Preprocessed directory: $PREPROCESSED_DIR"
log_message "Log file: $LOG_FILE"

# Step 0: Generate Entrez to Ensembl mapping
log_message "=== Stage 0: Generating Entrez to Ensembl Mapping ==="
ENTREZ_MAPPING="${METADATA_DIR}/entrez_to_ensembl_mapping.csv"

# Check if --force flag was passed
FORCE_FLAG=""
if [ "$1" == "--force" ] || [ "$1" == "--force-mapping" ]; then
    FORCE_FLAG="--force"
    log_message "Force flag detected. Will regenerate mapping files."
fi

run_command "python ${SCRIPTS_DIR}/entrez-to-ensembl-mapping.py --output ${ENTREZ_MAPPING} --species human ${FORCE_FLAG}"

# Check if mapping generation was successful
if [ $? -eq 0 ]; then
    log_message "Entrez to Ensembl mapping generated or verified successfully!"
else
    log_message "Warning: Entrez to Ensembl mapping generation failed. Check the log file for details."
    log_message "Will continue without this mapping, which may affect ENCODE gene ID mapping quality."
fi


# Create temp directory for intermediate files
mkdir -p "${OUTPUT_DIR}/temp"
# Step 1: Initial Data Conversion (existing)
log_message "=== Stage 1: Initial Data Conversion ==="
run_command "python ${SCRIPTS_DIR}/standardize_datasets.py \\
    --encode-dir \"${BASE_DIR}/encode/raw_data/cell_lines\" \\
    --encode-entex-dir \"${BASE_DIR}/encode/entex\" \\
    --entex-metadata-file \"${BASE_DIR}/encode/metadata/entex_metadata.json\" \\
    --mage-dir \"${BASE_DIR}/mage\" \\
    --adni-dir \"${BASE_DIR}/adni_microarray\" \\
    --gtex-file \"${BASE_DIR}/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz\" \\
    --metadata-dir \"${METADATA_DIR}\" \\
    --output-dir \"$OUTPUT_DIR\""

if [ $? -ne 0 ]; then
    log_message "Stage 1 (Initial Data Conversion) failed. Check the log file for details. Exiting."
    exit 1
fi

# # Check if Stage 1 ran successfully
# if [ $? -eq 0 ]; then
#     log_message "Stage 1 completed successfully!"
# else
#     log_message "Stage 1 failed. Check the log file for details."
#     exit 1
# fi


# # Use save wrapper to properly save AnnData files
# if [ -f "${OUTPUT_DIR}/temp/encode_standardized_v1.h5ad" ]; then
#     save_anndata_safely "${OUTPUT_DIR}/temp/encode_standardized_v1.h5ad" "${OUTPUT_DIR}/encode_standardized_v1.h5ad"
# fi
# if [ -f "${OUTPUT_DIR}/temp/gtex_standardized_v1.h5ad" ]; then
#     save_anndata_safely "${OUTPUT_DIR}/temp/gtex_standardized_v1.h5ad" "${OUTPUT_DIR}/gtex_standardized_v1.h5ad"
# fi
# if [ -f "${OUTPUT_DIR}/temp/mage_standardized_v1.h5ad" ]; then
#     save_anndata_safely "${OUTPUT_DIR}/temp/mage_standardized_v1.h5ad" "${OUTPUT_DIR}/mage_standardized_v1.h5ad"
# fi
# if [ -f "${OUTPUT_DIR}/temp/adni_standardized_v1.h5ad" ]; then
#     save_anndata_safely "${OUTPUT_DIR}/temp/adni_standardized_v1.h5ad" "${OUTPUT_DIR}/adni_standardized_v1.h5ad"
# fi

# Step 1.5: Generate Gene ID Reference Mapping with Entrez mapping
log_message "=== Stage 1.5: Generating Gene ID Reference Mapping ==="

# Check if Entrez mapping exists and add it to the command if it does
if [ -f "$ENTREZ_MAPPING" ]; then
    log_message "Using Entrez to Ensembl mapping for enhanced gene ID mapping"
    run_command "python ${SCRIPTS_DIR}/gene_id_mapping_reference.py \\
        --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
        --entex-dir \"${BASE_DIR}/encode/entex\" \\
        --entrez-mapping \"${ENTREZ_MAPPING}\" \\
        --output \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
        ${FORCE_FLAG}"
else
    log_message "Warning: Entrez to Ensembl mapping not found. Proceeding without it."
    run_command "python ${SCRIPTS_DIR}/gene_id_mapping_reference.py \\
        --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
        --entex-dir \"${BASE_DIR}/encode/entex\" \\
        --output \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
        ${FORCE_FLAG}"
fi

# Check if Gene ID Reference Mapping ran successfully
if [ $? -eq 0 ]; then
    log_message "Gene ID Reference Mapping completed successfully!"
else
    log_message "Gene ID Reference Mapping failed. Check the log file for details."
    exit 1
fi

# Step 1.6: Generate ENCODE ID to Ensembl mapping
log_message "=== Stage 1.6: Generating ENCODE ID to Ensembl Mapping ==="
run_command "python ${SCRIPTS_DIR}/generate_encode_mapping.py \\
    --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
    --output-dir \"${METADATA_DIR}/gene_mapping\" \\
    ${FORCE_FLAG}"

# Check if ENCODE ID mapping generation was successful
if [ $? -eq 0 ]; then
    log_message "ENCODE ID mapping generated successfully!"
else
    log_message "Warning: ENCODE ID mapping generation failed. Check the log file for details."
    log_message "Will continue without this mapping, which may affect ENCODE gene ID mapping quality."
fi

# Step 2: Enhanced Metadata Standardization
log_message "=== Stage 2: Enhanced Metadata Standardization ==="
run_command "python ${SCRIPTS_DIR}/standardize_metadata.py \\
    --data-dir \"$OUTPUT_DIR\" \\
    --output-dir \"$OUTPUT_DIR\" \\
    --metadata-dir \"${METADATA_DIR}\""

# Check if Stage 2 ran successfully
if [ $? -eq 0 ]; then
    log_message "Stage 2 completed successfully!"
else
    log_message "Stage 2 failed. Check the log file for details."
    exit 1
fi

# Step 2.5: Preprocess Datasets for Consistent Gene IDs
log_message "=== Stage 2.5: Preprocessing Datasets for Consistent Gene IDs ==="
run_command "python ${SCRIPTS_DIR}/preprocess_dataset_gene_ids.py \\
    --data-dir \"$OUTPUT_DIR\" \\
    --reference-mapping \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
    --output-dir \"$PREPROCESSED_DIR\" \\
    ${FORCE_FLAG}"

# Check if Preprocessing ran successfully
if [ $? -eq 0 ]; then
    log_message "Dataset Preprocessing completed successfully!"
else
    log_message "Dataset Preprocessing failed. Check the log file for details."
    exit 1
fi


# Fix placeholder gene IDs in preprocessed datasets
log_message "=== Step 2.6: Fixing placeholder gene IDs in preprocessed datasets ==="
for dataset in encode gtex mage adni; do
    preprocessed_file="${PREPROCESSED_DIR}/${dataset}_standardized_preprocessed.h5ad"
    if [ -f "$preprocessed_file" ]; then
        log_message "Fixing placeholder IDs in ${dataset} dataset"
        run_command "python ${SCRIPTS_DIR}/fix_placeholder_ids.py $preprocessed_file $preprocessed_file.fixed"
        if [ $? -eq 0 ]; then
            # Replace the original file with the fixed file
            mv "$preprocessed_file.fixed" "$preprocessed_file"
            log_message "Placeholder IDs fixed in ${dataset} dataset"
        else
            log_message "Warning: Failed to fix placeholder IDs in ${dataset} dataset"
        fi
    fi
done

# Step 2.6: Analyze ENCODE mapping quality
log_message "=== Stage 2.6: Analyzing ENCODE Gene Mapping Quality ==="
run_command "python - <<EOF
import scanpy as sc
import pandas as pd
import numpy as np

# Load the ENCODE preprocessed dataset
encode_path = '${PREPROCESSED_DIR}/encode_standardized_preprocessed.h5ad'
try:
    adata = sc.read_h5ad(encode_path)
    print(f\"ENCODE dataset shape: {adata.shape}\")

    # Check ensembl_id column
    if 'ensembl_id' in adata.var.columns:
        # Count non-empty ensembl_ids
        non_empty = sum(1 for x in adata.var['ensembl_id'] if x and str(x).strip() != '')
        percentage = non_empty / len(adata.var) * 100
        print(f\"ENCODE genes with mapped Ensembl IDs: {non_empty}/{len(adata.var)} ({percentage:.2f}%)\")
        
        # Sample of mapped genes
        mapped_ids = adata.var.loc[adata.var['ensembl_id'] != '', 'ensembl_id'].iloc[:5].tolist()
        print(f\"Sample mapped Ensembl IDs: {mapped_ids}\")
        
        # Mapping source distribution
        if 'mapping_source' in adata.var.columns:
            source_counts = adata.var['mapping_source'].value_counts()
            for source, count in source_counts.items():
                source_percentage = count / len(adata.var) * 100
                print(f\"Mapping source '{source}': {count} ({source_percentage:.2f}%)\")
    else:
        print(\"Error: ENCODE dataset does not have an 'ensembl_id' column!\")
    
    # Save mapping stats to file
    mapping_stats = {
        'total_genes': len(adata.var),
        'mapped_genes': non_empty,
        'mapping_percentage': percentage,
        'mapping_sources': {k: int(v) for k, v in source_counts.items()} if 'source_counts' in locals() else {}
    }
    
    import json
    with open('${OUTPUT_DIR}/encode_mapping_stats.json', 'w') as f:
        json.dump(mapping_stats, f, indent=2)
    
    print(f\"ENCODE mapping stats saved to ${OUTPUT_DIR}/encode_mapping_stats.json\")
    
except Exception as e:
    print(f\"Error analyzing ENCODE dataset: {e}\")
    print(\"Continuing with pipeline execution...\")
EOF"



# # Step 3b: Create combined dataset with ALL genes (new approach with improved gene ID mapping)
# log_message "=== Stage 3b: Creating Sparse Combined Dataset with ALL Genes ==="
# run_command "python ${SCRIPTS_DIR}/create_combined_dataset_all_genes_sparse.py \\
#     --input-dir \"$PREPROCESSED_DIR\" \\
#     --reference-mapping \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
#     --output-file \"${OUTPUT_DIR}/combined_all_genes_sparse_standardized.h5ad\" \\
#     --include-datasets \"encode,gtex,mage,adni,entex\" \\
#     ${FORCE_FLAG}"

# # Check if all-genes combined dataset creation was successful
# if [ $? -eq 0 ]; then
#     log_message "Combined dataset (all genes) created successfully!"
# else
#     log_message "Combined dataset (all genes) creation failed. Check the log file for details."
#     # Continue execution even if this fails
# fi

# # Step 3c: Analyze the combined dataset
# log_message "=== Stage 3c: Analyzing Combined Dataset ==="
# run_command "python - <<EOF
# import scanpy as sc
# import pandas as pd
# import numpy as np

# # Load the combined dataset
# combined_path = '${OUTPUT_DIR}/combined_all_genes_sparse_standardized.h5ad'
# try:
#     adata = sc.read_h5ad(combined_path)
#     print(f\"Combined dataset shape: {adata.shape}\")

#     # Check dataset sample distribution
#     print(\"\\nSample distribution by dataset:\")
#     dataset_counts = adata.obs['dataset'].value_counts()
#     for dataset, count in dataset_counts.items():
#         percentage = count / len(adata.obs) * 100
#         print(f\"  - {dataset}: {count} ({percentage:.2f}%)\")

#     # Check gene presence by dataset
#     print(\"\\nGene presence by dataset:\")
#     dataset_gene_counts = {}
#     for dataset in ['adni', 'encode', 'entex', 'gtex', 'mage']:
#         count = sum(1 for x in adata.var['present_in_datasets'] if dataset in str(x))
#         percentage = count / len(adata.var) * 100
#         print(f\"  - {dataset}: {count}/{len(adata.var)} ({percentage:.2f}%)\")
#         dataset_gene_counts[dataset] = count

#     # Check ensembl ID vs. other ID format distribution
#     ensembl_count = sum(1 for x in adata.var_names if str(x).startswith('ENSG'))
#     spike_in_count = sum(1 for x in adata.var_names if str(x).startswith('gSpikein'))
#     other_count = len(adata.var) - ensembl_count - spike_in_count

#     print(\"\\nGene ID format distribution:\")
#     print(f\"  - Ensembl IDs: {ensembl_count}/{len(adata.var)} ({ensembl_count/len(adata.var)*100:.2f}%)\")
#     print(f\"  - Spike-in controls: {spike_in_count}/{len(adata.var)} ({spike_in_count/len(adata.var)*100:.2f}%)\")
#     print(f\"  - Other IDs: {other_count}/{len(adata.var)} ({other_count/len(adata.var)*100:.2f}%)\")
    
#     # Save stats to file
#     combined_stats = {
#         'total_samples': len(adata.obs),
#         'total_genes': len(adata.var),
#         'dataset_distribution': {k: int(v) for k, v in dataset_counts.items()},
#         'gene_presence': dataset_gene_counts,
#         'gene_id_distribution': {
#             'ensembl': int(ensembl_count),
#             'spike_in': int(spike_in_count),
#             'other': int(other_count)
#         }
#     }
    
#     import json
#     with open('${OUTPUT_DIR}/combined_dataset_stats.json', 'w') as f:
#         json.dump(combined_stats, f, indent=2)
    
#     print(f\"Combined dataset stats saved to ${OUTPUT_DIR}/combined_dataset_stats.json\")
    
# except Exception as e:
#     print(f\"Error analyzing combined dataset: {e}\")
#     print(\"Continuing with pipeline execution...\")
# EOF"

# # Step 4: Run validation on all standardized datasets
# log_message "=== Stage 4: Validating Standardized Datasets ==="
# VALIDATION_OUTPUT="${OUTPUT_DIR}/validation_report_${TIMESTAMP}.json"
# run_command "python ${SCRIPTS_DIR}/validate_standardized_datasets.py \\
#     --input-dir \"$OUTPUT_DIR\" \\
#     --output-file \"$VALIDATION_OUTPUT\" \\
#     --file-pattern \"*_standardized*.h5ad\""

# # Check if validation ran successfully
# if [ $? -eq 0 ]; then
#     log_message "Validation completed successfully!"
#     log_message "Validation report saved to $VALIDATION_OUTPUT"
# else
#     log_message "Validation failed. Check the log file for details."
#     # Continue execution even if this fails
# fi


# # Also create a symbolic link to the latest run
# log_message "=== Creating symbolic link to latest run ==="
# LATEST_LINK="${BASE_DIR}/standardized_data/latest"
# run_command "ln -sfn \"$OUTPUT_DIR\" \"$LATEST_LINK\""

# log_message "Complete pipeline executed successfully!"
# log_message "Results saved to $OUTPUT_DIR"
# log_message "Symbolic link to latest run: $LATEST_LINK"

log_message "Complete pipeline executed successfully (Stopped after Stage 2.6)!"
log_message "Individual preprocessed datasets saved to $PREPROCESSED_DIR"
log_message "Log file: $LOG_FILE"
log_message "Individual corrected preprocessed files are available in $PREPROCESSED_DIR"

# # Print a summary of the ENCODE mapping improvement
# log_message "=== ENCODE Mapping Summary ==="
# if [ -f "${OUTPUT_DIR}/encode_mapping_stats.json" ]; then
#     # Extract mapping percentage from the stats file
#     MAPPING_PERCENTAGE=$(grep -o '"mapping_percentage":[^,}]*' "${OUTPUT_DIR}/encode_mapping_stats.json" | cut -d':' -f2)
#     log_message "ENCODE gene mapping percentage: ${MAPPING_PERCENTAGE}%"
# fi

# if [ -f "${OUTPUT_DIR}/combined_dataset_stats.json" ]; then
#     # Extract ENCODE gene presence
#     ENCODE_GENES=$(grep -o '"encode":[^,}]*' "${OUTPUT_DIR}/combined_dataset_stats.json" | head -1 | cut -d':' -f2)
#     TOTAL_GENES=$(grep -o '"total_genes":[^,}]*' "${OUTPUT_DIR}/combined_dataset_stats.json" | head -1 | cut -d':' -f2)
    
#     if [ ! -z "$ENCODE_GENES" ] && [ ! -z "$TOTAL_GENES" ]; then
#         ENCODE_PERCENTAGE=$(echo "scale=2; 100 * $ENCODE_GENES / $TOTAL_GENES" | bc)
#         log_message "ENCODE genes in combined dataset: ${ENCODE_GENES}/${TOTAL_GENES} (${ENCODE_PERCENTAGE}%)"
#     fi
# fi

# log_message "Run the pipeline again with --force flag to regenerate all files"
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/setup_workspace.sh`

````
#!/bin/bash

# Setup script for RNA-seq analysis workspace in Kubernetes environment
# Creates directory structure and installs required packages

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPT_DIR="${BASE_DIR}/scripts"

echo "Setting up RNA-seq analysis workspace..."

# Determine Python path
PYTHON_PATH=$(which python)
echo "Using Python at: ${PYTHON_PATH}"
$PYTHON_PATH --version

# Install pip if needed
echo "Ensuring pip is installed..."
if ! $PYTHON_PATH -m pip --version > /dev/null 2>&1; then
    echo "Installing pip..."
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
    $PYTHON_PATH get-pip.py
    rm get-pip.py
fi

# Install required system dependencies
echo "Installing system dependencies..."
apt-get update && apt-get install -y \
    build-essential \
    gcc \
    python3-dev \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev

# Install required Python packages using the correct Python interpreter
echo "Installing Python packages..."
$PYTHON_PATH -m pip install --no-cache-dir \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    requests \
    scipy \
    pyBigWig \
    scanpy \
    h5py \
    tqdm \
    tables \
    statsmodels \
    openpyxl

# Install anndata separately to ensure it completes successfully
echo "Installing anndata package..."
$PYTHON_PATH -m pip install --no-cache-dir anndata

# Verify installations
echo "Verifying installations..."
$PYTHON_PATH -c "import numpy; print(f'NumPy {numpy.__version__} installed successfully')"
$PYTHON_PATH -c "import pandas; print(f'Pandas {pandas.__version__} installed successfully')"
$PYTHON_PATH -c "import matplotlib; print(f'Matplotlib {matplotlib.__version__} installed successfully')"
$PYTHON_PATH -c "import seaborn; print(f'Seaborn {seaborn.__version__} installed successfully')"
$PYTHON_PATH -c "import requests; print(f'Requests {requests.__version__} installed successfully')"
$PYTHON_PATH -c "import scipy; print(f'SciPy {scipy.__version__} installed successfully')"
$PYTHON_PATH -c "import pyBigWig; print('pyBigWig installed successfully')"
$PYTHON_PATH -c "import h5py; print(f'h5py {h5py.__version__} installed successfully')"
$PYTHON_PATH -c "import tqdm; print(f'tqdm {tqdm.__version__} installed successfully')"
$PYTHON_PATH -c "import openpyxl; print(f'openpyxl {openpyxl.__version__} installed successfully')"

# Verify anndata installation specifically
echo "Verifying anndata installation..."
$PYTHON_PATH -c "
try:
    import anndata
    print(f'anndata {anndata.__version__} installed successfully')
except ImportError as e:
    print(f'Warning: anndata package not installed properly. Error: {e}')
    exit(1)
"

echo "Setup complete!"
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/standardize_datasets.py`

````
#!/usr/bin/env python3
"""
Multi-Dataset Standardization Pipeline

This script processes RNA-seq and expression data from multiple sources 
(ENCODE, GTEx, MAGE, ADNI, ENTEx) into standardized AnnData objects.
It uses JSON configuration files for metadata and ontology mappings.

Usage:
  python standardize_datasets.py --encode-dir /path/to/encode/data \
                                 --gtex-file /path/to/gtex.gct.gz \
                                 --mage-dir /path/to/mage/data \
                                 --adni-dir /path/to/adni/data \
                                 --metadata-dir /path/to/metadata/json \
                                 --output-dir /path/to/output

  # For ENTEx processing:
  python standardize_datasets.py --metadata /path/to/entex_metadata.json \
                                 --dataset-type entex \
                                 --output-dir /path/to/output
"""


import os
import sys
import argparse
import pandas as pd
import numpy as np
import anndata as ad
import glob
import re
import logging
import time
import gzip
import json
from pathlib import Path
import scipy.sparse as sp 

# Import shared utilities
from rnaseq_utils import (
    standardize_ensembl_id,
    load_gencode_mapping,
    add_gencode_annotations,
    standardize_metadata,
    validate_metadata,
    prepare_for_anndata,
    load_mappings,
    load_json_mapping,
    load_dataset_specific_metadata,
    apply_dataset_specific_metadata,
)


# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("data_standardizer")

# Define paths and constants
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
DEFAULT_OUTPUT_DIR = BASE_DIR / "standardized_data"
DEFAULT_METADATA_DIR = BASE_DIR / "metadata/json"

# GTEx metadata files
GTEx_SAMPLE_ATTRS = Path(
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
)
GTEx_SUBJECT_ATTRS = Path(
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
)




def ensure_serializable(obj):
    """
    Recursively convert an object to HDF5-compatible native Python types.
    Handles nested lists, tuples, dicts, numpy arrays, scalars, pandas NA, etc.
    Converts problematic structures (like list of dicts) to JSON strings.
    """
    # --- Handle None first ---
    if obj is None:
        return None

    # --- Handle basic types ---
    if isinstance(obj, (str, int, float, bool)):
        return obj

    # --- Handle numpy/pandas scalars early ---
    # Check specific numpy types before generic hasattr('item')
    if isinstance(obj, np.integer): return int(obj)
    if isinstance(obj, np.floating):
        return None if np.isnan(obj) else float(obj) # Handle NaN directly
    if isinstance(obj, np.bool_): return bool(obj)
    if isinstance(obj, np.void): return None # Or appropriate representation

    # --- Handle numpy arrays specifically ---
    if isinstance(obj, np.ndarray):
        if obj.size == 0: return [] # Empty list for empty array
        # Attempt conversion to list of basic types first
        try:
            if obj.dtype == 'object':
                 # Recursively serialize elements if object type
                 return [ensure_serializable(x) for x in obj.tolist()]
            else:
                 # For numeric/bool arrays, handle internal NaNs then convert
                 if np.issubdtype(obj.dtype, np.floating):
                     # Replace NaN with None representation
                     return [None if np.isnan(x) else float(x) for x in obj.tolist()]
                 elif np.issubdtype(obj.dtype, np.integer):
                     return [int(x) for x in obj.tolist()]
                 elif np.issubdtype(obj.dtype, np.bool_):
                      return [bool(x) for x in obj.tolist()]
                 else:
                     # Other basic types (like fixed-length strings) might be listable directly
                     return obj.tolist()
        except Exception as e:
             logger.warning(f"Could not convert numpy array (dtype: {obj.dtype}) to list: {e}. Trying JSON string.")
             # Fallback to JSON string for complex arrays
             try:
                  # Need a custom encoder for numpy types within JSON
                  class NumpyEncoder(json.JSONEncoder):
                      def default(self, obj):
                          if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                                              np.int16, np.int32, np.int64, np.uint8,
                                              np.uint16, np.uint32, np.uint64)):
                              return int(obj)
                          elif isinstance(obj, (np.float_, np.float16, np.float32,
                                                np.float64)):
                              if np.isnan(obj): return None
                              if np.isinf(obj): return str(obj)
                              return float(obj)
                          elif isinstance(obj, (np.bool_)):
                              return bool(obj)
                          elif isinstance(obj, (np.void)):
                              return None
                          elif isinstance(obj, np.ndarray):
                              # Crucial: Recursively call ensure_serializable here too
                              return ensure_serializable(obj.tolist())
                          return super(NumpyEncoder, self).default(obj)
                  return json.dumps(obj, cls=NumpyEncoder)
             except Exception as json_e:
                  logger.error(f"Could not serialize numpy array to JSON: {json_e}. Falling back to str().")
                  return str(obj)

    # --- Handle simple Iterables (list, tuple, set) ---
    if isinstance(obj, (list, tuple, set)):
        temp_list = list(obj)
        # Check if it *contains* any dictionaries - more robust
        if any(isinstance(item, dict) for item in temp_list):
            logger.debug("List contains dictionaries, converting elements and dumping to JSON string.")
            try:
                # Ensure inner dicts are fully serializable first
                serializable_inner_list = [ensure_serializable(item) for item in temp_list]
                return json.dumps(serializable_inner_list, default=str) # Add default=str for safety
            except Exception as e:
                logger.warning(f"Could not serialize list containing dicts to JSON: {e}. Falling back to str().")
                return str(obj) # Fallback
        else:
            # Otherwise, just serialize elements recursively
            return [ensure_serializable(x) for x in temp_list]

    # --- Handle Dictionaries ---
    if isinstance(obj, dict):
        # Recursively serialize keys (as strings) and values in dicts
        serializable_dict = {}
        for k, v in obj.items():
            serializable_dict[str(k)] = ensure_serializable(v)
        return serializable_dict

    # --- Handle other specific types AFTER basic/numpy/list/dict ---
    if isinstance(obj, pd.Timestamp): return obj.isoformat()
    if isinstance(obj, pd.Series):
        return ensure_serializable(obj.tolist()) # Convert to list first
    if isinstance(obj, pd.Categorical):
        # Convert categories to string list
        return [str(x) for x in obj.tolist()]
    if isinstance(obj, pd.DataFrame):
         logger.warning(f"Attempting to serialize DataFrame in uns. Converting to JSON string.")
         try:
             # Convert to list of records (dicts), ensuring values are serializable
             dict_list = [ensure_serializable(row.to_dict()) for _, row in obj.iterrows()]
             return json.dumps(dict_list, default=str) # Use default=str for safety
         except Exception as e:
              logger.error(f"Failed to convert DataFrame to JSON string: {e}. Using str().")
              return str(obj)

    # --- Robust pd.isna check (LAST check before final string conversion) ---
    # This should only get scalars or types pd.isna understands if logic above is correct
    try:
        # Avoid calling pd.isna on types it explicitly errors on (like multi-element arrays)
        # The checks above should have handled numpy arrays already.
        if not isinstance(obj, (np.ndarray)) and pd.isna(obj):
            return None
    except ValueError as ve:
        # Catch "The truth value of an array is ambiguous" if it still occurs
        logger.warning(f"pd.isna failed for type {type(obj)}: {ve}. Treating as non-NA.")
        pass # Let it fall through to string conversion

    # --- Fallback for any remaining types ---
    logger.warning(f"Converting unrecognized type {type(obj)} to string for serialization.")
    return str(obj)

def identify_columns(df):
    """Identify gene_id and TPM columns across different file formats"""

    # Initialize column mappings
    gene_id_col = None
    tpm_col = None

    # Check column names (case-insensitive)
    columns_lower = [col.lower() for col in df.columns]

    # Check for gene identifier column
    if "gene_id" in df.columns:
        gene_id_col = "gene_id"
    elif "target_id" in df.columns:
        gene_id_col = "target_id"

    # Check for TPM column (case-insensitive)
    if "tpm" in columns_lower:
        tpm_col = df.columns[columns_lower.index("tpm")]

    return gene_id_col, tpm_col


def get_encode_cell_info(metadata_dir=None):
    """Get ENCODE cell line info from JSON."""
    if metadata_dir:
        json_path = os.path.join(metadata_dir, "encode_metadata.json")
        if os.path.exists(json_path):
            try:
                with open(json_path, "r") as f:
                    metadata = json.load(f)
                    # Extract cell line info from the JSON
                    if "dataset_info" in metadata and "cell_lines" in metadata["dataset_info"]:
                        return metadata["dataset_info"]["cell_lines"]
            except Exception as e:
                logger.warning(f"Error loading ENCODE cell info: {e}")

    # Provide a minimal fallback dictionary instead of referring to the removed ENCODE_CELL_INFO
    logger.warning("Could not load cell line info from JSON, using minimal fallback values")
    return {
        "A549": {"tissue": "lung", "organism": "human"},
        "K562": {"tissue": "bone marrow", "organism": "human"},
        "HepG2": {"tissue": "liver", "organism": "human"},
        # Add minimal entries for other cell lines if needed
    }


def create_standard_anndata(data_df, obs_df, var_df, dataset_info):
    # Ensure all obs columns can be properly saved
    for col in obs_df.columns:
        if obs_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting obs column {col} to string")
            obs_df[col] = obs_df[col].astype(str)

    # Ensure all var columns can be properly saved
    for col in var_df.columns:
        if var_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting var column {col} to string")
            var_df[col] = var_df[col].astype(str)

    # Convert all dictionary values to strings
    safe_dataset_info = {}
    for key, value in dataset_info.items():
        if isinstance(value, dict):
            safe_dataset_info[key] = {k: str(v) if v is not None else "" for k, v in value.items()}
        elif value is None:
            safe_dataset_info[key] = ""
        else:
            safe_dataset_info[key] = str(value)
    """
    Create standardized AnnData object with consistent structure.
    Now with enhanced metadata including ontology mappings and validation.

    Parameters:
    -----------
    data_df : pandas.DataFrame
        Expression data (samples x genes)
    obs_df : pandas.DataFrame
        Observation metadata (samples)
    var_df : pandas.DataFrame
        Variable metadata (genes)
    dataset_info : dict
        Dataset information for uns slot

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    # Ensure all obs columns can be properly saved
    for col in obs_df.columns:
        if obs_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting obs column {col} to string")
            obs_df[col] = obs_df[col].astype(str)

    # Ensure all var columns can be properly saved
    for col in var_df.columns:
        if var_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting var column {col} to string")
            var_df[col] = var_df[col].astype(str)

    # Create AnnData object
    adata = ad.AnnData(X=data_df.values, obs=obs_df, var=var_df)

    # Ensure all uns values are serializable
    for key in dataset_info.keys():
        if isinstance(dataset_info[key], (pd.Series, pd.DataFrame, np.ndarray)):
            logger.debug(
                f"Converting uns[{key}] from {type(dataset_info[key])} to standard Python type"
            )
            dataset_info[key] = (
                dataset_info[key].tolist()
                if hasattr(dataset_info[key], "tolist")
                else str(dataset_info[key])
            )

    # Validate metadata before preparing for AnnData
    dataset_name = dataset_info.get("source", "unknown")
    validated_obs_df, validation_report = validate_metadata(obs_df, dataset_name)

    # Prepare data for AnnData
    validated_obs_df, var_df, data_df = prepare_for_anndata(validated_obs_df, var_df, data_df)

    # Create AnnData object (assuming data_df is samples x genes)
    logger.debug(
        f"Creating AnnData with X shape={data_df.values.shape}, obs shape={validated_obs_df.shape}, var shape={var_df.shape}"
    )
    adata = ad.AnnData(X=data_df.values, obs=validated_obs_df, var=var_df)

    # Verify metadata was transferred correctly
    if adata.var.shape[1] == 0:
        logger.error("ERROR: var DataFrame is empty after AnnData creation!")
        logger.error(
            f"Original var_df had shape {var_df.shape} with columns: {list(var_df.columns)}"
        )

        # Emergency fix - manually add the metadata if it was lost
        for col in var_df.columns:
            adata.var[col] = var_df[col].values

    if adata.obs.shape[1] == 0:
        logger.error("ERROR: obs DataFrame is empty after AnnData creation!")
        logger.error(
            f"Original obs_df had shape {validated_obs_df.shape} with columns: {list(validated_obs_df.columns)}"
        )

        # Emergency fix - manually add the metadata if it was lost
        for col in validated_obs_df.columns:
            adata.obs[col] = validated_obs_df[col].values

    # Add dataset information to uns
    # Ensure dataset_info is serializable
    dataset_info = ensure_serializable(dataset_info)
    adata.uns["dataset_info"] = dataset_info

    # Add reference genome info to uns
    adata.uns["reference_genome"] = "hg38"
    adata.uns["gencode_version"] = dataset_info.get("gencode_version", 24)

    # Add validation report to uns
    adata.uns["metadata_validation"] = validation_report

    # Add ontology mapping dictionaries to uns for reference
    mappings = load_mappings()
    adata.uns["ontology_mappings"] = {
        "tissue_to_uberon": mappings.get("TISSUE_TO_UBERON", {}),
        "assay_to_efo": mappings.get("ASSAY_TO_EFO", {}),
        "age_to_hsapdv": mappings.get("AGE_TO_HSAPDV", {}),
        "species_to_ncbi_taxon": mappings.get("SPECIES_TO_NCBI_TAXON", {}),
    }

    # Add processing timestamp
    adata.uns["processing_date"] = pd.Timestamp.now().strftime("%Y-%m-%d")

    return adata



def save_anndata(adata, file_path):
    """Save AnnData object to file with validation and improved serialization."""
    try:
        logger.info(f"Preparing to save AnnData to {file_path}")

        # --- Pre-Save Validation and Cleanup ---
        if adata.n_vars == 0:
            logger.error("Cannot save AnnData: var DataFrame is empty!")
            return False
        if 'gene_id' not in adata.var.columns and len(adata.var_names) > 0:
             logger.warning("Reconstructing potentially missing 'gene_id' column in var from index.")
             adata.var['gene_id'] = adata.var_names
        elif adata.var.shape[1] == 0 and len(adata.var_names) > 0:
             logger.error("Cannot save AnnData: var DataFrame has no columns! Attempting reconstruction.")
             adata.var['gene_id'] = adata.var_names


        # Ensure index names don't conflict with columns
        if adata.obs.index.name is not None and adata.obs.index.name in adata.obs.columns:
            new_name = f"original_obs_{adata.obs.index.name}"
            logger.warning(f"Fixing index/column name conflict in obs: Renaming column '{adata.obs.index.name}' to '{new_name}'")
            adata.obs = adata.obs.rename(columns={adata.obs.index.name: new_name})

        if adata.var.index.name is not None and adata.var.index.name in adata.var.columns:
            new_name = f"original_var_{adata.var.index.name}"
            logger.warning(f"Fixing index/column name conflict in var: Renaming column '{adata.var.index.name}' to '{new_name}'")
            adata.var = adata.var.rename(columns={adata.var.index.name: new_name})

        # Ensure obs/var indices are strings
        adata.obs.index = adata.obs.index.astype(str)
        adata.var.index = adata.var.index.astype(str)

        # --- Explicitly convert object columns in obs/var to string ---
        logger.info("Converting object columns in obs/var to string before saving...")
        for df_name, df in [('obs', adata.obs), ('var', adata.var)]:
            temp_df = df.copy() # Work on a copy to avoid modifying original during iteration
            for col in temp_df.columns:
                if temp_df[col].dtype == 'object':
                    try:
                        # Convert complex objects to string representation safely
                        temp_df[col] = temp_df[col].apply(lambda x: str(x) if isinstance(x, (dict, list)) else x).astype(str).fillna('')
                        logger.debug(f"Converted {df_name} object column '{col}' to string.")
                    except Exception as e:
                        logger.error(f"Failed to convert {df_name} object column '{col}' to string: {e}. Trying str() directly.")
                        try:
                            temp_df[col] = temp_df[col].astype(str).fillna('')
                        except Exception as e2:
                             logger.error(f"Direct str() conversion failed for {df_name} column '{col}': {e2}. Skipping column.")
                             # Optionally drop the column: temp_df = temp_df.drop(columns=[col])
            # Assign potentially modified dataframe back
            if df_name == 'obs':
                 adata.obs = temp_df
            else:
                 adata.var = temp_df
        logger.info("Object column conversion complete.")
        # --- End object column conversion ---

        # --- Serialize .uns ---
        logger.info("Ensuring adata.uns is serializable...")
        original_uns = adata.uns.copy() # Keep a copy of the original uns
        try:
            serializable_uns = ensure_serializable(original_uns)
            adata.uns = serializable_uns
            logger.info("adata.uns converted to serializable format.")
        except Exception as e:
             logger.error(f"Failed during uns serialization: {e}", exc_info=True)
             logger.warning("Attempting to save without full uns serialization.")
             adata.uns = original_uns # Restore original on failure? Maybe not safe.
             # If serialization fails, it's likely the write will fail too.
             # Consider returning False here or letting the write fail.
             # For now, let's let the write attempt proceed but log the warning.
        # --- End .uns serialization ---


        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        # --- Save the AnnData object ---
        logger.info(f"Writing AnnData to {file_path}...")
        save_successful = False # Initialize as False
        try:
            adata.write_h5ad(file_path)
            logger.info(f"Successfully wrote AnnData to {file_path}")
            save_successful = True
        except Exception as write_e:
             logger.error(f"Error writing AnnData to {file_path}: {write_e}", exc_info=True)
             # Explicitly set save_successful to False (already is, but for clarity)
             save_successful = False
        # --- End Save ---

        # Restore original uns in the in-memory object after saving attempt
        adata.uns = original_uns

        # --- Post-Save Verification ---
        if save_successful:
            logger.info(f"Verifying saved file: {file_path}")
            try:
                test_load = ad.read_h5ad(file_path)
                logger.info(
                    f"Verification successful: Loaded AnnData with shape={test_load.shape}, "
                    f"var columns={list(test_load.var.columns)}, obs columns={len(test_load.obs.columns)}"
                )
                if test_load.var.shape[1] == 0:
                     logger.error("Verification FAILED: Loaded AnnData has empty var DataFrame!")
                     return False # Treat as failure
                if test_load.obs.shape[1] == 0:
                     logger.error("Verification FAILED: Loaded AnnData has empty obs DataFrame!")
                     return False # Treat as failure
                return True # Return True only on successful save and verification

            except Exception as load_e:
                logger.error(f"Verification FAILED: Error loading saved file {file_path}: {load_e}")
                return False # Saving failed if we can't reload it
        else:
            # If save_successful is False from the write attempt
            logger.error(f"Save attempt failed for {file_path} due to write error.")
            return False

    except Exception as e:
        logger.error(f"Unhandled error in save_anndata for {file_path}: {e}", exc_info=True)
        # Attempt to restore original uns even on error
        if 'original_uns' in locals() and 'adata' in locals():
             adata.uns = original_uns
        return False # Return False on error
    
def export_metadata_to_json(metadata_df, output_file):
    """
    Export standardized metadata to a JSON file.

    Parameters:
    -----------
    metadata_df : pandas.DataFrame
        Standardized metadata DataFrame
    output_file : str
        Path to save the JSON file

    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Convert DataFrame to dictionary
        metadata_dict = {}

        # Create a record for each sample
        for idx, row in metadata_df.iterrows():
            sample_id = idx if not pd.isna(idx) else f"sample_{len(metadata_dict)}"
            metadata_dict[sample_id] = row.to_dict()

            # Convert numpy and pandas types to Python native types
            for key, value in metadata_dict[sample_id].items():
                if isinstance(value, (np.integer, np.floating)):
                    metadata_dict[sample_id][key] = value.item()
                elif pd.isna(value):
                    metadata_dict[sample_id][key] = None

        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Write to JSON file
        with open(output_file, "w") as f:
            json.dump(metadata_dict, f, indent=2)

        logger.info(f"Exported metadata to {output_file}")
        return True

    except Exception as e:
        logger.error(f"Error exporting metadata: {e}")
        import traceback

        logger.error(traceback.format_exc())
        return False


# ====================================================================================
# ENCODE-specific Processing
# ====================================================================================
def load_entex_metadata(metadata_file=None):
    """
    Load ENTEx metadata from JSON file.

    Parameters:
    -----------
    metadata_file : str or Path, optional
        Path to the ENTEx metadata JSON file

    Returns:
    --------
    dict
        Dictionary containing ENTEx sample metadata and donor information
    """
    try:
        if metadata_file and os.path.exists(metadata_file):
            metadata_path = metadata_file
        else:
            # Use default path if not specified
            metadata_path = Path(
                "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/metadata/entex_metadata.json"
            )

        if not os.path.exists(metadata_path):
            logger.warning(f"ENTEx metadata file not found: {metadata_path}")
            return {"donor_map": {}, "entex_metadata": [], "sample_lookup": {}}

        logger.info(f"Loading ENTEx metadata from {metadata_path}")
        with open(metadata_path, "r") as f:
            metadata = json.load(f)

        # Initialize sample_lookup dictionary
        sample_lookup = {}
        donor_map = metadata.get("donor_map", {})

        # First, load samples from the entex_metadata array if it exists
        if "entex_metadata" in metadata and isinstance(metadata["entex_metadata"], list):
            for sample in metadata["entex_metadata"]:
                sample_id = sample.get("sample_id")
                if sample_id:
                    sample_lookup[sample_id] = sample
            logger.info(f"Loaded {len(sample_lookup)} samples from entex_metadata array")

        # Next, check for individual sample entries at the root level (like ENCFF* keys)
        for key, value in metadata.items():
            if key.startswith("ENCFF") and isinstance(value, dict):
                # This looks like a sample entry directly in the root
                sample_id = key
                sample_lookup[sample_id] = value

                # Ensure it has sample_id field for consistency
                if "sample_id" not in value:
                    sample_lookup[sample_id]["sample_id"] = sample_id

        # Also incorporate any existing sample_lookup dictionary
        if "sample_lookup" in metadata and isinstance(metadata["sample_lookup"], dict):
            for sample_id, sample_info in metadata["sample_lookup"].items():
                sample_lookup[sample_id] = sample_info

        # Update the metadata with the complete sample_lookup
        metadata["sample_lookup"] = sample_lookup

        # If there's no donor_map, try to build it from samples
        if not donor_map:
            donor_map = {}
            for sample_id, info in sample_lookup.items():
                donor_id = info.get("donor", info.get("subject_id", info.get("donor_id")))
                if donor_id and donor_id not in donor_map:
                    donor_map[donor_id] = {"sex": info.get("sex", ""), "age": info.get("age", "")}
            metadata["donor_map"] = donor_map

        logger.info(f"Loaded metadata for {len(sample_lookup)} ENTEx samples")
        return metadata

    except Exception as e:
        logger.error(f"Error loading ENTEx metadata: {e}")
        import traceback

        logger.error(traceback.format_exc())
        return {"donor_map": {}, "entex_metadata": [], "sample_lookup": {}}


def read_encode_tpm_file(file_path):
    """Read an ENCODE TPM file - optimized version."""
    try:
        # Use the low_memory=False option to avoid DtypeWarning
        df = pd.read_csv(file_path, sep="\t", low_memory=False)

        # Use the identify_columns function to identify gene_id and TPM columns
        gene_id_col, tpm_col = identify_columns(df)

        if gene_id_col is not None and tpm_col is not None:
            # Set gene_id as index and handle numeric IDs
            df[gene_id_col] = df[gene_id_col].astype(str)
            # Extract only the columns we need and set index
            result_df = df[[gene_id_col, tpm_col]].set_index(gene_id_col)
            # Rename column to make it consistent
            result_df = result_df.rename(columns={tpm_col: "TPM"})
            return result_df
        else:
            logger.warning(f"Could not identify gene_id and TPM columns in {file_path}")
            logger.warning(f"Available columns: {list(df.columns)}")
            return None

    except Exception as e:
        logger.warning(f"Error reading TPM file {file_path}: {e}")
        return None


def extract_encode_metadata(file_path, entex_metadata=None, metadata_dir=None):
    """
    Extract metadata from ENCODE file path, now using JSON configurations.

    Parameters:
    -----------
    file_path : str
        Path to the ENCODE file
    entex_metadata : dict, optional
        ENTEx metadata with sample_lookup dictionary
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    dict
        Dictionary of metadata
    """
    metadata = {}

    # Get ENCODE cell line info
    encode_cell_info = get_encode_cell_info(metadata_dir)

    # Get file basename and extract file ID
    file_name = os.path.basename(file_path)
    file_id = file_name.split(".")[0]  # e.g., ENCFF008RKC

    # Store original sample ID
    metadata["original_sample_id"] = file_id

    # Check if this is an ENTEx file with available metadata
    if entex_metadata and "sample_lookup" in entex_metadata and "entex" in file_path.lower():
        if file_id in entex_metadata["sample_lookup"]:
            sample_info = entex_metadata["sample_lookup"][file_id]
            donor_id = sample_info.get(
                "donor_id", sample_info.get("donor", sample_info.get("subject_id"))
            )

            # Add basic metadata from lookup
            metadata["tissue"] = sample_info.get("tissue", "")
            metadata["subject_id"] = donor_id
            metadata["assay_type"] = sample_info.get("assay_type", "total")
            metadata["genome_annotation"] = sample_info.get("genome_annotation", "")
            metadata["genome_assembly"] = sample_info.get("genome_assembly", "")
            metadata["data_source"] = "ENTEx"

            # Add donor information if available
            if donor_id and donor_id in entex_metadata.get("donor_map", {}):
                donor_info = entex_metadata["donor_map"][donor_id]
                metadata["sex"] = donor_info.get("sex", "")
                metadata["age"] = donor_info.get("age", "")

            logger.debug(f"Found ENTEx metadata for sample {file_id}")

        else:
            # Extract minimal metadata from the directory structure
            logger.warning(
                f"No metadata found for ENTEx sample {file_id}, extracting from directory"
            )

            # Extract tissue from directory name (e.g., "ENTEx_adrenal gland_unknown")
            path_parts = file_path.split(os.sep)
            entex_parts = [part for part in path_parts if part.startswith("ENTEx_")]

            if entex_parts:
                # Parse tissue from directory name
                tissue_parts = entex_parts[0].split("_")
                if len(tissue_parts) > 1:
                    metadata["tissue"] = tissue_parts[1]

                # Add ENTEx flag
                metadata["data_source"] = "ENTEx"

            # Try to extract donor ID from path
            donor_parts = [part for part in path_parts if part.startswith("entex_tc_")]
            if donor_parts:
                metadata["subject_id"] = donor_parts[0].replace("entex_tc_", "")

    else:
        # Regular ENCODE cell line processing
        # Extract cell line from path - improved detection
        path_parts = file_path.split(os.sep)
        # Check each cell line more thoroughly
        cell_line_dir = []
        for part in path_parts:
            # Check exact match
            if part in encode_cell_info:
                cell_line_dir.append(part)
                continue
            
            # Check for cell line as a substring
            for cl in encode_cell_info.keys():
                if cl in part:
                    cell_line_dir.append(part)
                    break

        if cell_line_dir:
            # Handle directories like 'A549_polyA_plus'
            cell_line_full = cell_line_dir[0]

            # Extract cell line name
            cell_line = None
            for cl in encode_cell_info.keys():
                if cell_line_full.startswith(cl):
                    cell_line = cl
                    break

            if cell_line:
                metadata["cell_line"] = cell_line

                # Add basic cell line information
                if cell_line in encode_cell_info:
                    for key, value in encode_cell_info[cell_line].items():
                        metadata[key] = value

                # Extract RNA extraction method if in directory name
                if "_polyA_plus" in cell_line_full:
                    metadata["extraction_method"] = "polyA_plus"
                elif "_polyA_minus" in cell_line_full:
                    metadata["extraction_method"] = "polyA_minus"
                elif "_total" in cell_line_full:
                    metadata["extraction_method"] = "total"

    # Add standard fields
    metadata["data_type"] = "RNA-seq"
    metadata["expression_unit"] = "TPM"
    
    # Ensure tissue is never None/nan for cell lines
    if "cell_line" in metadata and not metadata.get("tissue"):
        cell_line = metadata["cell_line"]
        if cell_line in encode_cell_info and "tissue" in encode_cell_info[cell_line]:
            metadata["tissue"] = encode_cell_info[cell_line]["tissue"]
        else:
            # Fallback mapping for common cell lines
            fallback_mapping = {
                "A549": "lung",
                "HepG2": "liver", 
                "K562": "blood",
                "HeLa": "cervix",
                "IMR-90": "lung",
                "GM12878": "blood",
                "H1-hESC": "embryonic stem cell"
            }
            if cell_line in fallback_mapping:
                metadata["tissue"] = fallback_mapping[cell_line]

    return metadata


def process_encode_data(
    input_dir, entex_dir=None, output_file=None, entex_metadata_file=None, metadata_dir=None
):
    """
    Process ENCODE RNA-seq data into standardized AnnData, using JSON configurations.

    Parameters:
    -----------
    input_dir : str
        Directory containing ENCODE cell line TPM files
    entex_dir : str, optional
        Directory containing ENCODE ENTEx TPM files
    output_file : str
        Path to save the standardized AnnData
    entex_metadata_file : str, optional
        Path to ENTEx metadata JSON file
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    logger.info(f"Processing ENCODE data from {input_dir}")

    # Load ENTEx metadata
    entex_metadata = (
        load_entex_metadata(entex_metadata_file) if (entex_dir or entex_metadata_file) else None
    )

    # Find all TSV files in cell lines
    tsv_files = []
    if input_dir and os.path.exists(input_dir):
        tsv_files = glob.glob(os.path.join(input_dir, "**", "*.tsv"), recursive=True)

    # Also look for ENTEx files if directory is provided
    if entex_dir:
        logger.info(f"Checking ENTEx directory: {entex_dir}")
        entex_files = glob.glob(os.path.join(entex_dir, "**", "*.tsv"), recursive=True)
        logger.info(f"Found {len(entex_files)} TSV files in ENTEx directory")

        # Check metadata file
        if entex_metadata:
            logger.info(
                f"ENTEx metadata loaded with {len(entex_metadata.get('sample_lookup', {}))} samples"
            )
        else:
            logger.warning("No ENTEx metadata available")

        # Filter for gene quantification files for ENTEx
        if entex_metadata and "sample_lookup" in entex_metadata:
            # Group candidate files by donor and tissue.
            candidate_groups = {}
            for file_path in entex_files:
                file_name = os.path.basename(file_path)
                file_id = file_name.split(".")[0]

                if file_id in entex_metadata["sample_lookup"]:
                    sample_info = entex_metadata["sample_lookup"][file_id]
                    quant_method = sample_info.get("quantification_method", "").lower()
                    # Only consider gene quantification files (exclude transcript quantifications)
                    if "gene" not in quant_method or "transcript" in quant_method:
                        continue

                    # Create group key: use donor (or subject_id) and tissue (normalized to lower case)
                    donor = sample_info.get("donor") or sample_info.get("subject_id", "")
                    tissue = sample_info.get("tissue", "").lower()
                    group_key = (donor, tissue)
                    candidate_groups.setdefault(group_key, []).append(file_path)

            filtered_samples = []
            # Process each (donor, tissue) group
            for group_key, files in candidate_groups.items():
                if len(files) == 1:
                    # Only one candidate – keep it.
                    filtered_samples.append(files[0])
                    logger.debug(f"Keeping single file for group {group_key}: {files[0]}")
                else:
                    # More than one file: try to filter based on annotations.
                    preferred_files = []
                    hg38_files = []
                    for file_path in files:
                        file_name = os.path.basename(file_path)
                        file_id = file_name.split(".")[0]
                        sample_info = entex_metadata["sample_lookup"][file_id]
                        genome_annotation = sample_info.get("genome_annotation", "").upper()
                        genome_assembly = sample_info.get("genome_assembly", "").lower()
                        # Check if file meets both preferred criteria: V24 and hg38.
                        if genome_annotation == "V24" and genome_assembly == "hg38":
                            preferred_files.append(file_path)
                            logger.debug(f"Group {group_key}: file {file_id} meets V24 and hg38")
                        # Otherwise, check if it meets at least the hg38 criterion.
                        elif genome_assembly == "hg38":
                            hg38_files.append(file_path)
                            logger.debug(f"Group {group_key}: file {file_id} meets hg38 only")

                    if preferred_files:
                        # If at least one meets both criteria, use these.
                        filtered_samples.extend(preferred_files)
                    elif hg38_files:
                        # If none meet both, but some meet hg38, use these.
                        filtered_samples.extend(hg38_files)
                    else:
                        # Otherwise, fallback to keeping one candidate (choose the first).
                        filtered_samples.append(files[0])
                        logger.debug(
                            f"Group {group_key}: no file meets hg38; keeping first file: {files[0]}"
                        )

            logger.info(
                f"Selected {len(filtered_samples)} ENTEx gene quantification file(s) after grouping by donor and tissue"
            )
            tsv_files.extend(filtered_samples)

        else:
            # If no metadata filtering is possible, apply basic filtering to file paths
            gene_files = [
                f for f in entex_files if "gene" in f.lower() and "transcript" not in f.lower()
            ]
            logger.info(
                f"Selected {len(gene_files)} ENTEx files that appear to be gene quantifications"
            )
            tsv_files.extend(gene_files)

    if not tsv_files:
        logger.error(f"No TSV files found")
        return None

    logger.info(f"Found {len(tsv_files)} TPM files to process")

    # Process each file
    sample_dfs = {}  # Data frames for each sample
    sample_metadata = {}  # Dictionary for metadata, keyed by sample_id

    # Track all gene IDs
    all_gene_ids = set()

    for file_path in tsv_files:
        # Read TPM data
        tpm_df = read_encode_tpm_file(file_path)
        if tpm_df is None:
            continue

        # Extract sample ID
        file_name = os.path.basename(file_path)
        sample_id = file_name.split(".")[0]

        # Rename TPM column to sample ID
        col_name = tpm_df.columns[0]
        tpm_df = tpm_df.rename(columns={col_name: sample_id})

        # Add to collections
        sample_dfs[sample_id] = tpm_df

        # Extract and store metadata with ENTEx support
        metadata = extract_encode_metadata(file_path, entex_metadata, metadata_dir)
        sample_metadata[sample_id] = metadata

        # Track gene IDs
        all_gene_ids.update(tpm_df.index)

    if not sample_dfs:
        logger.error("No valid TPM data found")
        return None

    # Standardize gene IDs and create expression matrix
    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs")

    # Map original IDs to standardized IDs
    gene_id_mapping = {gene_id: standardize_ensembl_id(gene_id) for gene_id in all_gene_ids}

    # Get unique standardized IDs
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")

    # OPTIMIZED: Create a numpy array to hold the data
    # Shape: (unique_std_ids x sample_ids)
    num_genes = len(unique_std_ids)
    sample_ids = list(sample_dfs.keys())
    num_samples = len(sample_ids)

    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    sample_id_to_idx = {id: idx for idx, id in enumerate(sample_ids)}

    # Pre-allocate the expression matrix with zeros (as floats)
    expr_matrix = np.zeros((num_genes, num_samples), dtype=np.float32)

    # Fill the matrix efficiently
    logger.info("Creating unified expression matrix")
    for sample_id, tpm_df in sample_dfs.items():
        sample_idx = sample_id_to_idx[sample_id]

        # Create temporary dictionary to collect values for this sample
        sample_data = {}

        # Process each gene in this sample
        for gene_id, tpm_val in zip(tpm_df.index, tpm_df[sample_id]):
            if gene_id in gene_id_mapping:
                std_id = gene_id_mapping[gene_id]

                # If multiple original IDs map to same standardized ID, use maximum value
                if std_id in sample_data:
                    sample_data[std_id] = max(sample_data[std_id], tpm_val)
                else:
                    sample_data[std_id] = tpm_val

        # Update the matrix all at once
        for std_id, value in sample_data.items():
            if std_id in std_id_to_idx:  # Handle case if std_id isn't in our target set
                gene_idx = std_id_to_idx[std_id]
                expr_matrix[gene_idx, sample_idx] = value

    # Create the DataFrame from the matrix
    unified_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)

    # Create observation metadata DataFrame - ensuring it matches the sample IDs in unified_df
    # This is critical - the rows in obs_df must exactly match the columns in unified_df
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            # Make a copy of the metadata to avoid modifying the original
            metadata_copy = sample_metadata[sample_id].copy()

            # If there's a sample_id column, rename it to avoid conflict with index
            if "sample_id" in metadata_copy:
                metadata_copy["original_sample_id"] = metadata_copy.pop("sample_id")

            # Add the sample_id to each record
            metadata_copy["_sample_id"] = sample_id
            obs_data.append(metadata_copy)
        else:
            # If we somehow have a sample in the matrix but not in metadata,
            # create minimal metadata for it
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({"_sample_id": sample_id, "dataset": "ENCODE"})

    # Create obs_df with explicit index to ensure alignment
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index("_sample_id", inplace=True)

    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        logger.error(
            f"Metadata has {len(obs_df)} samples, expression matrix has {len(unified_df.columns)} samples"
        )

        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]
        logger.info(f"Fixed alignment - now using {len(obs_df)} samples")

    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df["gene_id"] = var_df.index

    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add original IDs to variable metadata
    var_df["original_ids"] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])))

    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize metadata
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "ENCODE", mappings)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        encode_metadata = load_dataset_specific_metadata(metadata_dir, "encode")
        if encode_metadata:
            # We'll apply this later to the AnnData object
            logger.info("Found ENCODE-specific metadata")

    # Dataset info for uns
    dataset_info = {
        "source": "ENCODE",
        "gencode_version": 24,
        "data_type": "RNA-seq",
        "expression_unit": "TPM",
        "samples": len(obs_df),
        "genes": len(var_df),
    }

    # Final dimension check before creating AnnData
    logger.info(f"Final dimensions check: {len(obs_df)} samples x {len(var_df)} genes")
    logger.info(f"Expression matrix shape: {unified_df.T.shape}")

    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        encode_metadata = load_dataset_specific_metadata(metadata_dir, "encode")
        if encode_metadata:
            adata = apply_dataset_specific_metadata(adata, encode_metadata)

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(
            f"Successfully processed ENCODE data with {adata.n_obs} samples and {adata.n_vars} genes"
        )

    return adata


def process_entex_data(metadata_file, output_file, metadata_dir=None):
    """
    Process ENTEx RNA-seq data into standardized AnnData using metadata.

    Parameters:
    -----------
    metadata_file : str
        Path to ENTEx metadata JSON file
    output_file : str
        Path to save the standardized AnnData
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    logger.info(f"Processing ENTEx data using metadata from {metadata_file}")

    # Load the metadata
    metadata = load_entex_metadata(metadata_file)
    if not metadata or not metadata.get("sample_lookup"):
        logger.error("No valid ENTEx metadata found")
        return None

    sample_lookup = metadata.get("sample_lookup", {})
    logger.info(f"Loaded metadata for {len(sample_lookup)} ENTEx samples")

    # Collect all file paths and process them
    sample_dfs = {}  # Data frames for each sample
    sample_metadata = {}  # Dictionary for metadata, keyed by sample_id
    all_gene_ids = set()  # Track all gene IDs
    processed_count = 0

    for file_id, sample_info in sample_lookup.items():
        # Get file path from metadata
        file_path = sample_info.get("file_path")

        if not file_path:
            logger.warning(f"No file path specified for {file_id}, attempting to locate file")
            # Try to find the file in standard locations
            possible_paths = [
                f"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/entex/{file_id}.tsv",
                f"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/entex/{file_id}.tsv",
            ]

            for path in possible_paths:
                if os.path.exists(path):
                    file_path = path
                    # Update the metadata for future runs
                    sample_info["file_path"] = file_path
                    logger.info(f"Found file for {file_id} at {file_path}")
                    break

            if not file_path:
                logger.warning(f"File not found for {file_id}")
                continue

        if not os.path.exists(file_path):
            logger.warning(f"File not found for {file_id}: {file_path}")
            continue

        # Read TPM data
        tpm_df = read_encode_tpm_file(file_path)
        if tpm_df is None:
            logger.warning(f"Could not read TPM data from {file_path}")
            continue

        # Rename TPM column to sample ID
        col_name = tpm_df.columns[0]
        tpm_df = tpm_df.rename(columns={col_name: file_id})

        # Add to collections
        sample_dfs[file_id] = tpm_df

        # Create metadata dictionary
        donor_id = sample_info.get("donor", sample_info.get("subject_id"))

        metadata_dict = {
            "original_sample_id": file_id,
            "tissue": sample_info.get("tissue", ""),
            "subject_id": donor_id,
            "assay_type": sample_info.get("assay_type", "total"),
            "genome_annotation": sample_info.get("genome_annotation", ""),
            "genome_assembly": sample_info.get("genome_assembly", ""),
            "data_source": "ENTEx",
            "data_type": "RNA-seq",
            "expression_unit": "TPM",
        }

        # Add donor information if available
        donor_map = metadata.get("donor_map", {})
        if donor_id and donor_id in donor_map:
            donor_info = donor_map[donor_id]
            metadata_dict["sex"] = donor_info.get("sex", "")
            metadata_dict["age"] = donor_info.get("age", "")

        # Store metadata
        sample_metadata[file_id] = metadata_dict

        # Track gene IDs
        all_gene_ids.update(tpm_df.index)
        processed_count += 1

        # Log progress
        if processed_count % 10 == 0:
            logger.info(f"Processed {processed_count} ENTEx samples")

    if not sample_dfs:
        logger.error("No valid ENTEx TPM data found")
        return None

    logger.info(f"Successfully processed {processed_count} ENTEx samples")

    # Standardize gene IDs and create expression matrix
    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs")

    # Map original IDs to standardized IDs
    gene_id_mapping = {gene_id: standardize_ensembl_id(gene_id) for gene_id in all_gene_ids}

    # Get unique standardized IDs
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")

    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    sample_ids = list(sample_dfs.keys())
    sample_id_to_idx = {id: idx for idx, id in enumerate(sample_ids)}

    # Pre-allocate the expression matrix
    num_genes = len(unique_std_ids)
    num_samples = len(sample_ids)
    expr_matrix = np.zeros((num_genes, num_samples), dtype=np.float32)

    # Fill the matrix efficiently
    logger.info("Creating unified expression matrix")
    for sample_id, tpm_df in sample_dfs.items():
        sample_idx = sample_id_to_idx[sample_id]

        # Create temporary dictionary to collect values for this sample
        sample_data = {}

        # Process each gene in this sample
        for gene_id, tpm_val in zip(tpm_df.index, tpm_df[sample_id]):
            if gene_id in gene_id_mapping:
                std_id = gene_id_mapping[gene_id]

                # If multiple original IDs map to same standardized ID, use maximum value
                if std_id in sample_data:
                    sample_data[std_id] = max(sample_data[std_id], tpm_val)
                else:
                    sample_data[std_id] = tpm_val

        # Update the matrix all at once
        for std_id, value in sample_data.items():
            if std_id in std_id_to_idx:
                gene_idx = std_id_to_idx[std_id]
                expr_matrix[gene_idx, sample_idx] = value

    # Create the DataFrame from the matrix
    unified_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)

    # Create observation metadata DataFrame
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            metadata_copy = sample_metadata[sample_id].copy()
            metadata_copy["_sample_id"] = sample_id
            obs_data.append(metadata_copy)
        else:
            # If we somehow have a sample in the matrix but not in metadata,
            # create minimal metadata for it
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({"_sample_id": sample_id, "dataset": "ENTEx"})

    # Create obs_df with explicit index
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index("_sample_id", inplace=True)

    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        # Fix alignment
        obs_df = obs_df.loc[unified_df.columns]

    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df["gene_id"] = var_df.index

    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add original IDs to variable metadata
    var_df["original_ids"] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])))

    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize metadata
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "ENTEx", mappings)

    # Dataset info for uns
    dataset_info = {
        "source": "ENTEx",
        "gencode_version": 24,
        "data_type": "RNA-seq",
        "expression_unit": "TPM",
        "samples": len(obs_df),
        "genes": len(var_df),
        "tissues": len(set(obs_df["tissue"])) if "tissue" in obs_df.columns else 0,
        "donors": len(set(obs_df["subject_id"])) if "subject_id" in obs_df.columns else 0,
    }

    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        entex_metadata = load_dataset_specific_metadata(metadata_dir, "entex")
        if entex_metadata:
            adata = apply_dataset_specific_metadata(adata, entex_metadata)

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(
            f"Successfully processed ENTEx data with {adata.n_obs} samples and {adata.n_vars} genes"
        )

    return adata


# ====================================================================================
# GTEx-specific Processing - Enhanced Version
# ====================================================================================


def load_gtex_metadata():
    """Load and standardize GTEx metadata with improved tissue mapping."""
    logger.info("Loading GTEx metadata")

    try:
        # Check if files exist
        if not GTEx_SAMPLE_ATTRS.exists():
            logger.error(f"GTEx sample attributes file not found: {GTEx_SAMPLE_ATTRS}")
            return pd.DataFrame()

        if not GTEx_SUBJECT_ATTRS.exists():
            logger.error(f"GTEx subject attributes file not found: {GTEx_SUBJECT_ATTRS}")
            return pd.DataFrame()

        # Load sample attributes
        logger.info(f"Loading sample attributes from {GTEx_SAMPLE_ATTRS}")
        sample_attrs = pd.read_csv(GTEx_SAMPLE_ATTRS, sep="\t")

        # Load subject attributes
        logger.info(f"Loading subject attributes from {GTEx_SUBJECT_ATTRS}")
        subject_attrs = pd.read_csv(GTEx_SUBJECT_ATTRS, sep="\t")

        # Extract SUBJID from SAMPID
        sample_attrs["SUBJID"] = sample_attrs["SAMPID"].str.extract(r"(GTEX-[A-Z0-9]+)")

        # Merge sample and subject attributes
        merged_attrs = pd.merge(sample_attrs, subject_attrs, on="SUBJID", how="left")

        # Standardize column names
        column_mapping = {
            "SAMPID": "sample_id",
            "SUBJID": "subject_id",
            "SMTSD": "tissue",  # Detailed tissue site
            "SMTS": "broad_tissue",  # Broad tissue type
            "SEX": "sex",
            "AGE": "age",
            "SMGEBTCH": "batch",
            "SMRIN": "rna_integrity_number",
            "SMTSISCH": "ischemic_time",
            "SMNABTCH": "array_batch",
        }

        # Rename columns that exist in the DataFrame
        existing_cols = [col for col in column_mapping.keys() if col in merged_attrs.columns]
        merged_attrs = merged_attrs.rename(
            columns={col: column_mapping[col] for col in existing_cols}
        )

        # Ensure we have a tissue column
        if "tissue" not in merged_attrs.columns and "broad_tissue" in merged_attrs.columns:
            merged_attrs["tissue"] = merged_attrs["broad_tissue"]
            logger.info("Using broad_tissue as tissue column")

        # Set sample ID as index
        merged_attrs.set_index("sample_id", inplace=True)

        logger.info(f"Loaded metadata for {len(merged_attrs)} GTEx samples")

        # Extract tissue statistics
        if "tissue" in merged_attrs.columns:
            tissue_counts = merged_attrs["tissue"].value_counts()
            logger.info(f"Found {len(tissue_counts)} unique tissue types")
            logger.info(f"Top 5 tissues by sample count:")
            for tissue, count in tissue_counts.head(5).items():
                logger.info(f"  - {tissue}: {count} samples")

        return merged_attrs

    except Exception as e:
        logger.error(f"Error loading GTEx metadata: {e}")
        import traceback

        logger.error(traceback.format_exc())
        return pd.DataFrame()


# --- Replace the process_gtex_data function with this memory-optimized version ---
def process_gtex_data(input_file, output_file, metadata_dir=None):
    """
    Process GTEx RNA-seq data into standardized AnnData - MEMORY OPTIMIZED V2.
    Reads GCT line by line, processes genes, and builds final NumPy matrix directly.
    """
    logger.info(f"Processing GTEx data from {input_file} (Memory Optimized V2)")

    # --- 1. Load Metadata ---
    metadata_df = load_gtex_metadata()
    if metadata_df.empty:
        logger.error("Failed to load GTEx metadata. Cannot proceed.")
        return None
    logger.info(f"Loaded metadata for {len(metadata_df)} potential GTEx samples.")

    # --- 2. Read GCT Header and Prepare Structures ---
    sample_ids_in_gct = []
    n_genes_in_gct = 0
    try:
        with gzip.open(input_file, "rt") as f:
            next(f) # Skip version line
            dims = next(f).strip().split()
            n_genes_in_gct = int(dims[0])
            header_line = next(f).strip().split("\t")
            sample_ids_in_gct = header_line[2:] # Skip Name and Description
        logger.info(f"GCT header indicates {n_genes_in_gct} genes and {len(sample_ids_in_gct)} samples.")
    except Exception as e:
        logger.error(f"Error reading GCT header: {e}")
        return None

    # --- 3. Filter Samples: Find intersection between metadata and GCT file ---
    common_samples_set = set(metadata_df.index).intersection(set(sample_ids_in_gct))
    if not common_samples_set:
        logger.error("No common samples found between metadata and GCT file.")
        return None

    final_sample_ids = sorted(list(common_samples_set))
    n_final_samples = len(final_sample_ids)
    logger.info(f"Processing {n_final_samples} common samples.")

    obs_df = metadata_df.loc[final_sample_ids]
    logger.info("Standardizing observation metadata...")
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "GTEx", mappings)
    if "sample_id" in obs_df.columns:
        obs_df = obs_df.rename(columns={"sample_id": "original_sample_id"})
    obs_df['sample_id'] = obs_df.index

    gct_sample_indices = {sid: i for i, sid in enumerate(sample_ids_in_gct)}
    sample_indices_to_keep = [gct_sample_indices[sid] for sid in final_sample_ids]

    # --- 4. First Pass: Identify unique standardized genes and map original IDs ---
    logger.info("First pass: Identifying unique standardized genes...")
    std_gene_id_map = {} # Map original_id -> std_id
    all_std_gene_ids = set() # Collect unique std_ids
    original_id_to_std_map = {} # Map std_id -> list of original_ids

    line_count = 0
    with gzip.open(input_file, "rt") as f:
        next(f); next(f); next(f) # Skip headers
        for line in f:
            line_count += 1
            try:
                original_gene_id = line.strip().split("\t", 1)[0]
                std_gene_id = standardize_ensembl_id(original_gene_id)
                if std_gene_id:
                    std_gene_id_map[original_gene_id] = std_gene_id
                    all_std_gene_ids.add(std_gene_id)
                    if std_gene_id not in original_id_to_std_map:
                        original_id_to_std_map[std_gene_id] = []
                    original_id_to_std_map[std_gene_id].append(original_gene_id)
            except IndexError:
                logger.warning(f"Skipping malformed line {line_count}")
                continue
            if line_count % 10000 == 0:
                logger.info(f"First pass processed {line_count}/{n_genes_in_gct} genes...")

    final_std_gene_ids = sorted(list(all_std_gene_ids))
    n_final_genes = len(final_std_gene_ids)
    std_gene_id_to_idx = {gid: i for i, gid in enumerate(final_std_gene_ids)}
    logger.info(f"Identified {n_final_genes} unique standardized genes.")

    # --- 5. Second Pass: Populate Expression Matrix ---
    logger.info("Second pass: Populating expression matrix...")
    # Initialize matrix with a value indicating "not set yet" (e.g., -1, or NaN if float handling is robust)
    # Using -1 allows easy detection and correct max calculation with non-negative TPMs.
    final_expr_matrix = np.full((n_final_genes, n_final_samples), -1.0, dtype=np.float32)

    line_count = 0
    with gzip.open(input_file, "rt") as f:
        next(f); next(f); next(f) # Skip headers
        for line in f:
            line_count += 1
            try:
                fields = line.strip().split("\t")
                original_gene_id = fields[0]
                std_gene_id = std_gene_id_map.get(original_gene_id) # Use precomputed map

                if not std_gene_id: continue # Skip if not standardizable

                std_idx = std_gene_id_to_idx.get(std_gene_id) # Get the final row index
                if std_idx is None: continue # Should not happen if map is correct

                expression_values_all = fields[2:]
                expression_values_filtered = np.array([float(expression_values_all[i]) for i in sample_indices_to_keep], dtype=np.float32)

                # Update the matrix row using np.maximum, handling the initial -1 value
                current_row = final_expr_matrix[std_idx, :]
                # Only update where the new value is greater OR the current value is the initial -1
                update_mask = (expression_values_filtered > current_row) | (current_row == -1.0)
                final_expr_matrix[std_idx, update_mask] = expression_values_filtered[update_mask]

            except Exception as line_e:
                logger.warning(f"Skipping line {line_count} during second pass due to error: {line_e}")
                continue
            if line_count % 5000 == 0:
                logger.info(f"Second pass processed {line_count}/{n_genes_in_gct} genes...")

    # Replace any remaining -1 values (genes with no expression data after filtering?) with 0
    final_expr_matrix[final_expr_matrix == -1.0] = 0.0
    logger.info("Finished populating expression matrix.")

    # --- 6. Create Var DataFrame ---
    var_df = pd.DataFrame(index=final_std_gene_ids)
    var_df['gene_id'] = var_df.index
    var_df['original_ids'] = [";".join(original_id_to_std_map.get(gid, [])) for gid in final_std_gene_ids]
    logger.info("Adding GENCODE annotations...")
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # --- 7. Create AnnData Object ---
    adata = ad.AnnData(X=sp.csr_matrix(final_expr_matrix.T), obs=obs_df, var=var_df) # Transpose X, store as sparse
    logger.info(f"Created AnnData object with shape {adata.shape}")

    # --- 8. Add UNS Metadata ---
    adata.uns["dataset_info"] = {
        "source": "GTEx", "version": "v10", "gencode_version": 24,
        "data_type": "RNA-seq", "expression_unit": "TPM",
        "samples": adata.n_obs, "genes": adata.n_vars,
        "tissue_count": len(obs_df["tissue"].unique()) if "tissue" in obs_df.columns else 0,
        "subject_count": len(obs_df["subject_id"].unique()) if "subject_id" in obs_df.columns else 0,
    }
    adata.uns["reference_genome"] = "hg38" # Harmonized
    adata.uns["gencode_version"] = 24     # Harmonized
    if metadata_dir:
        gtex_metadata = load_dataset_specific_metadata(metadata_dir, "gtex")
        if gtex_metadata:
            if 'gencode_version' in gtex_metadata: adata.uns['original_gencode_version'] = str(gtex_metadata['gencode_version'])
            if 'reference_genome' in gtex_metadata: adata.uns['original_reference_genome'] = str(gtex_metadata['reference_genome'])
            adata = apply_dataset_specific_metadata(adata, gtex_metadata)
    adata.uns['harmonized_gencode_version'] = '24'
    adata.uns['harmonized_reference_genome'] = 'hg38'
    validated_obs_df, validation_report = validate_metadata(adata.obs, "GTEx")
    adata.obs = validated_obs_df # Update obs
    adata.uns["metadata_validation"] = validation_report
    adata.uns["processing_date"] = pd.Timestamp.now().strftime("%Y-%m-%d")

    # --- 9. Save AnnData ---
    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed GTEx data with {adata.n_obs} samples and {adata.n_vars} genes (Memory Optimized V2)")
    else:
        logger.error("Failed to save memory-optimized GTEx AnnData V2.")
        return None

    return adata

# ====================================================================================
# MAGE-specific Processing
# ====================================================================================


def process_mage_dir(file_path, donor_id, tissue, sample_dfs, sample_metadata, all_gene_ids):
    """
    Process a single MAGE expression file with improved handling for file format.
    Also handles column names properly per standardized format.
    """
    try:
        file_name = os.path.basename(file_path)
        sample_id = f"{donor_id}_{os.path.splitext(file_name)[0]}"

        logger.debug(f"Processing MAGE file: {file_path}")

        # Check if file exists
        if not os.path.exists(file_path):
            logger.warning(f"File does not exist: {file_path}")
            return False

        # Check file size
        file_size = os.path.getsize(file_path)
        if file_size == 0:
            logger.warning(f"Empty file: {file_path}")
            return False

        # Try to read the file - MAGE files are tab-separated even though they have .csv extension
        try:
            # First check the delimiter by reading a few lines
            with open(file_path, "r") as f:
                first_line = f.readline().strip()

            # Determine delimiter based on content
            if "\t" in first_line:
                delimiter = "\t"
            elif "," in first_line:
                delimiter = ","
            else:
                # Default to tab if can't determine
                delimiter = "\t"

            logger.debug(f"Using delimiter '{delimiter}' for file {file_path}")

            # Now read with the appropriate delimiter
            df = pd.read_csv(file_path, delimiter=delimiter)

            # Check if empty after reading
            if df.empty:
                logger.warning(f"Empty data after reading: {file_path}")
                return False

            logger.debug(f"Read data with columns: {df.columns.tolist()}")

            # Check for expected columns
            if "gene_id" in df.columns and "TPM" in df.columns:
                # File has the expected structure, but TPM is actually normalized microarray intensity
                # Rename TPM to normalized_intensity if processing ADNI data
                if "ADNI" in file_path:
                    df = df.rename(columns={"TPM": "normalized_intensity"})
                    expr_col = "normalized_intensity"
                    expr_type = "normalized_intensity"
                else:  # For other datasets like MAGE, keep TPM
                    expr_col = "TPM"
                    expr_type = "TPM"

                # Remove FPKM column if it exists (as per your requirement)
                if "FPKM" in df.columns:
                    df = df.drop("FPKM", axis=1)

                # Create simplified DataFrame with gene_id and expression values
                simple_df = pd.DataFrame()
                simple_df[expr_col] = df[expr_col]
                simple_df.index = df["gene_id"]

                # Store the data
                sample_dfs[sample_id] = simple_df

                # Collect gene IDs
                all_gene_ids.update(df["gene_id"])

                # Check if gene IDs are Ensembl IDs
                is_ensembl = (
                    all(str(gene_id).startswith("ENSG") for gene_id in df["gene_id"].iloc[:5])
                    if len(df) >= 5
                    else False
                )

                # Prepare metadata
                metadata = {
                    "sample_id": sample_id,
                    "donor_id": donor_id,
                    "subject_id": donor_id,
                    "tissue": tissue,
                    "dataset": "MAGE",
                    "data_type": "RNA-seq",
                    "expression_unit": expr_type,  # This will be TPM for RNA-seq or normalized_intensity for microarray
                    "is_gencode": "gencode" in file_name.lower(),
                    "is_ensembl": is_ensembl,
                }

                # Store metadata
                sample_metadata[sample_id] = metadata
                logger.debug(f"Successfully processed {file_path} with {len(simple_df)} genes")
                return True
            else:
                # Unexpected columns - try to identify gene ID and expression columns
                logger.warning(f"Unexpected columns in {file_path}: {df.columns.tolist()}")

                # Look for gene ID column - try common names
                gene_id_col = None
                for col in df.columns:
                    if "gene" in col.lower() and "id" in col.lower():
                        gene_id_col = col
                        break

                # Look for expression column - try common names
                expr_col = None
                for col in df.columns:
                    if col.upper() in [
                        "TPM",
                        "FPKM",
                        "COUNTS",
                        "EXPRESSION",
                        "NORMALIZED_INTENSITY",
                    ]:
                        expr_col = col
                        break

                # If not found, try to use first and second columns as fallback
                if gene_id_col is None and len(df.columns) >= 1:
                    gene_id_col = df.columns[0]
                    logger.debug(f"Using first column '{gene_id_col}' as gene ID column")

                if expr_col is None and len(df.columns) >= 2:
                    expr_col = df.columns[1]
                    logger.debug(f"Using second column '{expr_col}' as expression column")

                # If we have identified columns, process the file
                if gene_id_col is not None and expr_col is not None:
                    # Create simplified DataFrame
                    simple_df = pd.DataFrame()
                    simple_df[expr_col] = df[expr_col]
                    simple_df.index = df[gene_id_col]

                    # Store the data
                    sample_dfs[sample_id] = simple_df

                    # Collect gene IDs
                    all_gene_ids.update(df[gene_id_col])

                    # Check if gene IDs are Ensembl IDs
                    is_ensembl = (
                        all(str(gene_id).startswith("ENSG") for gene_id in df[gene_id_col].iloc[:5])
                        if len(df) >= 5
                        else False
                    )

                    # Determine expression type based on column name
                    if "TPM" in expr_col.upper():
                        expr_type = "TPM"
                    elif "FPKM" in expr_col.upper():
                        expr_type = "FPKM"
                    elif "COUNT" in expr_col.upper():
                        expr_type = "counts"
                    elif "NORMALIZED" in expr_col.upper() or "INTENSITY" in expr_col.upper():
                        expr_type = "normalized_intensity"
                    else:
                        expr_type = "TPM"  # Default

                    # Prepare metadata
                    metadata = {
                        "sample_id": sample_id,
                        "donor_id": donor_id,
                        "subject_id": donor_id,
                        "tissue": tissue,
                        "dataset": "MAGE",
                        "data_type": "RNA-seq",
                        "expression_unit": expr_type,
                        "is_gencode": "gencode" in file_name.lower(),
                        "is_ensembl": is_ensembl,
                    }

                    # Store metadata
                    sample_metadata[sample_id] = metadata
                    logger.debug(f"Successfully processed {file_path} with {len(simple_df)} genes")
                    return True
                else:
                    logger.warning(f"Could not identify required columns in {file_path}")
                    return False

        except Exception as e:
            logger.warning(f"Error reading file {file_path}: {e}")
            return False

    except Exception as e:
        logger.error(f"Error processing MAGE file {file_path}: {e}")
        return False


def process_mage_data(input_dir, output_file, metadata_dir=None):
    """
    Process MAGE RNA-seq data into standardized AnnData - updated to use
    JSON configurations.

    Parameters:
    -----------
    input_dir : str
        Directory containing MAGE expression files
    output_file : str
        Path to save the standardized AnnData
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    logger.info(f"Processing MAGE data from {input_dir}")

    # Track processing statistics
    processed_donors = 0
    processed_samples = 0

    # Create collections to store data
    sample_dfs = {}  # Data frames for each sample
    sample_metadata = {}  # Metadata for each sample
    all_gene_ids = set()  # Track all gene IDs

    # Find all donor directories - includes both NA* and HG* patterns
    donor_dirs = glob.glob(os.path.join(input_dir, "NA*")) + glob.glob(
        os.path.join(input_dir, "HG*")
    )

    if not donor_dirs:
        logger.error(f"No donor directories found in {input_dir}")
        return None

    logger.info(f"Found {len(donor_dirs)} donor directories")

    # Process each donor directory
    for donor_dir in donor_dirs:
        donor_id = os.path.basename(donor_dir)
        processed_donors += 1

        logger.info(f"Processing donor {donor_id}")

        # Check for lymphoblast directory which contains the actual CSV files
        lymphoblast_dir = os.path.join(donor_dir, "lymphoblast")
        if os.path.isdir(lymphoblast_dir):
            tissue_name = "lymphoblast"

            # Find all CSV files in the directory
            csv_files = [f for f in os.listdir(lymphoblast_dir) if f.endswith(".csv")]

            # Separate original and gencode files
            original_files = [f for f in csv_files if "_gencode_" not in f]
            gencode_files = [f for f in csv_files if "_gencode_" in f]

            # Process original files first (prioritize them)
            for file_name in original_files:
                file_path = os.path.join(lymphoblast_dir, file_name)
                if process_mage_dir(
                    file_path, donor_id, tissue_name, sample_dfs, sample_metadata, all_gene_ids
                ):
                    processed_samples += 1

            # Only process gencode files if no original files were successfully processed
            if not any(donor_id in sample_id for sample_id in sample_dfs.keys()):
                for file_name in gencode_files:
                    file_path = os.path.join(lymphoblast_dir, file_name)
                    if process_mage_dir(
                        file_path, donor_id, tissue_name, sample_dfs, sample_metadata, all_gene_ids
                    ):
                        processed_samples += 1
        else:
            # Try to find any tissue directories
            tissue_dirs = [d for d in glob.glob(os.path.join(donor_dir, "*")) if os.path.isdir(d)]

            for tissue_dir in tissue_dirs:
                tissue_name = os.path.basename(tissue_dir)

                # Find all CSV files in the directory
                csv_files = [f for f in os.listdir(tissue_dir) if f.endswith(".csv")]

                # Separate original and gencode files
                original_files = [f for f in csv_files if "_gencode_" not in f]
                gencode_files = [f for f in csv_files if "_gencode_" in f]

                # Process original files first
                for file_name in original_files:
                    file_path = os.path.join(tissue_dir, file_name)
                    if process_mage_dir(
                        file_path, donor_id, tissue_name, sample_dfs, sample_metadata, all_gene_ids
                    ):
                        processed_samples += 1

                # Only process gencode files if no original files were successfully processed for this tissue
                tissue_sample_ids = [
                    sample_id
                    for sample_id in sample_dfs.keys()
                    if donor_id in sample_id and tissue_name in sample_id
                ]
                if not tissue_sample_ids:
                    for file_name in gencode_files:
                        file_path = os.path.join(tissue_dir, file_name)
                        if process_mage_dir(
                            file_path,
                            donor_id,
                            tissue_name,
                            sample_dfs,
                            sample_metadata,
                            all_gene_ids,
                        ):
                            processed_samples += 1

    if not sample_dfs:
        logger.error("No valid expression data files found")
        return None

    logger.info(f"Processed {processed_donors} donors with {processed_samples} samples")

    # Standardize gene IDs and create expression matrix
    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs")

    # Map original IDs to standardized IDs
    gene_id_mapping = {gene_id: standardize_ensembl_id(gene_id) for gene_id in all_gene_ids}

    # Get unique standardized IDs
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")

    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    sample_ids = list(sample_dfs.keys())
    sample_id_to_idx = {id: idx for idx, id in enumerate(sample_ids)}

    # Pre-allocate the expression matrix
    num_genes = len(unique_std_ids)
    num_samples = len(sample_ids)
    expr_matrix = np.zeros((num_genes, num_samples), dtype=np.float32)

    # Fill the matrix efficiently
    logger.info("Creating unified expression matrix")
    for sample_id, expr_df in sample_dfs.items():
        sample_idx = sample_id_to_idx[sample_id]

        # Check if the dataframe is valid
        if expr_df.empty or len(expr_df.columns) == 0:
            logger.warning(f"Skipping empty dataframe for sample {sample_id}")
            continue

        # Create temporary dictionary to collect values for this sample
        sample_data = {}

        try:
            # Get the first column name
            first_col = expr_df.columns[0]

            # Process each gene in this sample
            for gene_id, expr_val in zip(expr_df.index, expr_df[first_col]):
                if gene_id in gene_id_mapping:
                    std_id = gene_id_mapping[gene_id]

                    # If multiple original IDs map to same standardized ID, use maximum value
                    if std_id in sample_data:
                        sample_data[std_id] = max(sample_data[std_id], expr_val)
                    else:
                        sample_data[std_id] = expr_val

            # Update the matrix
            for std_id, value in sample_data.items():
                if std_id in std_id_to_idx:
                    gene_idx = std_id_to_idx[std_id]
                    expr_matrix[gene_idx, sample_idx] = value
        except Exception as e:
            logger.error(f"Error processing sample {sample_id}: {e}")

    # Create the DataFrame from the matrix
    unified_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)

    # Create observation metadata DataFrame
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            metadata_copy = sample_metadata[sample_id].copy()
            metadata_copy["_sample_id"] = sample_id
            obs_data.append(metadata_copy)
        else:
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({"_sample_id": sample_id, "dataset": "MAGE"})

    # Create obs_df with explicit index
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index("_sample_id", inplace=True)

    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]

    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df["gene_id"] = var_df.index

    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add original IDs to variable metadata
    var_df["original_ids"] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])))

    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize metadata
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "MAGE", mappings)

    # Dataset info for uns
    dataset_info = {
        "source": "MAGE",
        "gencode_version": 24,  # Assuming same version as other datasets
        "data_type": "RNA-seq",
        "expression_unit": "TPM",  # Update if different
        "samples": len(obs_df),
        "genes": len(var_df),
        "donor_count": processed_donors,
        "tissue_count": len(set(obs_df["tissue"] if "tissue" in obs_df.columns else [])),
    }

    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        mage_metadata = load_dataset_specific_metadata(metadata_dir, "mage")
        if mage_metadata:
            adata = apply_dataset_specific_metadata(adata, mage_metadata)

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(
            f"Successfully processed MAGE data with {adata.n_obs} samples and {adata.n_vars} genes"
        )

    return adata


# ====================================================================================
# ADNI-specific Processing
# ====================================================================================

def preprocess_adni_file(file_path):
    """Preprocess ADNI file to fix escaped tab characters."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Check if file has escaped tabs
    if '\\t' in content:
        logger.info(f"Fixing escaped tabs in {file_path}")
        # Replace escaped tabs with actual tabs
        fixed_content = content.replace('\\t', '\t')
        
        # Create a temporary fixed file
        fixed_path = file_path + '.fixed'
        with open(fixed_path, 'w') as f:
            f.write(fixed_content)
        
        return fixed_path
    
    return file_path


def preprocess_adni_file(file_path):
    """Preprocess ADNI file to fix escaped tabs."""
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check if file has escaped tabs
        if '\t' in content:
            logger.info(f"Fixing escaped tabs in {file_path}")
            # Replace escaped tabs with actual tabs
            fixed_content = content.replace('\t', '	')
            
            # Create a temporary fixed file
            fixed_path = file_path + '.fixed'
            with open(fixed_path, 'w') as f:
                f.write(fixed_content)
            
            return fixed_path
    except Exception as e:
        logger.warning(f"Error preprocessing file {file_path}: {e}")
    
    return file_path

def process_adni_data(input_dir, output_file, dict_file=None, metadata_dir=None):
    """
    Process ADNI gene expression data into standardized AnnData with JSON configuration support.

    Parameters:
    -----------
    input_dir : str
        Directory containing ADNI sample directories with expression files
    output_file : str
        Path to save the standardized AnnData
    dict_file : str, optional
        Path to gene ID mapping dictionary file
    metadata_dir : str, optional
        Directory containing metadata JSON files

    Returns:
    --------
    anndata.AnnData
        Standardized AnnData object
    """
    logger.info(f"Processing ADNI data from {input_dir}")

    # Track processing statistics
    processed_subjects = 0
    processed_samples = 0

    # Create collections to store data
    sample_dfs = {}  # Data frames for each sample
    sample_metadata = {}  # Metadata for each sample
    all_gene_ids = set()  # Track all gene IDs

    # Find all subject directories
    subject_dirs = glob.glob(os.path.join(input_dir, "*_S_*"))

    if not subject_dirs:
        logger.error(f"No subject directories found in {input_dir}")
        return None

    logger.info(f"Found {len(subject_dirs)} subject directories")

    # Process each subject directory
    for subject_dir in subject_dirs:
        subject_id = os.path.basename(subject_dir)
        processed_subjects += 1

        logger.info(f"Processing subject {subject_id}")

        # Find all CSV files in the directory
        csv_files = glob.glob(os.path.join(subject_dir, "*.csv"))

        if not csv_files:
            logger.warning(f"No CSV files found for subject {subject_id}")
            continue

        # Process each CSV file
        for file_path in csv_files:
            file_name = os.path.basename(file_path)
            sample_id = f"{subject_id}_{os.path.splitext(file_name)[0]}"

            logger.debug(f"Processing file: {file_path}")

            try:

                # Preprocess the file to fix escaped tabs
                fixed_file_path = preprocess_adni_file(file_path)
                try:
                    # Try reading with tab delimiter first
                    df = pd.read_csv(fixed_file_path, sep='\t')
                    # If successful, clean up temporary file if it was created
                    if fixed_file_path != file_path and os.path.exists(fixed_file_path):
                        os.remove(fixed_file_path)
                except Exception as e:
                    # If tab delimiter fails, try comma
                    logger.warning(f"Failed to parse with tab delimiter: {e}")
                    df = pd.read_csv(fixed_file_path)
                    # Clean up temporary file
                    if fixed_file_path != file_path and os.path.exists(fixed_file_path):
                        os.remove(fixed_file_path)

                gene_id_col, tpm_col = identify_columns(df)                
                

                if gene_id_col is None or tpm_col is None:
                    logger.warning(f"Could not identify gene_id and TPM columns in {file_path}")
                    logger.warning(f"Available columns: {list(df.columns)}")
                    continue  # Skip this file

                # Create expression DataFrame using the identified columns
                expr_df = pd.DataFrame(df[tpm_col].values, index=df[gene_id_col])

                # Store the data
                sample_dfs[sample_id] = expr_df

                # Collect gene IDs
                all_gene_ids.update(df[gene_id_col])

                metadata = {
                    "sample_id": sample_id,
                    "subject_id": subject_id,
                    "dataset": "ADNI",
                    "data_type": "Microarray",
                    "expression_unit": "Normalized intensity",
                    "tissue": "blood",
                    "platform": "Affymetrix Human Genome U219 Array",
                    "processing": "gencode_v24" if "gencode_v24" in file_name else "unknown",
                }
                # Store metadata
                sample_metadata[sample_id] = metadata
                processed_samples += 1

                logger.debug(f"Successfully processed {file_path} with {len(expr_df)} genes")

            except Exception as e:
                logger.error(f"Error processing file {file_path}: {e}")
                continue

    if not sample_dfs:
        logger.error("No valid expression data files found")
        return None

    logger.info(f"Processed {processed_subjects} subjects with {processed_samples} samples")

    # Standardize gene IDs and create expression matrix
    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs")

    # Map original IDs to standardized IDs (removing version numbers)
    gene_id_mapping = {gene_id: standardize_ensembl_id(gene_id) for gene_id in all_gene_ids}

    # Get unique standardized IDs
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")

    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    sample_ids = list(sample_dfs.keys())
    sample_id_to_idx = {id: idx for idx, id in enumerate(sample_ids)}

    # Pre-allocate the expression matrix
    num_genes = len(unique_std_ids)
    num_samples = len(sample_ids)
    expr_matrix = np.zeros((num_genes, num_samples), dtype=np.float32)

    # Fill the matrix efficiently
    logger.info("Creating unified expression matrix")
    for sample_id, expr_df in sample_dfs.items():
        sample_idx = sample_id_to_idx[sample_id]

        # Check if the dataframe is valid
        if expr_df.empty or len(expr_df.columns) == 0:
            logger.warning(f"Skipping empty dataframe for sample {sample_id}")
            continue

        # Create temporary dictionary to collect values for this sample
        sample_data = {}

        try:
            # Get the first column name (TPM)
            first_col = expr_df.columns[0]

            # Process each gene in this sample
            for gene_id, expr_val in zip(expr_df.index, expr_df[first_col]):
                if gene_id in gene_id_mapping:
                    std_id = gene_id_mapping[gene_id]

                    # If multiple original IDs map to same standardized ID, use maximum value
                    if std_id in sample_data:
                        sample_data[std_id] = max(sample_data[std_id], expr_val)
                    else:
                        sample_data[std_id] = expr_val

            # Update the matrix
            for std_id, value in sample_data.items():
                if std_id in std_id_to_idx:
                    gene_idx = std_id_to_idx[std_id]
                    expr_matrix[gene_idx, sample_idx] = value
        except Exception as e:
            logger.error(f"Error processing sample {sample_id}: {e}")

    # Create the DataFrame from the matrix
    unified_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)

    # Create observation metadata DataFrame
    obs_data = []
    for sample_id in unified_df.columns:
        if sample_id in sample_metadata:
            metadata_copy = sample_metadata[sample_id].copy()
            metadata_copy["_sample_id"] = sample_id
            obs_data.append(metadata_copy)
        else:
            logger.warning(f"Creating minimal metadata for sample {sample_id}")
            obs_data.append({"_sample_id": sample_id, "dataset": "ADNI"})

    # Create obs_df with explicit index
    obs_df = pd.DataFrame(obs_data)
    obs_df.set_index("_sample_id", inplace=True)

    # Double-check alignment
    if not all(obs_df.index == unified_df.columns):
        logger.error("Metadata index does not match expression matrix columns")
        # Fix alignment if needed
        obs_df = obs_df.loc[unified_df.columns]

    # Create variable metadata DataFrame
    var_df = pd.DataFrame(index=unified_df.index)
    var_df["gene_id"] = var_df.index

    # Create a mapping from standardized IDs back to original IDs
    reverse_mapping = {}
    for orig_id, std_id in gene_id_mapping.items():
        if std_id not in reverse_mapping:
            reverse_mapping[std_id] = []
        reverse_mapping[std_id].append(orig_id)

    # Add original IDs to variable metadata
    var_df["original_ids"] = var_df.index.map(lambda x: ";".join(reverse_mapping.get(x, [])))

    # Load GENCODE mapping and add annotations
    gencode_mapping = load_gencode_mapping()
    var_df = add_gencode_annotations(var_df, gencode_mapping)

    # Standardize metadata
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, "ADNI", mappings)

    # Dataset info for uns
    dataset_info = {
        "source": "ADNI",
        "data_type": "RNA-seq",  # Updated from microarray since these appear to be RNA-seq TPM values
        "gencode_version": 24,
        "expression_unit": "TPM",
        "samples": len(obs_df),
        "genes": len(var_df),
        "subject_count": processed_subjects,
    }

    # Create AnnData object
    adata = create_standard_anndata(unified_df.T, obs_df, var_df, dataset_info)

    # Apply dataset-specific metadata if available
    if metadata_dir:
        adni_metadata = load_dataset_specific_metadata(metadata_dir, "adni")
        if adni_metadata:
            adata = apply_dataset_specific_metadata(adata, adni_metadata)

    # Save the standardized AnnData
    if save_anndata(adata, output_file):
        logger.info(
            f"Successfully processed ADNI data with {adata.n_obs} samples and {adata.n_vars} genes"
        )

    return adata


# ====================================================================================
# Main Processing Pipeline
# ====================================================================================
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Multi-Dataset Standardization Pipeline")
    parser.add_argument("--encode-dir", help="Directory containing ENCODE cell line TPM files")
    parser.add_argument("--encode-entex-dir", help="Directory containing ENCODE ENTEx TPM files")
    parser.add_argument("--gtex-file", help="Path to GTEx expression file (GCT format)")
    parser.add_argument("--mage-dir", help="Directory containing MAGE expression files")
    parser.add_argument(
        "--adni-dir", help="Directory containing ADNI sample directories with expression files"
    )
    parser.add_argument(
        "--metadata-dir",
        default=str(DEFAULT_METADATA_DIR),
        help="Directory containing metadata JSON files",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Output directory for standardized files",
    )

    # Add ENTEx-specific arguments
    parser.add_argument("--entex-metadata-file", help="Path to ENTEx metadata JSON file")
    parser.add_argument(
        "--dataset-type",
        choices=["encode", "entex", "gtex", "mage", "adni"],
        help="Type of dataset to process when using metadata",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    start_time = time.time()

    # Process ENTEx data if metadata file is provided and dataset-type is 'entex'
    if args.entex_metadata_file and args.dataset_type == "entex":
        entex_output = os.path.join(args.output_dir, "entex_standardized_v1.h5ad")
        logger.info(f"Processing ENTEx data from metadata file: {args.entex_metadata_file}")
        entex_data = process_entex_data(args.entex_metadata_file, entex_output, args.metadata_dir)
        if entex_data is not None:
            logger.info(
                f"Successfully processed ENTEx data: {entex_data.n_obs} samples and {entex_data.n_vars} genes"
            )
            # Check gene statistics
            logger.info("ENTEx gene statistics:")
            gencode_stats = entex_data.var["mapping_source"].value_counts()
            for source, count in gencode_stats.items():
                logger.info(f"  {source}: {count} genes ({count/entex_data.n_vars:.1%})")

    # Process standard datasets if their respective inputs are provided
    elif args.encode_dir:
        encode_output = os.path.join(args.output_dir, "encode_standardized_v1.h5ad")
        logger.info(f"Processing ENCODE data from {args.encode_dir}")
        encode_data = process_encode_data(
            args.encode_dir,
            args.encode_entex_dir,
            encode_output,
            args.entex_metadata_file,
            args.metadata_dir,
        )
        if encode_data is not None:
            logger.info(
                f"Successfully processed ENCODE data: {encode_data.n_obs} samples and {encode_data.n_vars} genes"
            )
            # Check gene statistics
            logger.info("ENCODE gene statistics:")
            gencode_stats = encode_data.var["mapping_source"].value_counts()
            for source, count in gencode_stats.items():
                logger.info(f"  {source}: {count} genes ({count/encode_data.n_vars:.1%})")

    if args.gtex_file:
        gtex_output = os.path.join(args.output_dir, "gtex_standardized_v1.h5ad")
        logger.info(f"Processing GTEx data from {args.gtex_file}")
        gtex_data = process_gtex_data(args.gtex_file, gtex_output, args.metadata_dir)
        if gtex_data is not None:
            logger.info(
                f"Successfully processed GTEx data: {gtex_data.n_obs} samples, {gtex_data.n_vars} genes"
            )
            # Check gene statistics
            logger.info("GTEx gene statistics:")
            gencode_stats = gtex_data.var["mapping_source"].value_counts()
            for source, count in gencode_stats.items():
                logger.info(f"  {source}: {count} genes ({count/gtex_data.n_vars:.1%})")

    if args.mage_dir:
        mage_output = os.path.join(args.output_dir, "mage_standardized_v1.h5ad")
        logger.info(f"Processing MAGE data from {args.mage_dir}")
        mage_data = process_mage_data(args.mage_dir, mage_output, args.metadata_dir)
        if mage_data is not None:
            logger.info(
                f"Successfully processed MAGE data: {mage_data.n_obs} samples, {mage_data.n_vars} genes"
            )

    # Process ADNI data
    if args.adni_dir:
        adni_output = os.path.join(args.output_dir, "adni_standardized_v1.h5ad")
        logger.info(f"Processing ADNI data from directory: {args.adni_dir}")
        adni_data = process_adni_data(args.adni_dir, adni_output, None, args.metadata_dir)
        if adni_data is not None:
            logger.info(
                f"Successfully processed ADNI data: {adni_data.n_obs} samples, {adni_data.n_vars} genes"
            )

    # Check for common genes if multiple datasets were processed
    processed_datasets = []
    if args.encode_dir and "encode_data" in locals() and encode_data is not None:
        processed_datasets.append(("ENCODE", encode_data))
    if args.gtex_file and "gtex_data" in locals() and gtex_data is not None:
        processed_datasets.append(("GTEx", gtex_data))
    if args.mage_dir and "mage_data" in locals() and mage_data is not None:
        processed_datasets.append(("MAGE", mage_data))
    if args.adni_dir and "adni_data" in locals() and adni_data is not None:
        processed_datasets.append(("ADNI", adni_data))

    if len(processed_datasets) > 1:
        logger.info("Analyzing gene overlap between datasets:")

        # Compare each pair of datasets
        for i in range(len(processed_datasets)):
            for j in range(i + 1, len(processed_datasets)):
                name1, data1 = processed_datasets[i]
                name2, data2 = processed_datasets[j]

                genes1 = set(data1.var_names)
                genes2 = set(data2.var_names)
                common_genes = genes1.intersection(genes2)

                logger.info(f"  {name1} ({len(genes1)} genes) vs {name2} ({len(genes2)} genes):")
                logger.info(f"    Common genes: {len(common_genes)}")
                logger.info(
                    f"    Overlap: {len(common_genes)/len(genes1):.1%} of {name1}, {len(common_genes)/len(genes2):.1%} of {name2}"
                )

                # Check overlap of highly expressed genes
                # Define highly expressed as genes in the top 10% by mean expression
                def get_top_genes(adata, percent=10):
                    """Calculate top expressed genes based on mean expression."""
                    mean_expr = adata.X.mean(axis=0)

                    # --- START MODIFICATION ---
                    # Convert numpy.matrix (from sparse mean) or ensure 1D array
                    if isinstance(mean_expr, np.matrix):
                        mean_expr = np.array(mean_expr).flatten() # Convert matrix to 1D array
                    elif mean_expr.ndim > 1:
                        mean_expr = mean_expr.flatten() # Flatten if it's somehow multi-dimensional
                    # --- END MODIFICATION ---

                    # Handle cases where all means might be zero or NaN, causing percentile issues
                    if np.all(np.isnan(mean_expr)) or np.all(mean_expr == 0):
                        logger.warning("Mean expression is all NaN or zero. Cannot determine top genes.")
                        return set() # Return empty set

                    # Calculate percentile on potentially NaN-filtered data if necessary
                    # Note: np.percentile handles NaNs by default in recent versions, but explicit handling might be safer
                    valid_mean_expr = mean_expr[~np.isnan(mean_expr)]
                    if len(valid_mean_expr) == 0:
                        logger.warning("Mean expression contains only NaNs. Cannot determine top genes.")
                        return set()

                    cutoff = np.percentile(valid_mean_expr, 100 - percent)

                    # Find indices above cutoff, handling NaNs in the original mean_expr
                    # We need to compare against the original mean_expr to get correct indices
                    top_indices = np.where(np.nan_to_num(mean_expr, nan=-np.inf) >= cutoff)[0]

                    # Ensure indices are valid before slicing var_names
                    valid_indices = [idx for idx in top_indices if idx < len(adata.var_names)]
                    if len(valid_indices) != len(top_indices):
                        logger.warning(f"Found {len(top_indices) - len(valid_indices)} invalid indices in get_top_genes.")

                    return set(adata.var_names[valid_indices])

                top_genes1 = get_top_genes(data1)
                top_genes2 = get_top_genes(data2)
                common_top = top_genes1.intersection(top_genes2)

                logger.info(
                    f"    Top 10% expressed genes: {name1} ({len(top_genes1)}), {name2} ({len(top_genes2)})"
                )
                logger.info(f"    Common top genes: {len(common_top)}")
                logger.info(
                    f"    Top gene overlap: {len(common_top)/len(top_genes1):.1%} of top {name1}, {len(common_top)/len(top_genes2):.1%} of top {name2}"
                )

        # If both ENCODE and GTEx were processed, create a gene compatibility report
        logger.info("Skipping gene compatibility report generation (commented out for performance).") # Add this line
        # if (
        #     "encode_data" in locals()
        #     and "gtex_data" in locals()
        #     and encode_data is not None
        #     and gtex_data is not None
        # ):
        #     report_file = os.path.join(args.output_dir, "gene_compatibility_report.csv")
        #     logger.info(f"Creating gene compatibility report: {report_file}")
        #
        #     # Create DataFrame with gene information from both datasets
        #     encode_genes = pd.DataFrame(index=encode_data.var_names)
        #     encode_genes["gene_name"] = encode_data.var["gene_name"]
        #     encode_genes["gene_type"] = encode_data.var["gene_type"]
        #     encode_genes["in_encode"] = True
        #
        #     gtex_genes = pd.DataFrame(index=gtex_data.var_names)
        #     gtex_genes["gene_name"] = gtex_data.var["gene_name"]
        #     gtex_genes["gene_type"] = gtex_data.var["gene_type"]
        #     gtex_genes["in_gtex"] = True
        #
        #     # Merge the DataFrames
        #     all_genes = encode_genes.join(
        #         gtex_genes, how="outer", lsuffix="_encode", rsuffix="_gtex"
        #     )
        #     # These lines trigger the FutureWarning, they are part of the report generation
        #     all_genes["in_encode"] = all_genes["in_encode"].fillna(False)
        #     all_genes["in_gtex"] = all_genes["in_gtex"].fillna(False)
        #
        #     # Fill in missing gene names using the other dataset
        #     all_genes["gene_name"] = all_genes["gene_name_encode"].combine_first(
        #         all_genes["gene_name_gtex"]
        #     )
        #     all_genes["gene_type"] = all_genes["gene_type_encode"].combine_first(
        #         all_genes["gene_type_gtex"]
        #     )
        #
        #     # Calculate average expression in each dataset
        #     all_genes["mean_expr_encode"] = 0.0
        #     all_genes["mean_expr_gtex"] = 0.0
        #
        #     # ******* THIS LOOP IS THE BOTTLENECK ********
        #     for gene in all_genes.index:
        #         if gene in encode_data.var_names:
        #             gene_idx = encode_data.var_names.get_loc(gene)
        #             # Ensure X is numpy array before mean calculation if sparse
        #             X_encode = encode_data.X.toarray() if sp.issparse(encode_data.X) else encode_data.X # Added sparse check
        #             all_genes.loc[gene, "mean_expr_encode"] = X_encode[:, gene_idx].mean()
        #         if gene in gtex_data.var_names:
        #             gene_idx = gtex_data.var_names.get_loc(gene)
        #             # Ensure X is numpy array before mean calculation if sparse
        #             X_gtex = gtex_data.X.toarray() if sp.issparse(gtex_data.X) else gtex_data.X # Added sparse check
        #             all_genes.loc[gene, "mean_expr_gtex"] = X_gtex[:, gene_idx].mean()
        #     # ******* END BOTTLENECK LOOP ********
        #
        #     # Save the report
        #     all_genes.to_csv(report_file)
        #     logger.info(f"Saved gene compatibility report with {len(all_genes)} genes")
        
    total_time = time.time() - start_time
    logger.info(f"Total processing time: {total_time:.2f} seconds")

    # Print validation summary for all processed datasets
    validation_summary = {}
    if processed_datasets:
        logger.info("=== Metadata Validation Summary ===")
        for name, data in processed_datasets:
            if "metadata_validation" in data.uns:
                report = data.uns["metadata_validation"]
                validation_summary[name] = report

                # Print summary
                logger.info(f"Dataset: {name}")
                logger.info(f"  Samples: {report['sample_count']}")
                logger.info(f"  Status: {report['validation_status']}")

                if report["unmapped_tissues"]:
                    logger.info(f"  Unmapped tissues: {len(report['unmapped_tissues'])}")

                if report["unmapped_assays"]:
                    logger.info(f"  Unmapped assays: {len(report['unmapped_assays'])}")

                if report["unmapped_ages"]:
                    logger.info(f"  Unmapped ages: {len(report['unmapped_ages'])}")

                logger.info(
                    f"  Missing fields: {', '.join(report['missing_fields']) if report['missing_fields'] else 'None'}"
                )

        logger.info("================================")


if __name__ == "__main__":
    main()

````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/standardize_metadata.py`

````
#!/usr/bin/env python3
"""
RNA-seq Metadata Standardization Script (V2 - JSON Prioritized)

This script standardizes metadata for RNA-seq datasets, adding consistent
annotations, ontology mappings, and applying dataset-specific metadata
from JSON configuration files. Inference is used primarily as a fallback.
"""

import os
import re
import logging
import argparse
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import numpy as np
import json # <-- Add json import

# Import shared utilities
from rnaseq_utils import (
    # We might not need standardize_metadata here anymore if Stage 1 is sufficient
    load_dataset_specific_metadata,
    apply_dataset_specific_metadata,
    load_mappings # Keep for ontology lookups if needed
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('metadata_standardization_v2')

# Define paths and constants
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
DEFAULT_METADATA_DIR = BASE_DIR / "metadata/json"
DEFAULT_OUTPUT_DIR = BASE_DIR / "standardized_data"

# --- TissueOntologyMapper Class (Keep as is from previous version) ---
class TissueOntologyMapper:
    """
    A class to handle mapping of tissue terms to standardized ontology terms.
    Loads mappings from JSON files (tissue_to_uberon.json by default).
    """
    def _load_json_mapping(self, mapping_file):
        """Load mappings from a JSON file"""
        try:
            import json
            with open(mapping_file, 'r') as f:
                mapping_data = json.load(f)
                # Handle different JSON structures (simple key:value or key:dict)
                for tissue, ontology_info in mapping_data.items():
                    if isinstance(ontology_info, dict) and 'id' in ontology_info:
                        self.tissue_mappings[tissue.lower()] = ontology_info['id'] # Store keys lowercase
                        self.mapping_confidence[tissue.lower()] = ontology_info.get('confidence', 'medium')
                    elif isinstance(ontology_info, str) and ontology_info.startswith("UBERON:"):
                        self.tissue_mappings[tissue.lower()] = ontology_info # Store keys lowercase
                        self.mapping_confidence[tissue.lower()] = 'medium' # Default confidence
                    else:
                        logger.debug(f"Skipping invalid mapping entry: {tissue}: {ontology_info}")

            logger.info(f"Loaded {len(self.tissue_mappings)} tissue mappings from {mapping_file}")
        except FileNotFoundError:
             logger.warning(f"Tissue mapping file not found: {mapping_file}. Skipping.")
        except Exception as e:
            logger.error(f"Error loading tissue mappings from JSON {mapping_file}: {e}")

    def __init__(self, default_mapping_dir=DEFAULT_METADATA_DIR):
        """
        Initialize the tissue mapper.
        Args:
            default_mapping_dir: Directory containing the default tissue_to_uberon.json
        """
        self.tissue_mappings = {}
        self.mapping_confidence = {}
        self.unmapped_tissues = set() # Track tissues we couldn't map

        # Load the default JSON mapping file
        default_mapping_file = Path(default_mapping_dir) / "tissue_to_uberon.json"
        if default_mapping_file.exists():
            self._load_json_mapping(default_mapping_file)
        else:
            logger.warning(f"Default tissue mapping file not found at {default_mapping_file}")

    def map_tissue(self, tissue_name):
        """
        Map a tissue name to an ontology ID.
        Returns (ontology_id, confidence) tuple.
        """
        # Handle edge cases
        if pd.isna(tissue_name) or not tissue_name or str(tissue_name).strip() == '':
            return '', 'none'

        # Normalize tissue name (lowercase, strip)
        tissue_lower = str(tissue_name).strip().lower()

        # Check exact lowercase match
        if tissue_lower in self.tissue_mappings:
            return (self.tissue_mappings[tissue_lower],
                    self.mapping_confidence.get(tissue_lower, 'medium'))

        # If not found directly, mark as unmapped for now
        self.unmapped_tissues.add(tissue_name) # Store original case
        return '', 'none'

    def map_tissues_in_adata(self, adata, tissue_field='tissue'):
        """Add ontology mappings to an AnnData object based on loaded mappings."""
        if tissue_field not in adata.obs.columns:
            logger.warning(f"Tissue field '{tissue_field}' not found in adata.obs. Skipping ontology mapping.")
            # Ensure columns exist even if mapping is skipped
            if 'tissue_ontology' not in adata.obs:
                adata.obs['tissue_ontology'] = ''
            if 'tissue_ontology_confidence' not in adata.obs:
                adata.obs['tissue_ontology_confidence'] = ''
            return adata

        logger.info(f"Applying tissue ontology mappings using field '{tissue_field}'")

        # Initialize columns if they don't exist, ensuring string type
        if 'tissue_ontology' not in adata.obs.columns:
            adata.obs['tissue_ontology'] = pd.Series(index=adata.obs.index, dtype='str')
        else: # Ensure existing column is string type and fill NaNs
             adata.obs['tissue_ontology'] = adata.obs['tissue_ontology'].astype(str).fillna('')

        if 'tissue_ontology_confidence' not in adata.obs.columns:
            adata.obs['tissue_ontology_confidence'] = pd.Series(index=adata.obs.index, dtype='str')
        else: # Ensure existing column is string type and fill NaNs
             adata.obs['tissue_ontology_confidence'] = adata.obs['tissue_ontology_confidence'].astype(str).fillna('none')

        # --- Efficient Mapping using pandas map ---
        # Create mapping series from dictionaries
        tissue_map_series = pd.Series(self.tissue_mappings)
        confidence_map_series = pd.Series(self.mapping_confidence)

        # Normalize the tissue column for mapping (lowercase, strip)
        normalized_tissue_col = adata.obs[tissue_field].astype(str).str.strip().str.lower()

        # Apply mapping
        adata.obs['tissue_ontology'] = normalized_tissue_col.map(tissue_map_series).fillna('')
        adata.obs['tissue_ontology_confidence'] = normalized_tissue_col.map(confidence_map_series).fillna('none')

        # Identify unmapped tissues (where ontology is empty but original tissue was not)
        unmapped_mask = (adata.obs['tissue_ontology'] == '') & (adata.obs[tissue_field].astype(str).str.strip() != '') & (adata.obs[tissue_field].notna())
        self.unmapped_tissues.update(adata.obs.loc[unmapped_mask, tissue_field].unique())

        # Log stats
        total_samples = adata.n_obs
        mapped_samples = (adata.obs['tissue_ontology'] != '').sum()
        mapping_percentage = (mapped_samples / total_samples) * 100 if total_samples > 0 else 0
        confidence_counts = adata.obs['tissue_ontology_confidence'].value_counts().to_dict()

        logger.info(f"Tissue ontology mapping complete: {mapped_samples}/{total_samples} ({mapping_percentage:.1f}%) mapped.")
        logger.info(f"Mapping confidence counts: {confidence_counts}")

        if self.unmapped_tissues:
             # Log only a sample if there are too many
            sample_unmapped = list(self.unmapped_tissues)[:20]
            logger.warning(f"Found {len(self.unmapped_tissues)} unmapped tissues. Examples: {sample_unmapped}")

        # Convert back to categorical if desired (optional)
        # adata.obs['tissue_ontology'] = pd.Categorical(adata.obs['tissue_ontology'])
        # adata.obs['tissue_ontology_confidence'] = pd.Categorical(adata.obs['tissue_ontology_confidence'])

        return adata

# --- Main Enhancement Function ---

def enhance_standardized_metadata(dataset_name, adata, output_file=None, metadata_dir=None):
    """
    Enhance metadata standardization by prioritizing dataset-specific JSON files.

    Args:
        dataset_name: Name of the dataset (e.g., 'encode', 'gtex')
        adata: AnnData object with dataset loaded from Stage 1 output.
        output_file: Optional file path to save the enhanced standardized data.
        metadata_dir: Directory containing dataset-specific metadata JSON files.

    Returns:
        adata with enhanced standardized metadata, or None if error occurs.
    """
    logger.info(f"Enhancing metadata for {dataset_name} dataset using JSON configurations.")

    # --- 1. Load and Apply Dataset-Specific JSON Metadata ---
    dataset_specific_metadata = None
    if metadata_dir:
        try:
            # Validate the function exists before calling
            if callable(load_dataset_specific_metadata):
                 dataset_specific_metadata = load_dataset_specific_metadata(metadata_dir, dataset_name)
            else:
                 logger.error("load_dataset_specific_metadata function not found in rnaseq_utils.")

            if dataset_specific_metadata:
                logger.info(f"Applying metadata from {dataset_name}_metadata.json")
                # Validate the function exists before calling
                if callable(apply_dataset_specific_metadata):
                    adata = apply_dataset_specific_metadata(adata, dataset_specific_metadata)
                else:
                    logger.error("apply_dataset_specific_metadata function not found in rnaseq_utils.")
            else:
                logger.warning(f"No specific metadata JSON found for {dataset_name}. Will rely on defaults and Stage 1 metadata.")
        except ValueError as ve:
             logger.error(f"Validation error loading metadata for {dataset_name}: {ve}. Proceeding without it.")
        except Exception as e:
             logger.error(f"Could not load or apply metadata JSON for {dataset_name}: {e}")


    # --- 2. Ensure Core Fields Exist and Set Defaults (if not set by JSON) ---
    logger.info("Ensuring core fields and harmonized versions are set...")

    # Harmonization versions (critical for validation)
    if 'harmonized_gencode_version' not in adata.uns or adata.uns.get('harmonized_gencode_version') in [None, 'None', 'unknown']:
        adata.uns['harmonized_gencode_version'] = '24' # Default
        logger.warning(f"Setting default harmonized_gencode_version to '24' for {dataset_name}")
    else: # Ensure format consistency (remove 'v')
        adata.uns['harmonized_gencode_version'] = str(adata.uns['harmonized_gencode_version']).replace('v','')

    if 'harmonized_reference_genome' not in adata.uns or adata.uns.get('harmonized_reference_genome') in [None, 'None', 'unknown']:
        adata.uns['harmonized_reference_genome'] = 'hg38' # Default
        logger.warning(f"Setting default harmonized_reference_genome to 'hg38' for {dataset_name}")
    else: # Allow GRCh38 as well
        if adata.uns['harmonized_reference_genome'] == 'GRCh38':
             adata.uns['harmonized_reference_genome'] = 'hg38' # Standardize to hg38 representation

    # Core obs fields
    if 'dataset' not in adata.obs.columns:
        adata.obs['dataset'] = dataset_name
    adata.obs['dataset'] = adata.obs['dataset'].astype('category') # Ensure categorical

    if 'sample_id' not in adata.obs.columns and adata.obs.index.name != 'sample_id':
        adata.obs['sample_id'] = adata.obs.index.astype(str)
    elif 'sample_id' in adata.obs.columns:
         adata.obs['sample_id'] = adata.obs['sample_id'].astype(str)


    # Basic biological fields (ensure they exist, even if empty or 'unknown')
    for field, default_val, dtype in [
        ('species', 'unknown', 'category'),
        ('species_ontology', '', 'str'),
        ('tissue', 'unknown', 'category'),
        ('tissue_ontology', '', 'str'), # Will be populated by mapper
        ('sex', 'unknown', 'category'),
        ('age', '', 'str'), # Keep age flexible, ontology is separate
        ('developmental_stage_ontology', '', 'str'), # Provide a default if needed
        ('data_type', 'unknown', 'category'),
        ('assay_ontology', '', 'str'),
        ('expression_unit', 'unknown', 'category') # Added expression unit
    ]:
        if field not in adata.obs.columns:
            logger.debug(f"Adding missing obs column '{field}' with default '{default_val}'")
            adata.obs[field] = default_val
        # Ensure correct dtype and handle NaNs appropriately
        if dtype == 'category':
            adata.obs[field] = adata.obs[field].astype('category')
        elif dtype == 'str':
             adata.obs[field] = adata.obs[field].astype(str).fillna(default_val if default_val != 'unknown' else '')


    # --- 3. Apply Tissue Ontology Mapping ---
    if 'tissue' in adata.obs.columns:
        logger.info(f"Applying tissue ontology mapping for {dataset_name}")
        # Assuming default mapping dir contains tissue_to_uberon.json
        tissue_mapper = TissueOntologyMapper(default_mapping_dir=metadata_dir or DEFAULT_METADATA_DIR)
        adata = tissue_mapper.map_tissues_in_adata(adata, tissue_field='tissue')
    else:
        logger.warning(f"'tissue' column not found in {dataset_name}. Skipping tissue ontology mapping.")
        if 'tissue_ontology' not in adata.obs: adata.obs['tissue_ontology'] = ''
        if 'tissue_ontology_confidence' not in adata.obs: adata.obs['tissue_ontology_confidence'] = ''


    # --- 4. Final Checks and Logging ---
    logger.info(f"Final metadata checks for {dataset_name}")

    # Ensure required uns fields are strings for saving
    for key in ['harmonized_gencode_version', 'harmonized_reference_genome',
                'original_gencode_version', 'original_reference_genome',
                'gencode_mapping_notes', 'genome_mapping_notes']:
        if key in adata.uns:
            adata.uns[key] = str(adata.uns[key]) if adata.uns[key] is not None else None

    # Log final state
    log_original_gencode = adata.uns.get('original_gencode_version', 'unknown')
    log_harmonized_gencode = adata.uns.get('harmonized_gencode_version', 'unknown')
    log_original_genome = adata.uns.get('original_reference_genome', 'unknown')
    log_harmonized_genome = adata.uns.get('harmonized_reference_genome', 'unknown')
    protocol_type = adata.uns.get('rna_seq_protocol', {}).get('protocol_type', 'unknown')
    protocol_confidence = adata.uns.get('rna_seq_protocol', {}).get('protocol_confidence', 'low')

    logger.info(f"Metadata enhancement complete for {dataset_name}")
    logger.info(f"  Samples: {adata.n_obs}")
    logger.info(f"  Genes: {adata.n_vars}")
    logger.info(f"  Reference Genome: {log_original_genome} (harmonized to {log_harmonized_genome})")
    logger.info(f"  GENCODE Version: v{log_original_gencode} (harmonized to v{log_harmonized_gencode})")
    logger.info(f"  Protocol Type: {protocol_type} (confidence: {protocol_confidence})")

    # --- 5. Save Standardized Data (if output file specified) ---
    if output_file:
        logger.info(f"Attempting to save standardized {dataset_name} dataset to {output_file}")

        # Ensure necessary columns added by TissueOntologyMapper are strings
        for col in ['tissue_ontology', 'tissue_ontology_confidence']:
            if col in adata.obs.columns:
                 adata.obs[col] = adata.obs[col].astype(str)

        try:
            # Use the robust saving logic including serialization of uns
            from standardize_datasets import save_anndata # Import save function locally
            if callable(save_anndata):
                if save_anndata(adata, output_file):
                    logger.info(f"Successfully saved enhanced dataset to {output_file}")
                else:
                    logger.error(f"Failed to save enhanced dataset {output_file}")
                    return None # Indicate failure
            else:
                logger.error("save_anndata function not found in standardize_datasets.py. Cannot save.")
                return None

        except ImportError:
             logger.error("Could not import save_anndata from standardize_datasets.py. Ensure it's available.")
             return None
        except Exception as e:
            logger.error(f"Critical error during final save of {dataset_name}: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None # Indicate failure

    return adata # Return the enhanced AnnData object

# --- Keep process_all_datasets and main function ---
# (No significant changes needed here, they just call enhance_standardized_metadata)

def process_all_datasets(data_dir, output_dir=None, metadata_dir=None):
    """
    Process all datasets in a directory using enhance_standardized_metadata.

    Args:
        data_dir: Directory containing datasets (expects *_standardized_v1.h5ad)
        output_dir: Optional directory to save standardized datasets (as *_standardized_v2.h5ad)
        metadata_dir: Directory containing dataset-specific metadata JSON files
    """
    data_dir = Path(data_dir)

    if not output_dir:
        output_dir = data_dir
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)

    # Prioritize v1 files from Stage 1
    h5ad_files = list(data_dir.glob("*_standardized_v1.h5ad"))

    if not h5ad_files:
        logger.warning(f"No *_standardized_v1.h5ad files found in {data_dir}. Looking for any *.h5ad...")
        h5ad_files = list(data_dir.glob("*.h5ad"))
        # Avoid processing v2 files if they exist
        h5ad_files = [f for f in h5ad_files if '_standardized_v2' not in f.name]

    if not h5ad_files:
        logger.error(f"No suitable h5ad files found in {data_dir} to process.")
        return

    logger.info(f"Found {len(h5ad_files)} h5ad files to process in Stage 2")

    processed_count = 0
    failed_count = 0
    for h5ad_file in h5ad_files:
        # Extract dataset name (e.g., 'encode' from 'encode_standardized_v1.h5ad')
        dataset_name = h5ad_file.stem.split('_')[0].lower()

        try:
            logger.info(f"--- Processing dataset: {dataset_name} ---")
            logger.info(f"Loading {dataset_name} dataset from {h5ad_file}")
            adata = sc.read_h5ad(h5ad_file)

            # Define output file path for v2
            output_file = output_dir / f"{dataset_name}_standardized_v2.h5ad"

            # Enhance metadata standardization
            enhanced_adata = enhance_standardized_metadata(
                dataset_name,
                adata,
                output_file=output_file, # Pass output path for saving
                metadata_dir=metadata_dir
            )

            if enhanced_adata is not None:
                 processed_count += 1
            else:
                 failed_count += 1
                 logger.error(f"Failed to enhance metadata for {dataset_name}")


        except Exception as e:
            failed_count += 1
            logger.error(f"Error processing {dataset_name} from {h5ad_file}: {e}")
            import traceback
            logger.error(traceback.format_exc())

    logger.info(f"--- Stage 2 Summary ---")
    logger.info(f"Successfully processed: {processed_count} datasets")
    logger.info(f"Failed to process: {failed_count} datasets")
    logger.info("----------------------")


def main():
    parser = argparse.ArgumentParser(description="Standardize metadata for RNA-seq datasets (V2 - JSON Prioritized)")
    parser.add_argument("--data-dir", required=True, help="Directory containing *_standardized_v1.h5ad files")
    parser.add_argument("--output-dir", help="Directory to save *_standardized_v2.h5ad datasets (defaults to data-dir)")
    # Removed mapping-dir as TissueOntologyMapper now uses default path
    parser.add_argument("--metadata-dir", default=str(DEFAULT_METADATA_DIR),
                         help="Directory containing dataset-specific metadata JSON files (e.g., encode_metadata.json)")
    parser.add_argument("--dataset", help="Process only this specific dataset name (e.g., 'encode')")

    args = parser.parse_args()

    output_dir = args.output_dir if args.output_dir else args.data_dir

    if args.dataset:
        # Process a single specified dataset
        dataset_name = args.dataset.lower()
        data_dir = Path(args.data_dir)
        # Look for the v1 file first
        h5ad_file = data_dir / f"{dataset_name}_standardized_v1.h5ad"

        if not h5ad_file.exists():
            logger.warning(f"V1 file not found for {dataset_name}, looking for any {dataset_name}*.h5ad...")
            found_files = list(data_dir.glob(f"{dataset_name}*.h5ad"))
            # Exclude v2 files specifically
            found_files = [f for f in found_files if '_standardized_v2' not in f.name]
            if not found_files:
                 logger.error(f"No suitable h5ad file found for dataset {dataset_name} in {data_dir}")
                 return
            h5ad_file = found_files[0] # Process the first one found

        logger.info(f"Processing single dataset: {dataset_name} from {h5ad_file}")
        try:
            adata = sc.read_h5ad(h5ad_file)
            output_file_path = Path(output_dir) / f"{dataset_name}_standardized_v2.h5ad"
            enhance_standardized_metadata(dataset_name, adata, output_file=output_file_path, metadata_dir=args.metadata_dir)
        except Exception as e:
            logger.error(f"Error processing single dataset {dataset_name}: {e}")
            import traceback
            logger.error(traceback.format_exc())
    else:
        # Process all datasets found in the directory
        process_all_datasets(args.data_dir, output_dir, args.metadata_dir)

if __name__ == "__main__":
    main()
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v1/validate_standardized_datasets.py`

````
#!/usr/bin/env python3
"""
Standalone Validation Script for Standardized RNA-seq Datasets

This script validates standardized RNA-seq datasets against the required standards
and generates a comprehensive report.

Usage:
  python validate_standardized_datasets.py --input-dir /path/to/standardized/data \
                                          --output-file /path/to/validation_report.json
"""

import os
import sys
import json
import argparse
import logging
import scanpy as sc
import pandas as pd
from pathlib import Path
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('dataset_validator')

def extract_dataset_type(file_path):
    """Extract dataset type from filename."""
    file_name = os.path.basename(file_path)
    parts = file_name.split('_')
    if len(parts) > 0:
        return parts[0].lower()
    return "unknown"

def validate_dataset(file_path):
    """
    Validate a standardized dataset against requirements.
    
    Args:
        file_path: Path to the standardized dataset file
        
    Returns:
        Dictionary with validation results
    """
    try:
        file_path = Path(file_path)
        dataset_type = extract_dataset_type(file_path)
        logger.info(f"Validating {dataset_type} dataset: {file_path}")
        
        if not file_path.exists():
            logger.error(f"File not found: {file_path}")
            return {
                'dataset': dataset_type,
                'file': str(file_path),
                'status': 'error',
                'message': 'File not found',
                'validations': {}
            }
        
        # Load the dataset
        adata = sc.read_h5ad(file_path)
        
        # Initialize validation results
        validation_results = {
            'dataset': dataset_type,
            'file': str(file_path),
            'status': 'passed',
            'n_samples': adata.n_obs,
            'n_genes': adata.n_vars,
            'validations': {}
        }
        
        # Check harmonized GENCODE version
        gencode_version = adata.uns.get('harmonized_gencode_version', None)
        if gencode_version is not None:
            gencode_version = str(gencode_version).replace('v', '')
            if gencode_version == '24':
                validation_results['validations']['gencode_version'] = {
                    'status': 'passed',
                    'value': gencode_version
                }
            else:
                validation_results['validations']['gencode_version'] = {
                    'status': 'failed',
                    'value': gencode_version,
                    'expected': '24'
                }
                validation_results['status'] = 'failed'
        else:
            validation_results['validations']['gencode_version'] = {
                'status': 'missing',
                'expected': '24'
            }
            validation_results['status'] = 'failed'
        
        # Check harmonized reference genome
        genome_version = adata.uns.get('harmonized_reference_genome', None)
        if genome_version is not None:
            if genome_version in ['hg38', 'GRCh38']:
                validation_results['validations']['reference_genome'] = {
                    'status': 'passed',
                    'value': genome_version
                }
            else:
                validation_results['validations']['reference_genome'] = {
                    'status': 'failed',
                    'value': genome_version,
                    'expected': 'hg38/GRCh38'
                }
                validation_results['status'] = 'failed'
        else:
            validation_results['validations']['reference_genome'] = {
                'status': 'missing',
                'expected': 'hg38/GRCh38'
            }
            validation_results['status'] = 'failed'
        
        # Check observation metadata fields
        metadata_fields = {
            'tissue': {
                'ontology_field': 'tissue_ontology',
                'ontology_prefix': 'UBERON:',
                'importance': 'critical'
            },
            'sex': {
                'values': ['male', 'female', 'unknown'],
                'importance': 'important'
            },
            'species': {
                'ontology_field': 'species_ontology',
                'ontology_prefix': 'NCBITaxon:',
                'importance': 'important'
            },
            'data_type': {
                'values': ['RNA-seq', 'microarray'],
                'importance': 'important'
            },
            'assay_ontology': {
                'ontology_prefix': 'EFO:',
                'importance': 'important'
            },
            'cell_type': {
                'ontology_prefix': 'CL:',
                'importance': 'important'
            }
        }


        
        for field, config in metadata_fields.items():
            if field in adata.obs.columns:
                # Check for missing values
                missing_count = adata.obs[field].isna().sum()
                missing_percentage = (missing_count / adata.n_obs) * 100
                
                # Check for ontology field
                if 'ontology_field' in config and config['ontology_field'] in adata.obs.columns:
                    # Check if ontology values are valid
                    ontology_field = config['ontology_field']
                    ontology_prefix = config.get('ontology_prefix', '')
                    
                    # Count values with correct prefix
                    valid_values = adata.obs[ontology_field].astype(str).str.startswith(ontology_prefix)
                    valid_count = valid_values.sum()
                    valid_percentage = (valid_count / adata.n_obs) * 100
                    
                    validation_results['validations'][field] = {
                        'status': 'passed' if valid_percentage >= 90 else 'warning' if valid_percentage >= 70 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'valid_percentage': float(valid_percentage),
                        'has_ontology': True
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
                
                # Check for enumerated values
                elif 'values' in config:
                    valid_values = adata.obs[field].isin(config['values'])
                    valid_count = valid_values.sum()
                    valid_percentage = (valid_count / adata.n_obs) * 100
                    
                    validation_results['validations'][field] = {
                        'status': 'passed' if valid_percentage >= 90 else 'warning' if valid_percentage >= 70 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'valid_percentage': float(valid_percentage),
                        'has_ontology': False
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
                
                # Simple presence check
                else:
                    validation_results['validations'][field] = {
                        'status': 'passed' if missing_percentage <= 10 else 'warning' if missing_percentage <= 30 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'has_ontology': False
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
            else:
                validation_results['validations'][field] = {
                    'status': 'missing',
                    'importance': config['importance']
                }
                
                if config['importance'] == 'critical':
                    validation_results['status'] = 'failed'
        
        # Additional validation: Check gene IDs format
        if adata.n_vars > 0:
            # Check if gene IDs follow Ensembl format
            sample_genes = adata.var_names[:10].tolist()
            ensembl_format = all(str(g).startswith('ENSG') for g in sample_genes)
            
            validation_results['validations']['gene_id_format'] = {
                'status': 'passed' if ensembl_format else 'failed',
                'value': 'Ensembl' if ensembl_format else 'Unknown'
            }
            
            if not ensembl_format and validation_results['status'] == 'passed':
                validation_results['status'] = 'warning'
        
        logger.info(f"Validation completed for {dataset_type}: {validation_results['status']}")
        return validation_results
        
    except Exception as e:
        logger.error(f"Error validating dataset {file_path}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {
            'dataset': extract_dataset_type(file_path),
            'file': str(file_path),
            'status': 'error',
            'message': str(e),
            'validations': {}
        }

def generate_report(validation_results, output_file):
    """
    Generate a validation report from the results.
    
    Args:
        validation_results: List of validation result dictionaries
        output_file: Path to save the report
        
    Returns:
        Path to the report file
    """
    try:
        output_file = Path(output_file)
        
        # Create summary statistics
        summary = {
            'timestamp': datetime.now().isoformat(),
            'datasets_validated': len(validation_results),
            'datasets_passed': sum(1 for r in validation_results if r['status'] == 'passed'),
            'datasets_warning': sum(1 for r in validation_results if r['status'] == 'warning'),
            'datasets_failed': sum(1 for r in validation_results if r['status'] == 'failed'),
            'datasets_error': sum(1 for r in validation_results if r['status'] == 'error'),
            'total_samples': sum(r.get('n_samples', 0) for r in validation_results),
            'total_genes': sum(r.get('n_genes', 0) for r in validation_results),
            'dataset_results': validation_results
        }
        
        # Ensure output directory exists
        output_file.parent.mkdir(exist_ok=True, parents=True)
        
        # Write report to file
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Also generate a more human-readable HTML report
        html_file = output_file.with_suffix('.html')
        
        # Simple HTML report
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>RNA-seq Standardization Validation Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                h1, h2, h3 {{ color: #333; }}
                .summary {{ background-color: #f5f5f5; padding: 10px; border-radius: 5px; }}
                .dataset {{ margin-bottom: 20px; padding: 10px; border: 1px solid #ddd; }}
                .passed {{ background-color: #dff0d8; }}
                .warning {{ background-color: #fcf8e3; }}
                .failed {{ background-color: #f2dede; }}
                .error {{ background-color: #f5f5f5; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <h1>RNA-seq Standardization Validation Report</h1>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            
            <div class="summary">
                <h2>Summary</h2>
                <p>Datasets validated: {summary['datasets_validated']}</p>
                <p>Datasets passed: {summary['datasets_passed']}</p>
                <p>Datasets with warnings: {summary['datasets_warning']}</p>
                <p>Datasets failed: {summary['datasets_failed']}</p>
                <p>Datasets with errors: {summary['datasets_error']}</p>
                <p>Total samples: {summary['total_samples']}</p>
                <p>Total genes: {summary['total_genes']}</p>
            </div>
            
            <h2>Dataset Results</h2>
        """
        
        # Add sections for each dataset
        for result in validation_results:
            status_class = result['status']
            dataset_name = result['dataset'].upper()
            
            html_content += f"""
            <div class="dataset {status_class}">
                <h3>{dataset_name} - {result['status'].upper()}</h3>
                <p>File: {result.get('file', 'N/A')}</p>
                <p>Samples: {result.get('n_samples', 'N/A')}</p>
                <p>Genes: {result.get('n_genes', 'N/A')}</p>
            """
            
            if result['status'] == 'error':
                html_content += f"""
                <p>Error: {result.get('message', 'Unknown error')}</p>
                """
            else:
                html_content += """
                <h4>Validation Results</h4>
                <table>
                    <tr>
                        <th>Field</th>
                        <th>Status</th>
                        <th>Details</th>
                    </tr>
                """
                
                for field, details in result.get('validations', {}).items():
                    status = details.get('status', 'unknown')
                    
                    if status == 'passed':
                        details_str = f"Value: {details.get('value', 'N/A')}"
                    elif status == 'failed':
                        details_str = f"Value: {details.get('value', 'N/A')}, Expected: {details.get('expected', 'N/A')}"
                    elif status == 'missing':
                        details_str = f"Expected: {details.get('expected', 'N/A')}"
                    else:
                        # For fields with percentages
                        if 'missing_percentage' in details:
                            details_str = f"Missing: {details['missing_percentage']:.1f}%"
                        
                        if 'valid_percentage' in details:
                            details_str += f", Valid: {details['valid_percentage']:.1f}%"
                    
                    html_content += f"""
                    <tr>
                        <td>{field}</td>
                        <td>{status.upper()}</td>
                        <td>{details_str}</td>
                    </tr>
                    """
                
                html_content += """
                </table>
                """
            
            html_content += """
            </div>
            """
        
        html_content += """
        </body>
        </html>
        """
        
        # Write HTML report
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Report saved to {output_file}")
        logger.info(f"HTML report saved to {html_file}")
        
        return output_file
        
    except Exception as e:
        logger.error(f"Error generating report: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Validate Standardized RNA-seq Datasets')
    
    parser.add_argument('--input-dir', required=True, help='Directory containing standardized datasets')
    parser.add_argument('--output-file', help='Path to save the validation report JSON file')
    parser.add_argument('--file-pattern', default='*_standardized_v2.h5ad', 
                       help='File pattern to match standardized datasets')
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Find standardized datasets
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        sys.exit(1)
    
    # Determine output file
    if args.output_file:
        output_file = Path(args.output_file)
    else:
        output_file = input_dir / f"validation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    # Find datasets to validate
    h5ad_files = list(input_dir.glob(args.file_pattern))
    if not h5ad_files:
        logger.error(f"No standardized datasets found matching pattern '{args.file_pattern}' in {input_dir}")
        sys.exit(1)
    
    logger.info(f"Found {len(h5ad_files)} datasets to validate")
    
    # Validate each dataset
    validation_results = []
    for h5ad_file in h5ad_files:
        result = validate_dataset(h5ad_file)
        validation_results.append(result)
    
    # Generate report
    report_file = generate_report(validation_results, output_file)
    if report_file:
        logger.info(f"Validation completed. Report saved to {report_file}")
    else:
        logger.error("Failed to generate validation report")
        sys.exit(1)

if __name__ == '__main__':
    main()
````

