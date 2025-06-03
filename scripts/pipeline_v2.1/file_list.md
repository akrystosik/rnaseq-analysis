# File List with Contents

The following files were found, along with their contents:

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/anndata_save_wrapper.py`

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
import json # For complex object serialization to string

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('anndata_save_wrapper')

def convert_to_serializable(obj):
    """Convert dict values to HDF5-compatible native Python types or JSON strings."""
    if isinstance(obj, dict):
        return {convert_to_serializable(k): convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)): # Handle tuples as lists
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, (np.integer, int)):
        return int(obj)
    elif isinstance(obj, (np.floating, float)):
        if np.isnan(obj):
            return None  # Represent NaN as None, HDF5 handles None
        return float(obj)
    elif isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        if obj.dtype == 'object':
            # If array of objects, recursively serialize each element
            return [convert_to_serializable(x) for x in obj.tolist()]
        elif obj.size == 0: # Handle empty array
            return []
        else:
            # For numeric/boolean arrays, tolist() is generally safe
            # NaNs in float arrays become Python float NaNs, handled by HDF5
            return obj.tolist()
    elif isinstance(obj, (pd.Series, pd.DataFrame)):
        logger.warning(f"Attempting to serialize pandas {type(obj).__name__} in uns. Converting to JSON string.")
        try:
            # For DataFrames, convert to list of records (dicts), then serialize
            if isinstance(obj, pd.DataFrame):
                return json.dumps([convert_to_serializable(record) for record in obj.to_dict(orient='records')], default=str)
            # For Series, convert to list, then serialize
            return json.dumps(convert_to_serializable(obj.tolist()), default=str)
        except Exception as e:
            logger.error(f"Could not serialize pandas {type(obj).__name__} to JSON: {e}. Falling back to str().")
            return str(obj)
    elif pd.isna(obj): # Check for pandas NA types (like pd.NA, NaT) AFTER specific type checks
        return None # Represent as None
    elif isinstance(obj, str):
        return obj
    else:
        # Fallback for other types: attempt to convert to string.
        # If it's a complex custom object, this might be a lossy conversion.
        logger.debug(f"Converting type {type(obj)} to string for serialization.")
        try:
            # For some specific complex types, a targeted conversion or JSON dump might be better
            # Example: if obj is some specific BioPython object, etc.
            # For now, a general attempt to JSON dump if it's not basic, then string
            if not isinstance(obj, (int, float, bool, str, bytes)) and obj is not None:
                 try:
                     return json.dumps(obj, default=str) # Try JSON first for unknown complex types
                 except TypeError:
                     return str(obj) # Final fallback to string
            return str(obj)
        except Exception as e_str:
            logger.error(f"Failed to convert object of type {type(obj)} to string: {e_str}. Returning empty string.")
            return ""


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
        import traceback
        logger.error(traceback.format_exc()) # Print full traceback for debugging
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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/create_combined_dataset_all_genes_sparse.py`

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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/entrez-to-ensembl-mapping.py`

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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/fix_placeholder_ids.py`

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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/gene_id_mapping_reference.py`

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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/generate_encode_mapping.py`

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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/preprocess_dataset_gene_ids.py`

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
    encode_mapping_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/encode_specific_mappings/encode_id_to_ensembl_mapping.csv'
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
    save_script = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/anndata_save_wrapper.py'
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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/rnaseq_utils.py`

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
import logging

# Set up logging
logger = logging.getLogger('rnaseq_utils')

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
GENCODE_MAPPING_FILE = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/gencode_v24_complete_mapping.csv")
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

NIH_RACE_CATEGORIES_LOWER = sorted(list(set([
    "american indian or alaska native", "asian", "black or african american",
    "native hawaiian or other pacific islander", "white", "more than one race",
    "multiethnic", "unknown or not reported", "other"
])))

NIH_ETHNICITY_CATEGORIES_LOWER = sorted(list(set([
    "hispanic or latino", "not hispanic or latino", "unknown or not reported"
])))

# Categories for the 'self_reported_ethnicity' string label column
# This can also be defined here, or constructed where needed if it combines the above.
# For consistency, let's define a base set.
SRE_BASE_CATEGORIES_LOWER = sorted(list(set(
    NIH_RACE_CATEGORIES_LOWER +
    ["hispanic or latino"] # Add if "hispanic or latino" can be a standalone primary category
    # Add other combined terms if you decide on them as standard string labels,
    # e.g., "white hispanic or latino"
)))

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
        'SPECIES_TO_NCBI_TAXON': load_json_mapping('species_to_taxon.json', {}),
        'ETHNICITY_TO_HANCESTRO': load_json_mapping('ethnicity_to_hancestro.json', {}) 
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

# In rnaseq_utils.py

def standardize_metadata(metadata_df, dataset_name, mappings=None):
    """
    Standardize metadata across datasets to ensure consistent fields and values.
    Maps fields to standard ontology terms where possible.
    
    Args:
        metadata_df: Metadata DataFrame
        dataset_name: Source dataset name
        mappings: Dictionary with mapping information
        
    Returns:
        Standardized metadata DataFrame with ontology mappings, or None if a critical error occurs.
    """
    try:
        if metadata_df is None:
            logger.error(f"standardize_metadata called with None input DataFrame for {dataset_name}. Returning None.")
            return None

        # Load mappings if not provided
        if mappings is None:
            logger.info("standardize_metadata: Mappings not provided, loading defaults.")
            mappings = load_mappings() # Assuming load_mappings() is defined in this file or imported
    
        # Make a copy to avoid modifying the original
        df = metadata_df.copy()
        
        # Get specific mappings - ensure these keys exist in `mappings` or handle missing keys
        tissue_to_uberon = mappings.get('TISSUE_TO_UBERON', {})
        assay_to_efo = mappings.get('ASSAY_TO_EFO', {})
        sex_standardization = mappings.get('SEX_STANDARDIZATION', {})
        species_to_taxon = mappings.get('SPECIES_TO_NCBI_TAXON', {})
        hancestro_map_loaded = mappings.get('ETHNICITY_TO_HANCESTRO', {}) # For HANCESTRO IDs

        # Define CORE_METADATA_FIELDS if not globally available in this context
        # These are fields that are fundamental and should exist.
        CORE_METADATA_FIELDS = [
            'sample_id', 'subject_id', 'sex', 'age', 'tissue', 
            'dataset', 'data_type', 'expression_unit', 'species'
            # Add race/ethnicity fields once they are consistently populated by upstream processes
        ]

        # Add dataset identifier
        df['dataset'] = dataset_name
        
        # Ensure all core fields exist (initialize with appropriate defaults if missing)
        for field in CORE_METADATA_FIELDS:
            if field not in df.columns:
                logger.debug(f"Standardize_metadata for {dataset_name}: Core field '{field}' missing. Adding with default.")
                if field in ['age', 'expression_unit', 'tissue']: # Fields that might be empty string
                    df[field] = "" 
                else: # Fields that usually have an 'unknown' category
                    df[field] = "unknown" 
        
        # Standardize species
        if 'species' in df.columns:
            df['species_ontology'] = df['species'].astype(str).str.lower().map(
                lambda x: species_to_taxon.get(x, species_to_taxon.get("unknown", ""))
            ).fillna("")
        else: # Should have been added by CORE_METADATA_FIELDS check
            df['species'] = "human" # Default
            df['species_ontology'] = species_to_taxon.get("human", "")

        # Standardize sex values
        if 'sex' in df.columns:
            # Ensure 'sex' column is string type before mapping for robustness
            df['sex'] = df['sex'].astype(str).str.lower().map(
                lambda x: sex_standardization.get(x, "unknown")
            ).fillna("unknown")
        else: # Should have been added by CORE_METADATA_FIELDS check
            df['sex'] = "unknown"


        # Standardize tissue to Uberon ontology
        if 'tissue' in df.columns:
            tissue_source_data = df['tissue']
            if isinstance(tissue_source_data, pd.DataFrame): # Defensive check
                logger.warning(
                    f"standardize_metadata for {dataset_name}: df['tissue'] is a DataFrame. Using its first column for 'tissue_original'. "
                    f"Columns found: {tissue_source_data.columns.tolist()}"
                )
                if not tissue_source_data.empty and len(tissue_source_data.columns) > 0:
                    df['tissue_original'] = tissue_source_data.iloc[:, 0].copy()
                    df['tissue'] = tissue_source_data.iloc[:, 0].copy() 
                else:
                    logger.error(f"standardize_metadata for {dataset_name}: df['tissue'] is an empty DataFrame or has no columns.")
                    df['tissue_original'] = "" 
                    df['tissue'] = ""          
            else:
                df['tissue_original'] = df['tissue'].copy() 

            if isinstance(df['tissue'], pd.DataFrame): # Final check
                logger.error(f"standardize_metadata for {dataset_name}: df['tissue'] is still a DataFrame. Aborting tissue processing.")
                df['tissue_ontology'] = "" 
            else:
                normalized_tissue = df['tissue'].astype(str).str.lower().str.strip()
                df['tissue'] = normalized_tissue # Update main 'tissue' column
                
                # Ensure tissue_to_uberon keys are lowercase for mapping
                tissue_to_uberon_lower = {k.lower(): v for k,v in tissue_to_uberon.items()}
                df['tissue_ontology'] = normalized_tissue.map(tissue_to_uberon_lower).fillna("")
                
                # Fallback substring matching for unmapped tissues (optional, can be noisy)
                # unmapped_mask = (df['tissue_ontology'] == "") & (df['tissue'] != "")
                # if unmapped_mask.any():
                #     # ... (substring matching logic, careful with this as it can mis-map) ...
                #     pass
        else: # if 'tissue' column doesn't exist, ensure ontology column does
             df['tissue_ontology'] = ""
             if 'tissue_original' not in df.columns: df['tissue_original'] = ""


        # Standardize data_type to EFO ontology
        if 'data_type' in df.columns:
            df['data_type_original'] = df['data_type'].copy()
            df['assay_ontology'] = df['data_type'].astype(str).str.lower().map(assay_to_efo).fillna("")
            
            if 'extraction_method' in df.columns:
                mask = (df['assay_ontology'] == "") & (df['extraction_method'].notna())
                df.loc[mask, 'assay_ontology'] = df.loc[mask, 'extraction_method'].astype(str).str.lower().map(assay_to_efo).fillna("")
        else: # Should have been added by CORE_METADATA_FIELDS
            df['data_type_original'] = ""
            df['assay_ontology'] = ""


        # Handle subject_id fallback
        if 'subject_id' not in df.columns or df['subject_id'].fillna('').eq('').all():
            if 'donor_id' in df.columns and not df['donor_id'].fillna('').eq('').all():
                logger.info(f"Using donor_id as subject_id for {dataset_name}")
                df['subject_id'] = df['donor_id']
            elif 'donor' in df.columns and not df['donor'].fillna('').eq('').all():
                logger.info(f"Using donor as subject_id for {dataset_name}")
                df['subject_id'] = df['donor']    
            else: # Ensure subject_id column exists even if empty
                df['subject_id'] = "unknown_subject_" + df.index.astype(str) # Fallback if truly missing
        
        # Log unmapped tissues (using map_tissue_to_ontology helper if it exists and is preferred)
        # For now, simple check based on tissue_ontology column
        if 'tissue' in df.columns and 'tissue_ontology' in df.columns:
            unmapped_tissues_final = df.loc[(df['tissue_ontology'] == "") & (df['tissue'] != ""), 'tissue'].dropna().unique()
            if len(unmapped_tissues_final) > 0:
                logger.warning(f"{dataset_name}: Unmapped tissues after processing: {', '.join(unmapped_tissues_final[:20])}")
                    

        # Standardize age and map to HsapDv
        if 'age' in df.columns:
            # --- BEGIN MODIFICATION ---
            age_source_data = df['age']
            if isinstance(age_source_data, pd.DataFrame): # Defensive check for age
                logger.warning(
                    f"standardize_metadata for {dataset_name}: df['age'] is a DataFrame. Using its first sub-column."
                )
                if not age_source_data.empty and len(age_source_data.columns) > 0:
                    # Ensure age_original captures the original problematic structure's first column if needed, or just the scalar
                    df['age_original'] = age_source_data.iloc[:, 0].copy().astype(str) # Keep original as series of strings
                    df['age'] = age_source_data.iloc[:, 0].copy().astype(str).replace('nan', '').replace('None', '').str.strip() # df['age'] becomes a Series
                else: # DataFrame column exists but is empty or has no sub-columns
                    logger.warning(f"standardize_metadata for {dataset_name}: df['age'] DataFrame is empty/malformed. Setting to empty string.")
                    df['age_original'] = ""
                    df['age'] = "" # df['age'] becomes a Series of empty strings
            else: # age column is already a Series
                df['age_original'] = df['age'].copy() # df['age'] is a Series
                df['age'] = df['age'].astype(str).replace('nan', '').replace('None', '').str.strip()
            # --- END MODIFICATION ---

            age_to_hsapdv_map = mappings.get('AGE_TO_HSAPDV', {})
            # df['age'] is now guaranteed to be a Series
            df['developmental_stage_ontology'] = df['age'].apply(
                lambda x: map_age_to_hsapdv(x, age_to_hsapdv_map) # map_age_to_hsapdv needs to be defined/imported
            ).fillna("")
        else: # age column doesn't exist
             logger.debug(f"Standardize_metadata: Adding missing 'age' column with default for {dataset_name}")
             df['age_original'] = ""
             df['age'] = "" # Ensure age column exists as string Series
             df['developmental_stage_ontology'] = ""            
            
        

        # --- Standardize Self-Reported Ethnicity and Race ---
        logger.debug(f"Standardizing ethnicity/race for {dataset_name}")


        # Ensure 'race' and 'is_hispanic_or_latino' columns exist, populated by dataset-specific logic or defaults here
        if 'race' not in df.columns: df['race'] = 'unknown or not reported'
        if 'is_hispanic_or_latino' not in df.columns: df['is_hispanic_or_latino'] = 'unknown or not reported'

        df['race_original'] = df['race'].astype(str).fillna('')
        df['race'] = df['race'].astype(str).str.lower().str.strip().fillna('unknown or not reported')
        df['race'] = df['race'].apply(lambda x: x if x in NIH_RACE_CATEGORIES_LOWER else 'unknown or not reported')
        
        df['ethnicity_original_is_hispanic'] = df['is_hispanic_or_latino'].astype(str).fillna('')
        df['is_hispanic_or_latino'] = df['is_hispanic_or_latino'].astype(str).str.lower().str.strip().fillna('unknown or not reported')
        df['is_hispanic_or_latino'] = df['is_hispanic_or_latino'].apply(lambda x: x if x in NIH_ETHNICITY_CATEGORIES_LOWER else 'unknown or not reported')

        sre_values = []
        for index, row in df.iterrows():
            race_val = str(row.get('race', 'unknown or not reported')) # Use .get for safety
            is_hispanic_val = str(row.get('is_hispanic_or_latino', 'unknown or not reported'))

            if race_val == "more than one race":
                sre_values.append("multiethnic")
            elif race_val == "unknown or not reported" and is_hispanic_val == "unknown or not reported":
                sre_values.append("unknown or not reported")
            elif is_hispanic_val == "hispanic or latino":
                if race_val != "unknown or not reported" and race_val != "multiethnic":
                    # Example: if race is "white", self_reported_ethnicity could be "white hispanic or latino"
                    # Or, if your HANCESTRO map has a term for "White" and another for "Hispanic or Latino",
                    # this combination might be mapped to "multiethnic" for the ontology ID.
                    # For now, let's create a combined label or default to a broader term.
                    # This decision depends on your HANCESTRO mapping strategy.
                    # A simple approach for now:
                    sre_values.append(f"{race_val} (hispanic or latino)") # Or simply "hispanic or latino" or race_val or "multiethnic"
                else: # race is unknown or multiethnic, but is Hispanic
                    sre_values.append("hispanic or latino")
            else: # Not Hispanic or Latino, or ethnicity unknown. Default to race value.
                sre_values.append(race_val)
        
        # Define categories for self_reported_ethnicity based on potential values
        sre_categories = sorted(list(set(NIH_RACE_CATEGORIES_LOWER + \
                                         [f"{r} (hispanic or latino)" for r in NIH_RACE_CATEGORIES_LOWER if r not in ["unknown or not reported", "multiethnic"]] + \
                                         ["hispanic or latino"])))

        df['self_reported_ethnicity'] = pd.Categorical(
            sre_values,
            categories=sre_categories,
            ordered=False
        ).fillna("unknown or not reported")

        hancestro_map_from_json = mappings.get('ETHNICITY_TO_HANCESTRO', {}) # Renamed for clarity
        hancestro_map_lower = {k.lower(): v for k, v in hancestro_map_from_json.items()}
        
        df['self_reported_ethnicity_ontology_term_id'] = df['self_reported_ethnicity'].astype(str).str.lower().map(hancestro_map_lower).fillna('unknown')

        # Validate against schema-allowed values for the ontology ID
        # The values in hancestro_map_lower.values() are the HANCESTRO IDs
        temp_allowed_ids = list(hancestro_map_lower.values())
        temp_allowed_ids.append("multiethnic")
        temp_allowed_ids.append("unknown")
        schema_allowed_ontology_ids = sorted(list(set(temp_allowed_ids))) # <<< ENSURE UNIQUE AND SORTED
        
        invalid_ontology_mask = ~df['self_reported_ethnicity_ontology_term_id'].isin(schema_allowed_ontology_ids)
        if invalid_ontology_mask.any():
            # ... (warning log) ...
            df.loc[invalid_ontology_mask, 'self_reported_ethnicity_ontology_term_id'] = 'unknown'
        
        # Ensure the final column is categorical with these allowed_ontology_ids
        df['self_reported_ethnicity_ontology_term_id'] = pd.Categorical(
            df['self_reported_ethnicity_ontology_term_id'],
            categories=schema_allowed_ontology_ids, # Use the unique, sorted list
            ordered=False
        )
        
        # --- Finalize Categorical Dtypes for all relevant fields ---
        # Ensure all fields intended to be categorical are set as such, with consistent categories where possible.
        # This list should include all fields that are meant to be categorical in the final obs.

        final_categorical_fields = [
            'sex', 'tissue', 'dataset', 'species', 'data_type', 
            'race', 'is_hispanic_or_latino', 'self_reported_ethnicity', 
            'self_reported_ethnicity_ontology_term_id',
            'tissue_ontology', 'assay_ontology', 'species_ontology', 'developmental_stage_ontology'
        ]
        if 'cell_type' in df.columns: final_categorical_fields.append('cell_type')
        if 'disease' in df.columns: final_categorical_fields.append('disease')

        for field in final_categorical_fields:
            if field in df.columns:
                current_categories_for_field = None # Use specific lists if available
                default_val_for_cat = 'unknown'
                
                if field == 'race': 
                    current_categories_for_field = NIH_RACE_CATEGORIES_LOWER # Uses module global
                    default_val_for_cat = 'unknown or not reported'
                elif field == 'is_hispanic_or_latino': 
                    current_categories_for_field = NIH_ETHNICITY_CATEGORIES_LOWER # Uses module global
                    default_val_for_cat = 'unknown or not reported'
                elif field == 'self_reported_ethnicity':
                    current_categories_for_field = SRE_BASE_CATEGORIES_LOWER # Uses module global
                    default_val_for_cat = 'unknown or not reported'
                elif field == 'self_reported_ethnicity_ontology_term_id':
                    # This list is dynamically built based on hancestro_map earlier in this function
                    current_categories_for_field = schema_allowed_ontology_ids 
                
                try:
                    # Ensure all values are strings before converting to category, handle NaNs
                    str_series = df[field].astype(str).replace('nan', default_val_for_cat).fillna(default_val_for_cat)
                    if current_categories_for_field:
                        # Ensure all values in series are part of defined categories
                        # And that current_categories_for_field itself has unique values
                        unique_cats = sorted(list(set(current_categories_for_field)))
                        valid_series = str_series.apply(lambda x: x if x in unique_cats else default_val_for_cat)
                        df[field] = pd.Categorical(valid_series, categories=unique_cats, ordered=False)
                    else: # Let pandas infer categories for other fields like tissue_ontology, etc.
                        df[field] = str_series.astype('category')
                except Exception as e_cat_final:
                    logger.error(f"Error converting field '{field}' to category for {dataset_name}: {e_cat_final}. Leaving as string.")
                    df[field] = df[field].astype(str).fillna(default_val_for_cat)
            else:
                logger.debug(f"Standardize_metadata: Adding missing categorical field '{field}' as 'unknown' for {dataset_name}.")
                df[field] = pd.Categorical(['unknown'] * len(df), categories=['unknown'])
        return df
    except Exception as e_std_meta:
        logger.error(f"CRITICAL ERROR within standardize_metadata for {dataset_name}: {e_std_meta}")
        import traceback
        logger.error(traceback.format_exc())
        return None # Ensure None is returned on critical failure

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
            adata.uns['metadata_sources'] = [str(adata.uns['metadata_sources'])]

        source_info = metadata['metadata_source'] # Get the actual value
        
        # Ensure all elements in source_info, if it's a list, are serializable
        # This is more for the 'extras' field which is a list of dicts.
        # For metadata_source, it's usually a single dict.
        
        if isinstance(source_info, dict):
            try:
                # Serialize the source_info dict itself
                serializable_source_info = ensure_serializable(source_info) # Use ensure_serializable from standardize_datasets.py
                source_str = json.dumps(serializable_source_info, default=str)
                adata.uns['metadata_sources'].append(source_str)
                logger.debug(f"Appended metadata_source as JSON string: {source_str[:100]}...")
            except TypeError as te: # Handles cases where ensure_serializable might not fully clean it for json.dumps
                logger.warning(f"Could not directly JSON dump complex metadata_source: {te}. Converting to string.")
                adata.uns['metadata_sources'].append(str(source_info))
            except Exception as e:
                logger.error(f"Could not serialize metadata_source dict to JSON: {e}. Appending as raw string.")
                adata.uns['metadata_sources'].append(str(source_info))
        elif isinstance(source_info, str): # If it's already a string
            adata.uns['metadata_sources'].append(source_info)
        else: # Other types, convert to string
            adata.uns['metadata_sources'].append(str(source_info))    
    

    logger.info("Dataset-specific metadata applied successfully")
    return adata

````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/run_rnaseq_pipeline.sh`

````
#!/bin/bash
# Run the complete RNA-seq standardization pipeline with improved gene ID mapping
# Uses separate base directories for inputs and outputs.

# --- Configuration ---
INPUT_BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
OUTPUT_BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq"

SCRIPTS_DIR="${INPUT_BASE_DIR}/scripts/pipeline_v2" # Scripts are read from input location

# --- INPUT DATA PATHS (relative to INPUT_BASE_DIR or absolute) ---
# Metadata JSON files (configuration for the pipeline itself)
PIPELINE_METADATA_JSON_DIR="${INPUT_BASE_DIR}/metadata/json" # General pipeline config JSONs
# Gene mapping resources (source files like GTF, NCBI downloaded maps)
SOURCE_GENE_MAPPING_RESOURCE_DIR="${INPUT_BASE_DIR}/metadata/gene_mapping"
DOWNLOADED_GTF_GZ="${SOURCE_GENE_MAPPING_RESOURCE_DIR}/gencode.v24.annotation.gtf.gz"
UNZIPPED_GTF_FILE_SOURCE="${SOURCE_GENE_MAPPING_RESOURCE_DIR}/gencode.v24.annotation.gtf" # Source if already unzipped

# Raw data locations
ENCODE_RAW_DATA_DIR="${INPUT_BASE_DIR}/encode/raw_data"
ENTEX_RAW_DATA_DIR="${INPUT_BASE_DIR}/encode/entex"
ENTEX_METADATA_FILE_INPUT="${INPUT_BASE_DIR}/encode/metadata/entex_metadata.json"
MAGE_RAW_DATA_DIR="${INPUT_BASE_DIR}/mage"
ADNI_RAW_DATA_DIR="${INPUT_BASE_DIR}/adni_microarray"
ADNI_DEMOGRAPHICS_FILE_INPUT="${INPUT_BASE_DIR}/metadataADNI/subject_demographics/PTDEMOG_25Apr2025.csv"
ADNI_DIAGNOSIS_FILE_INPUT="${INPUT_BASE_DIR}/metadataADNI/subject_demographics/DXSUM_30Apr2025.csv"
GTEX_RAW_FILE_INPUT="${INPUT_BASE_DIR}/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"
MAGE_1000G_PED_FILE_INPUT="${INPUT_BASE_DIR}/project_gene_regulation/data/MAGE/WGS/20130606_g1k_3202_samples_ped_population.txt"


# --- OUTPUT & GENERATED DATA PATHS (relative to OUTPUT_BASE_DIR) ---
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
# Main output directories for H5ADs
OUTPUT_DIR_H5ADS="${OUTPUT_BASE_DIR}/standardized_data/run_${TIMESTAMP}"
PREPROCESSED_DIR_H5ADS="${OUTPUT_BASE_DIR}/preprocessed_data/run_${TIMESTAMP}"
# Log directory
LOG_DIR_PIPELINE="${OUTPUT_BASE_DIR}/logs/pipeline_runs"
# Directory for *generated* metadata and gene mapping files (e.g., NCBI map, primary ref map)
GENERATED_METADATA_JSON_DIR="${OUTPUT_BASE_DIR}/metadata/json_generated"
GENERATED_GENE_MAPPING_RESOURCE_DIR="${OUTPUT_BASE_DIR}/metadata/gene_mapping_generated"

# Paths for generated mapping files
ENTREZ_TO_ENSEMBL_CSV_GENERATED="${GENERATED_METADATA_JSON_DIR}/entrez_to_ensembl_mapping.csv"
PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED="${GENERATED_GENE_MAPPING_RESOURCE_DIR}/gencode_v24_complete_mapping.csv"
UNZIPPED_GTF_FILE_OUTPUT="${GENERATED_GENE_MAPPING_RESOURCE_DIR}/gencode.v24.annotation.gtf" # Where unzipped version will be placed if re-generated
ENCODE_ID_MAP_OUTPUT_DIR_GENERATED="${GENERATED_GENE_MAPPING_RESOURCE_DIR}/encode_specific_mappings"


# --- Logging Setup ---
mkdir -p "$LOG_DIR_PIPELINE"
PIPELINE_LOG_FILE="${LOG_DIR_PIPELINE}/pipeline_run_${TIMESTAMP}.log"

# --- Helper Functions ---
log_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$PIPELINE_LOG_FILE"
}

run_command() {
    log_message "Running: $1"
    PYTHONPATH="${SCRIPTS_DIR}:${PYTHONPATH}" eval "$1" 2>&1 | tee -a "$PIPELINE_LOG_FILE"
    return ${PIPESTATUS[0]}
}

# --- Pipeline Start ---
log_message "=== Starting RNA-seq Pipeline Run: ${TIMESTAMP} ==="
log_message "INPUT BASE DIRECTORY: ${INPUT_BASE_DIR}"
log_message "OUTPUT BASE DIRECTORY: ${OUTPUT_BASE_DIR}"

# Create output structure
mkdir -p "$OUTPUT_DIR_H5ADS" "$PREPROCESSED_DIR_H5ADS" \
           "$GENERATED_METADATA_JSON_DIR" "$GENERATED_GENE_MAPPING_RESOURCE_DIR" \
           "${OUTPUT_DIR_H5ADS}/temp" # Temp dir for standardize_datasets.py

# --- Argument Parsing for --force ---
FORCE_FLAG=""
EFFECTIVE_FORCE_FOR_REF_MAP=""
if [ "$1" == "--force" ] || [ "$1" == "--force-all" ]; then
    FORCE_FLAG="--force"
    EFFECTIVE_FORCE_FOR_REF_MAP="--force"
    log_message "Force flag ('$1') detected. Applicable scripts will regenerate outputs."
elif [ "$1" == "--force-mapping" ]; then
    EFFECTIVE_FORCE_FOR_REF_MAP="--force"
    log_message "Force mapping flag detected. Gene mapping files will be regenerated."
fi

# --- PRE-RUN CHECKS FOR INPUTS ---
log_message "--- Performing Pre-run Input Checks (from INPUT_BASE_DIR) ---"
if [ ! -d "$SCRIPTS_DIR" ]; then log_message "ERROR: Scripts directory $SCRIPTS_DIR not found! Exiting."; exit 1; fi
if [ ! -f "$DOWNLOADED_GTF_GZ" ] && [ ! -f "$UNZIPPED_GTF_FILE_SOURCE" ]; then log_message "WARNING: GENCODE GTF (gzipped or source unzipped) not found at expected input path. Download or ensure correct path."; fi
if [ -n "$MAGE_RAW_DATA_DIR" ] && [ ! -d "$MAGE_RAW_DATA_DIR" ]; then log_message "WARNING: MAGE raw data directory $MAGE_RAW_DATA_DIR not found."; fi
if [ ! -f "$GTEX_RAW_FILE_INPUT" ]; then log_message "WARNING: GTEx raw file $GTEX_RAW_FILE_INPUT not found."; fi
# Add more checks for other critical inputs if desired
log_message "--- Pre-run Input Checks Complete ---"


# --- Stage 0: Prepare Entrez to Ensembl Mapping (NCBI Source) ---
# This will be written to the OUTPUT_BASE_DIR structure
log_message "--- Stage 0: Preparing Entrez to Ensembl Mapping ---"
mkdir -p "$(dirname "${ENTREZ_TO_ENSEMBL_CSV_GENERATED}")"
if [ "$EFFECTIVE_FORCE_FOR_REF_MAP" == "--force" ] || [ ! -f "$ENTREZ_TO_ENSEMBL_CSV_GENERATED" ]; then
    run_command "python \"${SCRIPTS_DIR}/entrez-to-ensembl-mapping.py\" --output \"${ENTREZ_TO_ENSEMBL_CSV_GENERATED}\" --species human --force"
    if [ $? -ne 0 ]; then log_message "ERROR: Entrez to Ensembl mapping generation failed. Exiting."; exit 1; fi
else
    log_message "Entrez to Ensembl mapping file ${ENTREZ_TO_ENSEMBL_CSV_GENERATED} already exists and no relevant force flag. Skipping generation."
fi

# --- Stage 0.5: Unzip GENCODE GTF ---
# Unzips from INPUT_BASE_DIR to OUTPUT_BASE_DIR if needed
log_message "--- Stage 0.5: Preparing GENCODE v24 GTF File ---"
mkdir -p "$(dirname "${UNZIPPED_GTF_FILE_OUTPUT}")"
if [ ! -f "$DOWNLOADED_GTF_GZ" ]; then
    log_message "ERROR: Gzipped GTF file $DOWNLOADED_GTF_GZ not found in input path."
    if [ ! -f "$UNZIPPED_GTF_FILE_OUTPUT" ] && [ ! -f "$UNZIPPED_GTF_FILE_SOURCE" ]; then
        log_message "ERROR: And no unzipped GTF file found in output or source paths. Exiting."; exit 1;
    elif [ -f "$UNZIPPED_GTF_FILE_SOURCE" ] && [ ! -f "$UNZIPPED_GTF_FILE_OUTPUT" ] ; then
        log_message "WARNING: Gzipped GTF $DOWNLOADED_GTF_GZ not found. Using source unzipped version $UNZIPPED_GTF_FILE_SOURCE and copying to output."
        cp "$UNZIPPED_GTF_FILE_SOURCE" "$UNZIPPED_GTF_FILE_OUTPUT"
    elif [ -f "$UNZIPPED_GTF_FILE_OUTPUT" ]; then
        log_message "WARNING: Gzipped GTF $DOWNLOADED_GTF_GZ not found, but unzipped version $UNZIPPED_GTF_FILE_OUTPUT exists in output. Proceeding."
    fi
elif [ "$EFFECTIVE_FORCE_FOR_REF_MAP" == "--force" ] || [ ! -f "$UNZIPPED_GTF_FILE_OUTPUT" ] || [ "$DOWNLOADED_GTF_GZ" -nt "$UNZIPPED_GTF_FILE_OUTPUT" ]; then
    log_message "Unzipping $DOWNLOADED_GTF_GZ to $UNZIPPED_GTF_FILE_OUTPUT..."
    gunzip -c "$DOWNLOADED_GTF_GZ" > "$UNZIPPED_GTF_FILE_OUTPUT"; if [ $? -ne 0 ]; then log_message "ERROR: Failed to unzip GTF. Exiting."; exit 1; fi
else
    log_message "Unzipped GTF ${UNZIPPED_GTF_FILE_OUTPUT} is up-to-date in output. Skipping unzip."
fi
# Ensure the GTF path used by subsequent scripts is the one in the output/generated directory
FINAL_UNZIPPED_GTF_TO_USE=$UNZIPPED_GTF_FILE_OUTPUT
if [ ! -f "$FINAL_UNZIPPED_GTF_TO_USE" ] && [ -f "$UNZIPPED_GTF_FILE_SOURCE" ]; then
    log_message "Unzipped GTF not in output, using source: $UNZIPPED_GTF_FILE_SOURCE"
    FINAL_UNZIPPED_GTF_TO_USE=$UNZIPPED_GTF_FILE_SOURCE
elif [ ! -f "$FINAL_UNZIPPED_GTF_TO_USE" ]; then
    log_message "ERROR: No unzipped GTF file available. Exiting."
    exit 1
fi


# --- Stage A: Generate THE Primary Gene Annotation/Reference Map ---
# This map will be written to OUTPUT_BASE_DIR structure
log_message "--- Stage A: Generating Primary Gene Annotation & Reference Map ---"
TEMP_DOWNLOAD_DIR_FOR_REF="${OUTPUT_BASE_DIR}/temp_gene_ref_downloads" # Temp dir under output
mkdir -p "$TEMP_DOWNLOAD_DIR_FOR_REF"
mkdir -p "$(dirname "${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}")"

run_command "python \"${SCRIPTS_DIR}/gene_id_mapping_reference.py\" \\
    --encode-dir \"${ENCODE_RAW_DATA_DIR}\" \\
    --entex-dir \"${ENTEX_RAW_DATA_DIR}\" \\
    --gencode-gtf \"${FINAL_UNZIPPED_GTF_TO_USE}\" \\
    --entrez-mapping \"${ENTREZ_TO_ENSEMBL_CSV_GENERATED}\" \\
    --output \"${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}\" \\
    --temp-dir \"${TEMP_DOWNLOAD_DIR_FOR_REF}\" \\
    ${EFFECTIVE_FORCE_FOR_REF_MAP}"
if [ $? -ne 0 ]; then log_message "ERROR: Stage A (Primary Gene Annotation Map generation to ${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}) failed. Exiting."; exit 1; fi
log_message "Stage A completed. Primary gene annotation & reference map: ${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}"


# --- Stage 1: Initial Data Conversion (Raw to Standardized v1 H5ADs) ---
# Reads from INPUT_BASE_DIR, Writes to OUTPUT_DIR_H5ADS
# standardize_datasets.py will need to be aware of PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED
# This means rnaseq_utils.GENCODE_MAPPING_FILE might need to be overridden or the script adapted
# For now, we assume standardize_datasets.py is modified to accept this path or uses a symlink trick.
# A simpler way is to make sure rnaseq_utils.py points to the generated one.
# Since rnaseq_utils.py is in INPUT_BASE_DIR, it might be tricky.
# Let's assume for now that `standardize_datasets.py` is modified to take the primary map path as an arg, or we set an env var.
# For this iteration, we'll ensure the `rnaseq_utils.py` used by `standardize_datasets.py` is configured to look for
# the *generated* primary mapping file. This is usually handled by GENCODE_MAPPING_FILE in rnaseq_utils.py.
# A robust way is to pass the path to the standardize_datasets.py script.
# *Self-correction: standardize_datasets.py via rnaseq_utils.load_gencode_mapping() uses a fixed path.
# The simplest is to ensure that fixed path (GENCODE_MAPPING_FILE in rnaseq_utils.py) points to *our generated one*.
# This means rnaseq_utils.py needs to be aware of the generated path.
# Let's make rnaseq_utils.py flexible or ensure the path it looks for IS the generated one.
# For now, let's copy the generated primary map to the path expected by rnaseq_utils.py if it's different.
# The `rnaseq_utils.py` expects it at `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/gencode_v24_complete_mapping.csv`
# If INPUT_BASE_DIR is the same, it works. If different, we need to handle this.
# For this separated path logic, it's best if `standardize_datasets.py` can take the gene map path as an argument.
# Assuming `standardize_datasets.py` (and by extension `rnaseq_utils.py`) will correctly use the
# `PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED` (this requires modification to those scripts to accept it as param, or env var).
# For now, the script has rnaseq_utils.GENCODE_MAPPING_FILE hardcoded.
# The most direct fix is to ensure the `run_command` for standardize_datasets.py uses a python path that sees an rnaseq_utils.py
# configured to point to PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED. This is complex.
# Alternative: symlink or copy PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED to the location expected by the *original* rnaseq_utils.py
EXPECTED_GENCODE_MAP_PATH_FOR_UTILS="${INPUT_BASE_DIR}/metadata/gene_mapping/gencode_v24_complete_mapping.csv"
if [ "$PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED" != "$EXPECTED_GENCODE_MAP_PATH_FOR_UTILS" ]; then
    log_message "Copying generated primary gene map ${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED} to expected location ${EXPECTED_GENCODE_MAP_PATH_FOR_UTILS} for Stage 1..."
    mkdir -p "$(dirname "$EXPECTED_GENCODE_MAP_PATH_FOR_UTILS")"
    cp "$PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED" "$EXPECTED_GENCODE_MAP_PATH_FOR_UTILS"
    # Also copy the JSON version if it exists
    GENERATED_JSON_MAP="${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED%.csv}.json"
    EXPECTED_JSON_MAP="${EXPECTED_GENCODE_MAP_PATH_FOR_UTILS%.csv}.json"
    if [ -f "$GENERATED_JSON_MAP" ]; then
      cp "$GENERATED_JSON_MAP" "$EXPECTED_JSON_MAP"
    fi
fi

log_message "--- Stage 1: Initial Data Conversion (Raw to Standardized v1 H5ADs) ---"
ADNI_DEMOGRAPHICS_ARG=""; if [ -f "$ADNI_DEMOGRAPHICS_FILE_INPUT" ]; then ADNI_DEMOGRAPHICS_ARG="--adni-demographics-file \"${ADNI_DEMOGRAPHICS_FILE_INPUT}\""; else log_message "WARNING: ADNI demographics file ${ADNI_DEMOGRAPHICS_FILE_INPUT} missing."; fi
ADNI_DIAGNOSIS_ARG=""; if [ -f "$ADNI_DIAGNOSIS_FILE_INPUT" ]; then ADNI_DIAGNOSIS_ARG="--adni-diagnosis-file \"${ADNI_DIAGNOSIS_FILE_INPUT}\""; else log_message "WARNING: ADNI diagnosis file ${ADNI_DIAGNOSIS_FILE_INPUT} missing."; fi

run_command "python \"${SCRIPTS_DIR}/standardize_datasets.py\" \\
    --encode-dir \"${ENCODE_RAW_DATA_DIR}/cell_lines\" \\
    --encode-entex-dir \"${ENTEX_RAW_DATA_DIR}\" \\
    --entex-metadata-file \"${ENTEX_METADATA_FILE_INPUT}\" \\
    --mage-dir \"${MAGE_RAW_DATA_DIR}\" \\
    --adni-dir \"${ADNI_RAW_DATA_DIR}\" \\
    ${ADNI_DEMOGRAPHICS_ARG} \\
    ${ADNI_DIAGNOSIS_ARG} \\
    --gtex-file \"${GTEX_RAW_FILE_INPUT}\" \\
    --metadata-dir \"${PIPELINE_METADATA_JSON_DIR}\" \\
    --output-dir \"${OUTPUT_DIR_H5ADS}\""
if [ $? -ne 0 ]; then log_message "ERROR: Stage 1 (Initial Data Conversion) failed. Exiting."; exit 1; fi
log_message "Stage 1 completed."


# --- Stage 1.6: Generate ENCODE-specific ID to Ensembl Mapping ---
# Output to GENERATED_GENE_MAPPING_RESOURCE_DIR
log_message "--- Stage 1.6: Generating ENCODE ID to Ensembl Mapping ---"
mkdir -p "$ENCODE_ID_MAP_OUTPUT_DIR_GENERATED"
run_command "python \"${SCRIPTS_DIR}/generate_encode_mapping.py\" \\
    --encode-dir \"${ENCODE_RAW_DATA_DIR}\" \\
    --output-dir \"${ENCODE_ID_MAP_OUTPUT_DIR_GENERATED}\" \\
    ${FORCE_FLAG}"
if [ $? -ne 0 ]; then log_message "WARNING: Stage 1.6 (ENCODE ID Mapping) failed."; else log_message "Stage 1.6 completed."; fi


# --- Stage 2: Enhanced Metadata Standardization (v1 to v2 H5ADs) ---
# Reads from OUTPUT_DIR_H5ADS, writes to OUTPUT_DIR_H5ADS
log_message "--- Stage 2: Enhanced Metadata Standardization ---"
run_command "python \"${SCRIPTS_DIR}/standardize_metadata.py\" \\
    --data-dir \"${OUTPUT_DIR_H5ADS}\" \\
    --output-dir \"${OUTPUT_DIR_H5ADS}\" \\
    --metadata-dir \"${PIPELINE_METADATA_JSON_DIR}\"" # Uses config JSONs from input dir
if [ $? -ne 0 ]; then log_message "ERROR: Stage 2 (Enhanced Metadata Standardization) failed. Exiting."; exit 1; fi
log_message "Stage 2 completed."


# --- Stage 2.5: Preprocess Datasets (Gene ID Standardization for each H5AD) ---
# Reads from OUTPUT_DIR_H5ADS, uses generated primary map, writes to PREPROCESSED_DIR_H5ADS
log_message "--- Stage 2.5: Preprocessing Datasets for Consistent Gene IDs ---"
run_command "python \"${SCRIPTS_DIR}/preprocess_dataset_gene_ids.py\" \\
    --data-dir \"${OUTPUT_DIR_H5ADS}\" \\
    --reference-mapping \"${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}\" \\
    --output-dir \"${PREPROCESSED_DIR_H5ADS}\" \\
    ${FORCE_FLAG}"
if [ $? -ne 0 ]; then log_message "ERROR: Stage 2.5 (Dataset Gene ID Preprocessing) failed. Exiting."; exit 1; fi
log_message "Stage 2.5 completed. Preprocessed files in: $PREPROCESSED_DIR_H5ADS"


# --- Stage 2.6: Fix Placeholder Gene IDs ---
log_message "--- Stage 2.6: Fixing Placeholder Gene IDs in Preprocessed Datasets ---"
for dataset_label in encode gtex mage adni entex; do
    preprocessed_h5ad_file="${PREPROCESSED_DIR_H5ADS}/${dataset_label}_standardized_preprocessed.h5ad"
    if [ -f "$preprocessed_h5ad_file" ]; then
        log_message "Attempting to fix placeholder IDs in: ${preprocessed_h5ad_file}"
        run_command "python \"${SCRIPTS_DIR}/fix_placeholder_ids.py\" \"$preprocessed_h5ad_file\" \"$preprocessed_h5ad_file.fixed\""
        if [ $? -eq 0 ]; then
            mv "$preprocessed_h5ad_file.fixed" "$preprocessed_h5ad_file"
            log_message "Placeholder IDs fixed for ${dataset_label} dataset."
        else
            log_message "WARNING: Failed to fix placeholder IDs for ${dataset_label} dataset."
        fi
    else
        log_message "INFO: No preprocessed file found for ${dataset_label} at ${preprocessed_h5ad_file}, skipping placeholder fix."
    fi
done
log_message "Stage 2.6 completed."


# --- Stage 2.7: Analyze ENCODE Gene Mapping Quality ---
log_message "--- Stage 2.7: Analyzing ENCODE Gene Mapping Quality ---"
ENCODE_PREPROCESSED_H5AD="${PREPROCESSED_DIR_H5ADS}/encode_standardized_preprocessed.h5ad"
ENCODE_MAPPING_STATS_JSON="${OUTPUT_DIR_H5ADS}/encode_mapping_stats.json" # Output to main H5AD output dir

if [ -f "$ENCODE_PREPROCESSED_H5AD" ]; then
    run_command "python - <<EOF
import scanpy as sc
import pandas as pd
import numpy as np
import json

encode_path = '$ENCODE_PREPROCESSED_H5AD'
try:
    adata = sc.read_h5ad(encode_path)
    print(f\"ENCODE dataset shape: {adata.shape}\")

    non_empty = 0
    percentage = 0.0
    source_counts_dict = {}

    if 'ensembl_id' in adata.var.columns:
        non_empty = sum(1 for x in adata.var['ensembl_id'] if x and str(x).strip() != '' and not str(x).startswith('PLACEHOLDER_'))
        if len(adata.var) > 0:
            percentage = non_empty / len(adata.var) * 100
        print(f\"ENCODE genes with valid (non-placeholder) mapped Ensembl IDs: {non_empty}/{len(adata.var)} ({percentage:.2f}%)\")

        validly_mapped_ids = adata.var.loc[(adata.var['ensembl_id'].notna()) & (adata.var['ensembl_id'] != '') & (~adata.var['ensembl_id'].astype(str).str.startswith('PLACEHOLDER_')), 'ensembl_id']
        print(f\"Sample valid mapped Ensembl IDs: {validly_mapped_ids.unique()[:5].tolist()}\")

        if 'mapping_source' in adata.var.columns:
            source_counts = adata.var['mapping_source'].value_counts()
            source_counts_dict = {k: int(v) for k, v in source_counts.items()}
            for source, count_val in source_counts.items():
                source_percentage_val = 0.0 
                if len(adata.var) > 0:
                    source_percentage_val = count_val / len(adata.var) * 100
                print(f\"Mapping source '{source}': {count_val} ({source_percentage_val:.2f}%)\")
    else:
        print(\"Error: ENCODE preprocessed dataset does not have an 'ensembl_id' column in .var!\")

    mapping_stats = {
        'total_genes': len(adata.var) if hasattr(adata, 'var') else 0,
        'mapped_genes_valid_ensembl': non_empty,
        'mapping_percentage_valid_ensembl': percentage,
        'mapping_sources': source_counts_dict
    }

    with open('$ENCODE_MAPPING_STATS_JSON', 'w') as f:
        json.dump(mapping_stats, f, indent=2)

    print(f\"ENCODE mapping stats saved to $ENCODE_MAPPING_STATS_JSON\")

except FileNotFoundError:
    print(f\"Error: ENCODE preprocessed file not found at {encode_path}. Cannot analyze mapping quality.\")
except Exception as e:
    print(f\"Error analyzing ENCODE dataset: {e}\")
    print(\"Continuing with pipeline execution...\")
EOF"
else
    log_message "INFO: ENCODE preprocessed file not found at $ENCODE_PREPROCESSED_H5AD. Skipping ENCODE mapping quality analysis."
fi
log_message "Stage 2.7 completed."


# --- Final Message ---
log_message "=== RNA-seq Pipeline Main Processing (up to Stage 2.7) Finished ==="
log_message "Main Output Directory (v1, v2 H5ADs): $OUTPUT_DIR_H5ADS"
log_message "Preprocessed Data Directory (gene ID standardized H5ADs): $PREPROCESSED_DIR_H5ADS"
log_message "Generated Gene Mapping Files in: ${GENERATED_GENE_MAPPING_RESOURCE_DIR}"
log_message "Primary Gene Annotation & Reference Map used (and generated if forced): ${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}"
log_message "Full Log: $PIPELINE_LOG_FILE"
log_message "Consider running Stage 3 (Combined Dataset) and Stage 4 (Validation) next if desired."
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/setup_workspace.sh`

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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/standardize_datasets.py`

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
    ensure_serializable,
    standardize_ensembl_id,
    load_gencode_mapping, # Still import here, but call once in main
    add_gencode_annotations,
    standardize_metadata,
    validate_metadata,
    prepare_for_anndata,
    load_mappings,
    load_json_mapping,
    load_dataset_specific_metadata,
    apply_dataset_specific_metadata,
    NIH_RACE_CATEGORIES_LOWER,
    NIH_ETHNICITY_CATEGORIES_LOWER,
    SRE_BASE_CATEGORIES_LOWER
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


def identify_columns(df):
    """Identify gene_id and expression value columns across different file formats"""
    gene_id_col = None
    expr_val_col = None

    columns_lower = [str(col).lower() for col in df.columns]

    # Priority for gene identifier column
    gene_id_candidates = ["gene_id", "target_id", "gene", "name", "id"] # Ordered by preference
    for candidate in gene_id_candidates:
        if candidate in columns_lower:
            gene_id_col = df.columns[columns_lower.index(candidate)]
            logger.debug(f"Identified gene ID column: '{gene_id_col}' (matched '{candidate}')")
            break

    if gene_id_col is None and len(df.columns) >=1: # Fallback if no specific name matches
        gene_id_col = df.columns[0]
        logger.debug(f"Using fallback first column '{gene_id_col}' as gene ID column.")


    # Priority for expression value column
    expr_val_candidates = ["tpm", "fpkm", "counts", "normalized_intensity", "expression", "value"]
    for candidate in expr_val_candidates:
        if candidate in columns_lower:
            expr_val_col = df.columns[columns_lower.index(candidate)]
            logger.debug(f"Identified expression value column: '{expr_val_col}' (matched '{candidate}')")
            break

    if expr_val_col is None and len(df.columns) >= 2: # Fallback if no specific name matches
        # Try to find a numeric column if gene_id_col is the first one
        if gene_id_col == df.columns[0] and pd.api.types.is_numeric_dtype(df[df.columns[1]]):
            expr_val_col = df.columns[1]
            logger.debug(f"Using fallback second column '{expr_val_col}' as expression value column.")
        elif len(df.select_dtypes(include=np.number).columns) > 0 : # Generic numeric column
            expr_val_col = df.select_dtypes(include=np.number).columns[0]
            logger.debug(f"Using first numeric column '{expr_val_col}' as expression value column.")


    return gene_id_col, expr_val_col

def get_encode_cell_info(metadata_dir=None):
    """Get ENCODE cell line info from JSON."""
    if not metadata_dir:
        metadata_dir = DEFAULT_METADATA_DIR

    json_path = Path(metadata_dir) / "encode_metadata.json"
    if json_path.exists():
        try:
            with open(json_path, "r") as f:
                metadata_config = json.load(f)
            if "dataset_info" in metadata_config and "cell_lines" in metadata_config["dataset_info"]:
                logger.info(f"Successfully loaded ENCODE cell line info from {json_path}")
                return metadata_config["dataset_info"]["cell_lines"]
            else:
                logger.warning(f"'dataset_info.cell_lines' not found in {json_path}. No ENCODE cell line specific info will be used beyond defaults.")
        except Exception as e:
            logger.error(f"Error loading ENCODE cell info from {json_path}: {e}")
    else:
        logger.error(f"ENCODE metadata file not found: {json_path}. Cannot retrieve cell line info.")
    return {}




def create_standard_anndata(data_df, obs_df, var_df, dataset_specific_cfg, dataset_name_for_log):
    """
    Create standardized AnnData object with consistent structure.
    Relies on dataset_specific_cfg for .uns population.
    """
    logger.info(f"Creating AnnData for {dataset_name_for_log}")

    obs_df, var_df, data_df = prepare_for_anndata(obs_df.copy(), var_df.copy(), data_df.copy())

    logger.debug(
        f"Creating AnnData with X shape={data_df.shape}, obs shape={obs_df.shape}, var shape={var_df.shape}"
    )
    try:
        # Try sparse matrix creation first for X, especially for RNA-seq TPM/counts
        if not sp.issparse(data_df.values): # Check if not already sparse
             adata = ad.AnnData(X=sp.csr_matrix(data_df.values), obs=obs_df, var=var_df)
        else:
             adata = ad.AnnData(X=data_df.values, obs=obs_df, var=var_df) # If already sparse or small enough
    except Exception as e:
        logger.error(f"Error creating AnnData object for {dataset_name_for_log}: {e}")
        try: # Fallback to dense if sparse failed for some reason (should be rare for typical numeric data)
            logger.warning("Retrying AnnData creation with dense X matrix due to previous error.")
            adata = ad.AnnData(X=data_df.values, obs=obs_df, var=var_df)
        except Exception as e_dense:
            logger.error(f"Dense AnnData creation also failed for {dataset_name_for_log}: {e_dense}")
            raise e_dense

    # Populate adata.uns using dataset_specific_cfg
    if dataset_specific_cfg:
        for key, value in dataset_specific_cfg.items():
            if key == "obs_columns": # These are defaults for obs_df, not for direct .uns storage here.
                continue
            adata.uns[key] = value # Will be serialized by save_anndata -> ensure_serializable

    # Ensure a basic dataset_info structure exists in uns, even if cfg is minimal
    if 'dataset_info' not in adata.uns:
        adata.uns['dataset_info'] = {}

    # Populate calculated/summary fields into dataset_info
    adata.uns['dataset_info']['source_dataset_label'] = dataset_specific_cfg.get("label", dataset_name_for_log)
    adata.uns['dataset_info']['processed_version'] = dataset_specific_cfg.get("version", "1.0") # Our internal pipeline version
    adata.uns['dataset_info']['original_dataset_version'] = dataset_specific_cfg.get("dataset_version", "unknown")
    adata.uns['dataset_info']['sample_count'] = adata.n_obs
    adata.uns['dataset_info']['gene_count'] = adata.n_vars
    adata.uns['dataset_info']['data_type'] = obs_df['data_type'].unique()[0] if 'data_type' in obs_df and obs_df['data_type'].nunique() == 1 else dataset_specific_cfg.get('data_type', 'unknown')
    adata.uns['dataset_info']['expression_unit'] = obs_df['expression_unit'].unique()[0] if 'expression_unit' in obs_df and obs_df['expression_unit'].nunique() == 1 else dataset_specific_cfg.get('obs_columns',{}).get('expression_unit','unknown')


    # Set project-wide harmonized standards directly in .uns
    adata.uns["harmonized_reference_genome"] = "hg38"
    adata.uns["harmonized_gencode_version"] = "v24"

    # Original versions from config are already top-level in adata.uns if they existed in dataset_specific_cfg
    # If not, set them as 'unknown' or based on config.
    adata.uns.setdefault("original_reference_genome", dataset_specific_cfg.get("reference_genome", "unknown"))
    adata.uns.setdefault("original_gencode_version", dataset_specific_cfg.get("gencode_version", "unknown"))

    # Validation report will be added in Stage 2 by standardize_metadata.py script
    # Ontology mappings also added in Stage 2
    adata.uns["stage1_processing_date"] = pd.Timestamp.now().strftime("%Y-%m-%d")

    return adata



def save_anndata(adata, file_path):
    """Save AnnData object to file with validation and improved serialization."""
    try:
        logger.info(f"Preparing to save AnnData to {file_path}")

        if adata.n_vars == 0:
            logger.error("Cannot save AnnData: var DataFrame is empty!")
            return False
        if 'gene_id' not in adata.var.columns and len(adata.var_names) > 0:
             logger.warning("Reconstructing potentially missing 'gene_id' column in var from index.")
             adata.var['gene_id'] = adata.var_names.astype(str)
        elif adata.var.shape[1] == 0 and len(adata.var_names) > 0: # No columns but index exists
             logger.error("var DataFrame has no columns! Attempting reconstruction with 'gene_id'.")
             adata.var['gene_id'] = adata.var_names.astype(str)

        if 'gene_id' in adata.var.columns: # Ensure string type
            adata.var['gene_id'] = adata.var['gene_id'].astype(str)

        if adata.obs.index.name is not None and adata.obs.index.name in adata.obs.columns:
            new_name = f"original_obs_index_{adata.obs.index.name}"
            logger.warning(f"Fixing index/column name conflict in obs: Renaming column '{adata.obs.index.name}' to '{new_name}'")
            adata.obs = adata.obs.rename(columns={adata.obs.index.name: new_name})

        if adata.var.index.name is not None and adata.var.index.name in adata.var.columns:
            new_name = f"original_var_index_{adata.var.index.name}"
            logger.warning(f"Fixing index/column name conflict in var: Renaming column '{adata.var.index.name}' to '{new_name}'")
            adata.var = adata.var.rename(columns={adata.var.index.name: new_name})

        adata.obs.index = adata.obs.index.astype(str)
        adata.var.index = adata.var.index.astype(str)

        logger.info("Converting object columns in obs/var to string before saving...")
        for df_attr_name in ['obs', 'var']:
            df = getattr(adata, df_attr_name)
            # Create a new DataFrame for modifications to avoid SettingWithCopyWarning on slices
            new_df_data = {}
            object_cols_converted = False
            for col in df.columns:
                if df[col].dtype == 'object':
                    # Safely convert dicts/lists to JSON strings, then everything to string
                    def safe_stringify(x):
                        if isinstance(x, (dict, list)):
                            try:
                                return json.dumps(x, default=str)
                            except TypeError:
                                return str(x) # Fallback for unjsonable complex objects
                        return str(x) if pd.notna(x) else '' # Convert others to string, NaN to empty string

                    new_df_data[col] = df[col].apply(safe_stringify).astype(str)
                    logger.debug(f"Converted {df_attr_name} object column '{col}' to string.")
                    object_cols_converted = True
                else:
                    new_df_data[col] = df[col] # Keep non-object columns as is

            if object_cols_converted:
                setattr(adata, df_attr_name, pd.DataFrame(new_df_data, index=df.index))
        logger.info("Object column conversion complete.")

        logger.info("Ensuring adata.uns is serializable...")
        # Create a new dictionary for the serialized .uns to avoid modifying during iteration
        serializable_uns_dict = {}
        for key, value in adata.uns.items():
            try:
                serializable_uns_dict[key] = ensure_serializable(value)
            except Exception as e_uns_item:
                logger.error(f"Could not serialize adata.uns['{key}'] (type: {type(value)}). Error: {e_uns_item}. Storing as string representation.")
                serializable_uns_dict[key] = str(value)
        adata.uns.clear()
        adata.uns.update(serializable_uns_dict)
        logger.info("adata.uns converted to serializable format.")

        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        logger.info(f"Writing AnnData to {file_path}...")
        adata.write_h5ad(file_path)
        logger.info(f"Successfully wrote AnnData to {file_path}")

        logger.info(f"Verifying saved file: {file_path}")
        test_load = ad.read_h5ad(file_path)
        logger.info(
            f"Verification successful: Loaded AnnData with shape={test_load.shape}, "
            f"var columns ({len(test_load.var.columns)}): {list(test_load.var.columns)}, "
            f"obs columns ({len(test_load.obs.columns)}): {list(test_load.obs.columns)}"
        )
        if test_load.n_vars > 0 and test_load.var.shape[1] == 0 :
             logger.error("Verification FAILED: Loaded AnnData has empty var DataFrame columns despite having genes!")
             return False
        if test_load.n_obs > 0 and test_load.obs.shape[1] == 0:
             logger.error("Verification FAILED: Loaded AnnData has empty obs DataFrame columns despite having samples!")
             return False
        return True

    except Exception as e:
        logger.error(f"Unhandled error in save_anndata for {file_path}: {e}", exc_info=True)
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
    """Read an ENCODE TPM file - uses identify_columns."""
    try:
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        gene_id_col, expr_val_col = identify_columns(df) # expr_val_col is the identified expression column

        if gene_id_col is not None and expr_val_col is not None:
            df[gene_id_col] = df[gene_id_col].astype(str)
            # Keep only gene ID and the identified expression column
            result_df = df[[gene_id_col, expr_val_col]].set_index(gene_id_col)
            # RENAME the identified expression column to a STANDARD intermediate name "expression_value"
            result_df = result_df.rename(columns={expr_val_col: "expression_value"})
            logger.debug(f"Read {file_path}, gene_col='{gene_id_col}', expr_col='{expr_val_col}', renamed to 'expression_value'. Shape: {result_df.shape}")
            return result_df
        else:
            logger.warning(f"Could not identify gene_id and/or expression value columns in {file_path}. Gene: {gene_id_col}, Expr: {expr_val_col}")
            logger.warning(f"Columns present: {df.columns.tolist()}")
            return None
    except Exception as e:
        logger.warning(f"Error reading TPM file {file_path}: {e}")
        return None


def extract_encode_metadata(file_path, entex_full_metadata=None, encode_specific_cfg=None):
    """
    Extract metadata for an ENCODE or ENTEx file.
    Prioritizes specific metadata from entex_full_metadata if it's an ENTEx sample.
    Uses encode_specific_cfg for cell line info and general ENCODE defaults.
    """
    metadata = {}
    file_name = os.path.basename(file_path)
    file_id = file_name.split(".")[0]
    metadata["original_sample_id"] = file_id

    encode_cell_line_info_from_json = (encode_specific_cfg or {}).get("dataset_info", {}).get("cell_lines", {})

    is_entex_sample = "entex" in str(file_path).lower() # Check path first
    if entex_full_metadata and entex_full_metadata.get("sample_lookup") and file_id in entex_full_metadata["sample_lookup"]:
        is_entex_sample = True # Confirm if found in ENTEx lookup

    if is_entex_sample:
        metadata["data_source"] = "ENTEx"
        sample_info = (entex_full_metadata.get("sample_lookup", {})).get(file_id, {})
        if sample_info:
            logger.debug(f"Using specific ENTEx metadata for {file_id}")
            donor_id = sample_info.get("donor_id", sample_info.get("donor", sample_info.get("subject_id")))
            metadata["tissue"] = sample_info.get("tissue", "")
            metadata["subject_id"] = donor_id
            metadata["assay_type"] = sample_info.get("assay_type", sample_info.get("library_preparation", "unknown"))
            metadata["genome_annotation_source"] = sample_info.get("genome_annotation", "unknown")
            metadata["genome_assembly_source"] = sample_info.get("genome_assembly", "unknown")
            if donor_id and donor_id in entex_full_metadata.get("donor_map", {}):
                donor_info_map = entex_full_metadata["donor_map"][donor_id]
                metadata["sex"] = donor_info_map.get("sex", "")
                metadata["age"] = str(donor_info_map.get("age", ""))
        else:
            logger.warning(f"File {file_id} appears to be ENTEx by path, but no specific metadata in lookup. Using defaults.")
    else: # ENCODE Cell Line
        metadata["data_source"] = "ENCODE"
        path_parts = Path(file_path).parts
        cell_line_name_from_path = None
        for part in reversed(path_parts):
            for cl_key in encode_cell_line_info_from_json.keys():
                if cl_key.lower() in part.lower():
                    cell_line_name_from_path = cl_key
                    break
            if cell_line_name_from_path: break

        if cell_line_name_from_path:
            metadata["cell_line"] = cell_line_name_from_path
            cell_specific_info = encode_cell_line_info_from_json.get(cell_line_name_from_path, {})
            # Update metadata with specific cell line info from JSON if available
            for key_to_update in ['tissue', 'sex', 'age', 'organism', 'subject_id', 'disease', 'ethnicity']:
                if key_to_update in cell_specific_info:
                    metadata[key_to_update] = cell_specific_info[key_to_update]

            relevant_path_part = next((p for p in path_parts if cell_line_name_from_path.lower() in p.lower()), "")
            if "_polya_plus" in relevant_path_part.lower() or "_polya+" in relevant_path_part.lower():
                metadata["extraction_method"] = "polyA_plus"
            elif "_polya_minus" in relevant_path_part.lower() or "_polya-" in relevant_path_part.lower():
                metadata["extraction_method"] = "polyA_minus"
            elif "_total" in relevant_path_part.lower():
                metadata["extraction_method"] = "total"
        else:
            logger.warning(f"Could not determine cell line for ENCODE sample {file_id} from path. Defaults will be used.")

    # Apply general ENCODE defaults from encode_specific_cfg if not already set
    if encode_specific_cfg and "obs_columns" in encode_specific_cfg:
        for key, default_val in encode_specific_cfg["obs_columns"].items():
            metadata.setdefault(key, default_val)

    # Final fallback for core fields if still missing
    for fld in ['tissue', 'subject_id', 'sex', 'age', 'data_type', 'expression_unit']:
        metadata.setdefault(fld, "unknown" if fld != 'age' else "")
    metadata['age'] = str(metadata['age']) # Ensure age is string

    return metadata

def process_encode_data(
    input_dir, entex_dir=None, output_file=None, entex_metadata_file=None, metadata_dir=None, gencode_map=None
): # Added gencode_map parameter
    logger.info(f"Processing ENCODE & ENTEx data. Cell line dir: {input_dir}, ENTEx dir: {entex_dir}")
    dataset_name_for_log = "ENCODE_ENTEx_Combined"

    encode_specific_cfg = load_dataset_specific_metadata(metadata_dir, "encode") or {}
    entex_full_metadata = load_entex_metadata(entex_metadata_file) if (entex_metadata_file or entex_dir) else {}

    all_tsv_files_paths = []
    if input_dir and os.path.exists(input_dir):
        all_tsv_files_paths.extend(glob.glob(os.path.join(input_dir, "**", "*.tsv"), recursive=True))
    if entex_dir and os.path.exists(entex_dir) and entex_dir != input_dir :
        all_tsv_files_paths.extend(glob.glob(os.path.join(entex_dir, "**", "*.tsv"), recursive=True))

    all_tsv_files_paths = sorted(list(set(all_tsv_files_paths)))

    final_tsv_files_to_process = []
    for file_path_str in all_tsv_files_paths:
        file_id = Path(file_path_str).name.split(".")[0]
        is_potentially_entex = "entex" in file_path_str.lower() or (entex_full_metadata.get("sample_lookup") and file_id in entex_full_metadata["sample_lookup"])

        if is_potentially_entex:
            sample_info = (entex_full_metadata.get("sample_lookup", {})).get(file_id)
            if sample_info:
                quant_method = sample_info.get("quantification_method", "").lower()
                genome_annotation = sample_info.get("genome_annotation", "").upper()
                is_preferred_quant = "gene" in quant_method and "transcript" not in quant_method
                is_preferred_annotation = genome_annotation in ["V24", "V29", ""] or not genome_annotation
                if is_preferred_quant and is_preferred_annotation:
                    final_tsv_files_to_process.append(file_path_str)
                else:
                    logger.debug(f"Skipping ENTEx file {file_id} (from {file_path_str}) due to non-preferred quant/annotation: {quant_method}, {genome_annotation}")
            else:
                logger.warning(f"File {file_id} (from {file_path_str}) seems like ENTEx by path but no specific metadata. Including.")
                final_tsv_files_to_process.append(file_path_str)
        else:
            final_tsv_files_to_process.append(file_path_str)

    logger.info(f"Total files after initial scan: {len(all_tsv_files_paths)}. After filtering for ENCODE/ENTEx: {len(final_tsv_files_to_process)}")

    if not final_tsv_files_to_process:
        logger.error(f"No TSV files to process for ENCODE/ENTEx.")
        return None

    sample_dfs = {}
    sample_metadata_list = []
    all_gene_ids = set()

    for file_path_str in final_tsv_files_to_process:
        tpm_df = read_encode_tpm_file(file_path_str)
        if tpm_df is None:
            continue
        file_name = os.path.basename(file_path_str)
        sample_id = file_name.split(".")[0]

        if "expression_value" in tpm_df.columns:
            tpm_df = tpm_df.rename(columns={"expression_value": sample_id})
        else: # Should ideally not happen if read_encode_tpm_file worked
            logger.error(f"For {sample_id}, 'expression_value' column was expected but not found after read_encode_tpm_file. Columns: {tpm_df.columns.tolist()}")
            continue

        sample_dfs[sample_id] = tpm_df
        current_metadata = extract_encode_metadata(file_path_str, entex_full_metadata, encode_specific_cfg)
        current_metadata["sample_id_from_file"] = sample_id
        sample_metadata_list.append(current_metadata)
        all_gene_ids.update(tpm_df.index)

    if not sample_dfs:
        logger.error(f"No valid TPM data found for {dataset_name_for_log} (ENCODE/ENTEx)")
        return None

    num_processed_encode_entex_files = len(sample_dfs)
    logger.info(f"Successfully read data from {num_processed_encode_entex_files} ENCODE/ENTEx files.")

    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs for {dataset_name_for_log}")
    gene_id_mapping = {gid: standardize_ensembl_id(gid) for gid in all_gene_ids}

    valid_mapped_values_encode = [str(gid) for gid in gene_id_mapping.values() if gid is not None]
    if not valid_mapped_values_encode:
        logger.error(f"{dataset_name_for_log}: No gene IDs could be meaningfully standardized. Aborting processing for this dataset.")
        return None
    unique_std_ids = sorted(list(set(valid_mapped_values_encode)))

    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs for {dataset_name_for_log}")

    std_id_to_idx = {gid: i for i, gid in enumerate(unique_std_ids)}

    obs_df_encode = pd.DataFrame(sample_metadata_list)
    if 'sample_id_from_file' not in obs_df_encode.columns or obs_df_encode.empty:
        logger.error(f"{dataset_name_for_log}: 'sample_id_from_file' key missing or obs_df is empty.")
        return None

    obs_df_encode = obs_df_encode.set_index('sample_id_from_file')

    final_sample_ids_encode = obs_df_encode.index.tolist()
    if obs_df_encode.index.has_duplicates:
        logger.warning(f"{dataset_name_for_log}: Correcting duplicate sample IDs in obs_df index. Original count: {len(obs_df_encode)}")
        obs_df_encode = obs_df_encode[~obs_df_encode.index.duplicated(keep='first')]
        final_sample_ids_encode = obs_df_encode.index.tolist()

    sample_dfs = {sid: df_val for sid, df_val in sample_dfs.items() if sid in final_sample_ids_encode}

    num_genes_encode = len(unique_std_ids)
    num_samples_encode = len(final_sample_ids_encode)

    if num_samples_encode == 0:
        logger.error(f"{dataset_name_for_log}: No samples remaining after metadata processing and filtering. Aborting.")
        return None
    if num_genes_encode == 0:
        logger.error(f"{dataset_name_for_log}: No genes remaining after standardization. Aborting.")
        return None

    logger.info(f"Creating unified expression matrix for {dataset_name_for_log} ({num_genes_encode} genes x {num_samples_encode} samples)")

    expr_matrix = np.full((num_genes_encode, num_samples_encode), np.nan, dtype=np.float32)
    sample_id_to_col_idx = {sid: i for i, sid in enumerate(final_sample_ids_encode)}

    logger.info("DEBUG: Checking sample_dfs structure before ENCODE/ENTEx aggregation loop...")
    for temp_sample_id_debug, temp_df_debug in sample_dfs.items():
        if temp_df_debug is None:
            logger.info(f"DEBUG: sample_dfs['{temp_sample_id_debug}'] is None")
            continue
        # Ensure the column name is indeed the sample_id (which it should be after rename)
        if temp_sample_id_debug not in temp_df_debug.columns:
            logger.warning(f"DEBUG WARNING: Column '{temp_sample_id_debug}' NOT FOUND in sample_dfs['{temp_sample_id_debug}']. Columns are: {temp_df_debug.columns.tolist()}")

    std_to_originals_map_encode = {}
    for original_gid, mapped_std_id_val in gene_id_mapping.items():
        if mapped_std_id_val is not None:
            str_mapped_std_id_val = str(mapped_std_id_val)
            std_to_originals_map_encode.setdefault(str_mapped_std_id_val, []).append(str(original_gid))

    for std_id, original_ids_for_this_std in std_to_originals_map_encode.items():
        if std_id not in std_id_to_idx:
            logger.debug(f"Standardized ID '{std_id}' not in final index for {dataset_name_for_log}. Skipping.")
            continue
        gene_idx = std_id_to_idx[std_id]
        for sample_id, col_idx in sample_id_to_col_idx.items():
            sample_df = sample_dfs.get(sample_id)
            if sample_df is None: continue
            if sample_id not in sample_df.columns: # This check is crucial
                logger.error(f"CRITICAL DEBUG in ENCODE: Column '{sample_id}' missing in its DataFrame. Columns: {sample_df.columns.tolist()}")
                continue

            relevant_orig_ids_for_this_std_id = [orig_id for orig_id in original_ids_for_this_std if orig_id in sample_df.index]
            if relevant_orig_ids_for_this_std_id:
                values = sample_df.loc[relevant_orig_ids_for_this_std_id, sample_id]
                if not values.empty and not values.isna().all():
                    expr_matrix[gene_idx, col_idx] = values.max()

    expr_matrix = np.nan_to_num(expr_matrix, nan=0.0)
    unified_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=final_sample_ids_encode)

    var_df_encode = pd.DataFrame(index=unique_std_ids)
    var_df_encode["gene_id"] = var_df_encode.index.astype(str)
    var_df_encode["original_ids"] = var_df_encode.index.map(lambda x: ";".join(std_to_originals_map_encode.get(x, [])))

    # gencode_map = load_gencode_mapping() # Removed: use passed gencode_map
    var_df_encode = add_gencode_annotations(var_df_encode, gencode_map)

    mappings = load_mappings()
    obs_df_encode = standardize_metadata(obs_df_encode, "ENCODE", mappings)
    if obs_df_encode is None: # Check if standardize_metadata failed
        logger.error("ENCODE: standardize_metadata returned None. Aborting ENCODE/ENTEx processing.")
        return None


    # ---- ADD DEBUG LOGS BEFORE CREATING ANNDATA FOR ENCODE ---
    logger.info(f"ENCODE: Type of unified_expr_df.T: {type(unified_expr_df.T)}, Shape: {unified_expr_df.T.shape if hasattr(unified_expr_df, 'T') else 'N/A'}")
    logger.info(f"ENCODE: Type of obs_df_encode: {type(obs_df_encode)}, Shape: {obs_df_encode.shape if hasattr(obs_df_encode, 'shape') else 'N/A'}")
    logger.info(f"ENCODE: Type of var_df_encode: {type(var_df_encode)}, Shape: {var_df_encode.shape if hasattr(var_df_encode, 'shape') else 'N/A'}")

    if unified_expr_df is None: logger.error("CRITICAL ENCODE: unified_expr_df is None!")
    if obs_df_encode is None: logger.error("CRITICAL ENCODE: obs_df_encode is None!")
    if var_df_encode is None: logger.error("CRITICAL ENCODE: var_df_encode is None!")
    # ---- END DEBUG LOGS ---

    adata = create_standard_anndata(unified_expr_df.T, obs_df_encode, var_df_encode, encode_specific_cfg, dataset_name_for_log)

    if encode_specific_cfg:
        adata = apply_dataset_specific_metadata(adata, encode_specific_cfg)

    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed {dataset_name_for_log} data: {adata.n_obs} samples and {adata.n_vars} genes")
    else:
        logger.error(f"Failed to save AnnData for {dataset_name_for_log}")
        return None
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
        sample_attrs = pd.read_csv(GTEx_SAMPLE_ATTRS, sep="\t", low_memory=False)

        # Load subject attributes
        logger.info(f"Loading subject attributes from {GTEx_SUBJECT_ATTRS}")
        subject_attrs = pd.read_csv(GTEx_SUBJECT_ATTRS, sep="\t", low_memory=False)

        # Extract SUBJID from SAMPID
        sample_attrs["SUBJID"] = sample_attrs["SAMPID"].str.extract(r"(GTEX-[A-Z0-9]+)")

        # Merge sample and subject attributes
        merged_attrs = pd.merge(sample_attrs, subject_attrs, on="SUBJID", how="left")

        # Standardize column names based on a mapping
        column_mapping = {
            "SAMPID": "sample_id", # Will become index
            "SUBJID": "subject_id",
            "SMTSD": "tissue",  # Detailed tissue site
            "SMTS": "broad_tissue",  # Broad tissue type
            "SEX": "sex_code",
            "AGE": "age_bracket",
            "RACE": "race_code_source",
            "ETHNCTY": "ethnicity_code_source",
        }

        cols_to_rename = {k: v for k, v in column_mapping.items() if k in merged_attrs.columns}
        merged_attrs.rename(columns=cols_to_rename, inplace=True)

        gtex_cfg = load_dataset_specific_metadata(DEFAULT_METADATA_DIR, "gtex") or {}

        if 'race_code_source' in merged_attrs.columns and gtex_cfg.get('source_race_to_standard_label_map'):
            race_map = {str(k): v for k,v in gtex_cfg['source_race_to_standard_label_map'].items()}
            merged_attrs['race'] = merged_attrs['race_code_source'].astype(str).map(race_map).fillna('unknown or not reported')
            logger.info("Mapped GTEx 'race_code_source' to 'race' column.")
        else:
            merged_attrs['race'] = 'unknown or not reported'
            logger.warning(f"GTEx: Could not map race codes. Source column 'race_code_source' or mapping in JSON missing.")

        if 'ethnicity_code_source' in merged_attrs.columns and gtex_cfg.get('source_ethnicity_to_standard_label_map'):
            eth_map = {str(k): v for k,v in gtex_cfg['source_ethnicity_to_standard_label_map'].items()}
            merged_attrs['is_hispanic_or_latino'] = merged_attrs['ethnicity_code_source'].astype(str).map(eth_map).fillna('unknown or not reported')
            logger.info("Mapped GTEx 'ethnicity_code_source' to 'is_hispanic_or_latino' column.")
        else:
            merged_attrs['is_hispanic_or_latino'] = 'unknown or not reported'
            logger.warning(f"GTEx: Could not map ethnicity codes. Source column 'ethnicity_code_source' or mapping in JSON missing.")

        if 'sex_code' in merged_attrs.columns:
            sex_map = {'1': 'male', '2': 'female'}
            merged_attrs['sex'] = merged_attrs['sex_code'].astype(str).map(sex_map).fillna('unknown')
        elif 'sex' not in merged_attrs.columns:
            merged_attrs['sex'] = 'unknown'

        if 'age_bracket' in merged_attrs.columns:
            merged_attrs['age'] = merged_attrs['age_bracket'].astype(str)
        elif 'age' not in merged_attrs.columns:
            merged_attrs['age'] = ''

        if "tissue" not in merged_attrs.columns and "broad_tissue" in merged_attrs.columns:
            merged_attrs["tissue"] = merged_attrs["broad_tissue"]
            logger.info("Using 'broad_tissue' as 'tissue' column for GTEx.")
        elif "tissue" not in merged_attrs.columns:
             merged_attrs["tissue"] = "unknown"

        if "sample_id" in merged_attrs.columns:
            merged_attrs.set_index("sample_id", inplace=True)
        else:
            logger.error("GTEx metadata merging failed to produce 'sample_id' column from SAMPID.")
            return pd.DataFrame()

        logger.info(f"Loaded and initially processed metadata for {len(merged_attrs)} GTEx samples")
        return merged_attrs

    except Exception as e:
        logger.error(f"Error loading GTEx metadata: {e}", exc_info=True)
        return pd.DataFrame()


def process_gtex_data(input_file, output_file, metadata_dir=None, gencode_map=None): # Added gencode_map parameter
    logger.info(f"Processing GTEx data from {input_file} (Memory Optimized V2)")
    dataset_name_for_log = "GTEx"

    metadata_df = load_gtex_metadata()
    if metadata_df.empty:
        logger.error("Failed to load GTEx metadata. Cannot proceed.")
        return None
    logger.info(f"Loaded metadata for {len(metadata_df)} potential GTEx samples.")

    sample_ids_in_gct = []
    n_genes_in_gct = 0
    try:
        with gzip.open(input_file, "rt") as f:
            next(f)
            dims = next(f).strip().split()
            n_genes_in_gct = int(dims[0])
            header_line = next(f).strip().split("\t")
            sample_ids_in_gct = header_line[2:]
        logger.info(f"GCT header indicates {n_genes_in_gct} genes and {len(sample_ids_in_gct)} samples.")
    except Exception as e:
        logger.error(f"Error reading GCT header: {e}")
        return None

    common_samples_set = set(metadata_df.index).intersection(set(sample_ids_in_gct))
    if not common_samples_set:
        logger.error("No common samples found between metadata and GCT file.")
        return None

    final_sample_ids = sorted(list(common_samples_set))
    n_final_samples = len(final_sample_ids)
    logger.info(f"Processing {n_final_samples} common samples.")
    obs_df = metadata_df.loc[final_sample_ids].copy()

    obs_df['expression_unit'] = 'TPM'
    obs_df['data_type'] = 'RNA-seq'

    logger.info("Standardizing observation metadata for GTEx...")
    mappings = load_mappings()
    obs_df = standardize_metadata(obs_df, dataset_name_for_log, mappings)
    if obs_df is None: # Check if standardize_metadata failed
        logger.error("GTEx: standardize_metadata returned None. Aborting GTEx processing.")
        return None
    obs_df['sample_id'] = obs_df.index.astype(str)

    gct_sample_indices = {sid: i for i, sid in enumerate(sample_ids_in_gct)}
    sample_indices_to_keep = [gct_sample_indices[sid] for sid in final_sample_ids]

    logger.info("First pass: Identifying unique standardized genes...")
    std_gene_id_map = {}
    all_std_gene_ids = set()
    original_id_to_std_map = {}

    line_count = 0
    with gzip.open(input_file, "rt") as f:
        next(f); next(f); next(f)
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
    if not final_std_gene_ids:
        logger.error("GTEx: No genes could be standardized. Aborting GTEx processing.")
        return None
    n_final_genes = len(final_std_gene_ids)
    std_gene_id_to_idx = {gid: i for i, gid in enumerate(final_std_gene_ids)}
    logger.info(f"Identified {n_final_genes} unique standardized genes.")

    logger.info("Second pass: Populating expression matrix...")
    final_expr_matrix = np.full((n_final_genes, n_final_samples), -1.0, dtype=np.float32)

    line_count = 0
    with gzip.open(input_file, "rt") as f:
        next(f); next(f); next(f)
        for line in f:
            line_count += 1
            try:
                fields = line.strip().split("\t")
                original_gene_id = fields[0]
                std_gene_id = std_gene_id_map.get(original_gene_id)

                if not std_gene_id: continue

                std_idx = std_gene_id_to_idx.get(std_gene_id)
                if std_idx is None: continue

                expression_values_all = fields[2:]
                expression_values_filtered = np.array([float(expression_values_all[i]) for i in sample_indices_to_keep], dtype=np.float32)

                current_row = final_expr_matrix[std_idx, :]
                update_mask = (expression_values_filtered > current_row) | (current_row == -1.0)
                final_expr_matrix[std_idx, update_mask] = expression_values_filtered[update_mask]

            except Exception as line_e:
                logger.warning(f"Skipping line {line_count} during second pass due to error: {line_e}")
                continue
            if line_count % 5000 == 0:
                logger.info(f"Second pass processed {line_count}/{n_genes_in_gct} genes...")

    final_expr_matrix[final_expr_matrix == -1.0] = 0.0
    logger.info("Finished populating expression matrix.")

    var_df = pd.DataFrame(index=final_std_gene_ids)
    var_df['gene_id'] = var_df.index.astype(str)
    var_df['original_ids'] = [";".join(original_id_to_std_map.get(gid, [])) for gid in final_std_gene_ids]
    logger.info("Adding GENCODE annotations...")
    # gencode_mapping = load_gencode_mapping() # Removed: use passed gencode_map
    var_df = add_gencode_annotations(var_df, gencode_map)

    gtex_specific_cfg = load_dataset_specific_metadata(metadata_dir, "gtex") or {}
    adata = create_standard_anndata(pd.DataFrame(final_expr_matrix.T, index=final_sample_ids, columns=final_std_gene_ids),
                                    obs_df, var_df, gtex_specific_cfg, dataset_name_for_log)
    logger.info(f"Created AnnData object with shape {adata.shape}")

    if gtex_specific_cfg:
        adata = apply_dataset_specific_metadata(adata, gtex_specific_cfg)

    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed GTEx data with {adata.n_obs} samples and {adata.n_vars} genes (Memory Optimized V2)")
    else:
        logger.error("Failed to save memory-optimized GTEx AnnData V2.")
        return None

    return adata

# ====================================================================================
# MAGE-specific Processing
# ====================================================================================

def process_mage_dir(file_path, donor_id, tissue, sample_dfs, sample_metadata_list, all_gene_ids):
    try:
        file_name = os.path.basename(file_path)
        sample_id = f"{donor_id}_{os.path.splitext(file_name)[0]}"
        logger.debug(f"Processing MAGE file: {file_path} as sample_id: {sample_id}")

        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            logger.warning(f"File does not exist or is empty: {file_path}")
            return False

        with open(file_path, "r") as f:
            first_line = f.readline().strip()
        delimiter = "\t" if "\t" in first_line else ","

        df = pd.read_csv(file_path, delimiter=delimiter, low_memory=False)
        if df.empty:
            logger.warning(f"Empty DataFrame after reading MAGE file: {file_path}")
            return False

        gene_id_col_name, expr_col_name = identify_columns(df)
        expr_type = "unknown_expression"
        if expr_col_name and "TPM" in expr_col_name.upper():
            expr_type = "TPM"

        if not (gene_id_col_name and expr_col_name and gene_id_col_name in df.columns and expr_col_name in df.columns):
            logger.warning(f"Could not reliably identify gene ID ('{gene_id_col_name}') and/or expression ('{expr_col_name}') columns in {file_path}. Columns: {df.columns.tolist()}")
            return False

        gene_id_series = df[gene_id_col_name].astype(str)
        expr_series = pd.to_numeric(df[expr_col_name], errors='coerce')

        valid_expr_mask = expr_series.notna()
        if not valid_expr_mask.any():
            logger.warning(f"No valid numeric expression values in column '{expr_col_name}' for {file_path}.")
            return False

        gene_id_series = gene_id_series[valid_expr_mask]
        expr_series = expr_series[valid_expr_mask]

        if gene_id_series.empty:
            logger.warning(f"Gene ID series became empty after filtering for valid expression in {file_path}.")
            return False

        simple_df = pd.DataFrame({sample_id: expr_series.values}, index=gene_id_series.values)

        if simple_df.index.has_duplicates:
            logger.warning(f"Duplicate gene IDs found in {file_path}. Aggregating by max.")
            simple_df = simple_df.groupby(simple_df.index).max()

        sample_dfs[sample_id] = simple_df
        all_gene_ids.update(simple_df.index)

        metadata = {
            "sample_id_from_file": sample_id,
            "donor_id": donor_id, "subject_id": donor_id, "tissue": tissue,
            "dataset": "MAGE", "data_type": "RNA-seq", "expression_unit": expr_type,
            "is_gencode": "gencode" in file_name.lower(),
            # Add placeholders for sex and ethnicity, to be filled from 1000G data
            "sex": "unknown",
            "race": "unknown or not reported",
            "is_hispanic_or_latino": "unknown or not reported",
            "population_code_1000g": "unknown"
        }
        sample_metadata_list.append(metadata)
        logger.debug(f"Successfully processed {file_path} with {len(simple_df)} genes for sample {sample_id}")
        return True

    except Exception as e:
        logger.error(f"Error processing MAGE file {file_path}: {e}", exc_info=True)
        return False


def process_mage_data(input_dir, output_file, metadata_dir=None, gencode_map=None): # Added gencode_map
    logger.info(f"Processing MAGE data from {input_dir}")
    dataset_name_for_log = "MAGE"
    mage_specific_cfg = load_dataset_specific_metadata(metadata_dir, "mage") or {}

    # --- Load 1000 Genomes Metadata ---
    genomes_metadata_file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/project_gene_regulation/data/MAGE/WGS/20130606_g1k_3202_samples_ped_population.txt'
    genomes_metadata_map = {}
    if os.path.exists(genomes_metadata_file_path):
        try:
            logger.info(f"MAGE: Loading 1000 Genomes metadata from {genomes_metadata_file_path}")
            ped_df = pd.read_csv(genomes_metadata_file_path, sep='\t', dtype=str)
            ped_df.columns = ped_df.columns.str.strip()

            required_ped_cols_for_mapping = ['SampleID', 'Sex', 'Population']
            missing_cols = [col for col in required_ped_cols_for_mapping if col not in ped_df.columns]

            if missing_cols:
                logger.error(f"MAGE: 1000G PED file missing required columns for mapping after stripping: {missing_cols}. Columns found: {ped_df.columns.tolist()}")
            else:
                for _, row in ped_df.iterrows():
                    sample_id_1000g = str(row['SampleID']).strip().upper()
                    sex_code = str(row['Sex']).strip()
                    population_code = str(row['Population']).strip()
                    genomes_metadata_map[sample_id_1000g] = {
                        'sex_code_1000g': sex_code,
                        'population_code_1000g': population_code
                    }
                logger.info(f"MAGE: Loaded metadata for {len(genomes_metadata_map)} samples from 1000 Genomes PED file.")
                # ---- MAGE DEBUGGING ----
                if genomes_metadata_map:
                    logger.info(f"MAGE DEBUG: First 5 keys in genomes_metadata_map: {list(genomes_metadata_map.keys())[:5]}")
                # ---- END MAGE DEBUGGING ----

        except Exception as e:
            logger.error(f"MAGE: Error loading or processing 1000 Genomes metadata file {genomes_metadata_file_path}: {e}", exc_info=True)
    else:
        logger.warning(f"MAGE: 1000 Genomes metadata file not found at {genomes_metadata_file_path}. Sex and ethnicity will be 'unknown'.")




    sample_dfs = {}
    sample_metadata_list = []
    all_gene_ids = set()
    processed_samples_count = 0
    donor_dirs = glob.glob(os.path.join(input_dir, "NA*")) + glob.glob(os.path.join(input_dir, "HG*"))

    if not donor_dirs:
        logger.error(f"No donor directories (NA* or HG*) found in MAGE input directory: {input_dir}")
        return None

    for donor_dir_path in donor_dirs:
        donor_id_original_case = os.path.basename(donor_dir_path)
        donor_id_lookup_key = donor_id_original_case.upper() # Normalize for lookup
        # ---- MAGE DEBUGGING ----
        if processed_samples_count < 5: # Log for the first few samples
            logger.info(f"MAGE DEBUG: Processing donor_id_original_case: '{donor_id_original_case}', lookup_key: '{donor_id_lookup_key}'")
        # ---- END MAGE DEBUGGING ----


        lymphoblast_dir = os.path.join(donor_dir_path, "lymphoblast")
        if not os.path.isdir(lymphoblast_dir):
            logger.warning(f"No 'lymphoblast' subdirectory for MAGE donor {donor_id_original_case}, skipping.")
            continue

        all_csv_files_in_dir = glob.glob(os.path.join(lymphoblast_dir, "*.csv"))
        preferred_file = None
        gencode_v24_pruned_files = [f for f in all_csv_files_in_dir if Path(f).name.endswith("_gencode_V24_pruned.csv")]

        if gencode_v24_pruned_files:
            preferred_file = gencode_v24_pruned_files[0]
            if len(gencode_v24_pruned_files) > 1:
                logger.warning(f"Multiple '*_gencode_V24_pruned.csv' files for MAGE donor {donor_id_original_case} in {lymphoblast_dir}. Using {Path(preferred_file).name}")
        else:
            simple_name_files = [f for f in all_csv_files_in_dir if Path(f).name == f"{donor_id_original_case}.csv"]
            if simple_name_files:
                preferred_file = simple_name_files[0]
                logger.info(f"Using '{Path(preferred_file).name}' for MAGE donor {donor_id_original_case} as no '*_gencode_V24_pruned.csv' was found.")
                if len(simple_name_files) > 1:
                     logger.warning(f"Multiple '{donor_id_original_case}.csv' files found (unexpected). Using the first one.")
            elif all_csv_files_in_dir:
                preferred_file = all_csv_files_in_dir[0]
                logger.warning(f"No '*_gencode_V24_pruned.csv' or '{donor_id_original_case}.csv' file for MAGE donor {donor_id_original_case}. Using first available CSV: {Path(preferred_file).name}.")

        if preferred_file:
            logger.info(f"Selected file for MAGE donor {donor_id_original_case}: {Path(preferred_file).name}")
            if process_mage_dir(preferred_file, donor_id_original_case, "lymphoblast", sample_dfs, sample_metadata_list, all_gene_ids):
                processed_samples_count +=1
                if sample_metadata_list:
                    entry_to_update = next((item for item in reversed(sample_metadata_list) if item.get('donor_id') == donor_id_original_case), None)
                    if entry_to_update:
                        sample_1000g_info = genomes_metadata_map.get(donor_id_lookup_key) # Use normalized key for lookup
                        if sample_1000g_info:
                            logger.info(f"MAGE SUCCESS: Found 1000G metadata for {donor_id_original_case} using key {donor_id_lookup_key}") # Add success log
                            sex_code_1000g = sample_1000g_info.get('sex_code_1000g')
                            pop_code_1000g = sample_1000g_info.get('population_code_1000g')
                            sex_map_1000g = {'1': 'male', '2': 'female'}
                            entry_to_update['sex'] = sex_map_1000g.get(sex_code_1000g, 'unknown')
                            pop_to_race_map = mage_specific_cfg.get('pop_to_race_map', {'ACB':'black or african american', 'ASW':'black or african american', 'ESN':'black or african american', 'GWD':'black or african american', 'LWK':'black or african american', 'MSL':'black or african american', 'YRI':'black or african american', 'CLM':'multiethnic', 'MXL':'multiethnic', 'PEL':'multiethnic', 'PUR':'multiethnic', 'CDX':'asian', 'CHB':'asian', 'CHS':'asian', 'JPT':'asian', 'KHV':'asian', 'CEU':'white', 'FIN':'white', 'GBR':'white', 'IBS':'white', 'TSI':'white', 'BEB':'asian', 'GIH':'asian', 'ITU':'asian', 'PJL':'asian', 'STU':'asian'})
                            pop_to_hispanic_map = mage_specific_cfg.get('pop_to_hispanic_map', {'CLM':'hispanic or latino', 'MXL':'hispanic or latino', 'PEL':'hispanic or latino', 'PUR':'hispanic or latino'})
                            default_race = mage_specific_cfg.get('default_race_for_unmapped_pop', 'unknown or not reported')
                            derived_race = pop_to_race_map.get(pop_code_1000g, default_race)
                            derived_is_hispanic = pop_to_hispanic_map.get(pop_code_1000g, 'not hispanic or latino')
                            if derived_race == default_race and derived_is_hispanic == 'hispanic or latino': derived_race = 'multiethnic'
                            entry_to_update['race'] = derived_race
                            entry_to_update['is_hispanic_or_latino'] = derived_is_hispanic
                            entry_to_update['population_code_1000g'] = pop_code_1000g
                        else:
                            logger.warning(f"MAGE: No 1000G metadata found for sample {donor_id_original_case} (lookup key: {donor_id_lookup_key}). Sex/Ethnicity remain defaults.")
                    else:
                        logger.error(f"MAGE: Could not find metadata entry for {donor_id_original_case} to enrich with 1000G data.")
        else:
            logger.warning(f"No suitable CSV file found for MAGE donor {donor_id_original_case} in {lymphoblast_dir}")

    if processed_samples_count == 0:
        logger.error(f"No MAGE samples were successfully processed from {input_dir}.")
        return None
    logger.info(f"Processed {processed_samples_count} MAGE samples.")

    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs for {dataset_name_for_log}")
    gene_id_mapping = {gid: standardize_ensembl_id(gid) for gid in all_gene_ids}

    valid_mapped_values_mage = [str(gid) for gid in gene_id_mapping.values() if gid is not None]
    if not valid_mapped_values_mage:
        logger.error(f"{dataset_name_for_log}: No gene IDs could be meaningfully standardized. Aborting.")
        return None
    unique_std_ids = sorted(list(set(valid_mapped_values_mage)))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs for {dataset_name_for_log}")

    std_id_to_idx = {gid: i for i, gid in enumerate(unique_std_ids)}

    obs_df_mage = pd.DataFrame(sample_metadata_list)
    if 'sample_id_from_file' not in obs_df_mage.columns or obs_df_mage.empty:
        logger.error(f"{dataset_name_for_log}: 'sample_id_from_file' key missing or obs_df is empty.")
        return None
    obs_df_mage = obs_df_mage.set_index('sample_id_from_file')
    obs_df_mage.index.name = "sample_id"

    final_sample_ids_mage = obs_df_mage.index.tolist()
    if obs_df_mage.index.has_duplicates:
        logger.warning(f"{dataset_name_for_log}: Correcting duplicate sample IDs in obs_df index.")
        obs_df_mage = obs_df_mage[~obs_df_mage.index.duplicated(keep='first')]
        final_sample_ids_mage = obs_df_mage.index.tolist()

    sample_dfs = {sid: df_val for sid, df_val in sample_dfs.items() if sid in final_sample_ids_mage}

    num_genes_mage = len(unique_std_ids)
    num_samples_mage = len(final_sample_ids_mage)

    if num_samples_mage == 0 or num_genes_mage == 0:
        logger.error(f"{dataset_name_for_log}: No samples ({num_samples_mage}) or genes ({num_genes_mage}) remaining. Aborting.")
        return None

    logger.info(f"Creating unified expression matrix for {dataset_name_for_log} ({num_genes_mage} genes x {num_samples_mage} samples) using optimized method.")

    # Optimized Matrix Creation for MAGE
    all_series = []
    for sample_id, s_df in sample_dfs.items():
        series = s_df[sample_id].rename(sample_id)
        all_series.append(series)

    if not all_series:
        logger.error(f"No data series collected for {dataset_name_for_log}. Aborting matrix creation.")
        unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_mage).astype(np.float32)

    else:
        logger.info(f"Concatenating {len(all_series)} sample series for {dataset_name_for_log}...")
        try:
            concat_df = pd.concat(all_series, axis=1, join='outer').astype(np.float32)
            logger.info(f"Concatenated DataFrame shape: {concat_df.shape}")

            std_to_originals_map_mage = {}
            for original_gid, mapped_std_id_val in gene_id_mapping.items():
                if mapped_std_id_val is not None:
                    std_to_originals_map_mage.setdefault(str(mapped_std_id_val), []).append(str(original_gid))

            original_to_std_series_map = {orig_id: std_id_val for std_id_val, orig_ids_list in std_to_originals_map_mage.items() for orig_id in orig_ids_list}
            mappable_original_ids = [idx for idx in concat_df.index if idx in original_to_std_series_map]

            if not mappable_original_ids:
                logger.warning(f"No mappable original gene IDs in concat_df for {dataset_name_for_log}. Resulting matrix might be mostly zeros or have unexpected shape.")
                unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_mage).astype(np.float32)
            else:
                filtered_concat_df = concat_df.loc[mappable_original_ids]
                std_ids_for_aggregation = filtered_concat_df.index.map(original_to_std_series_map)

                if std_ids_for_aggregation.isna().all():
                    logger.error(f"All standardized IDs for aggregation are NaN for {dataset_name_for_log}. Check mapping.")
                    unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_mage).astype(np.float32)
                else:
                    filtered_concat_df.index = std_ids_for_aggregation
                    logger.info(f"Grouping by standardized gene IDs and aggregating for {dataset_name_for_log}...")
                    unified_expr_df_grouped = filtered_concat_df.groupby(level=0).max()
                    logger.info(f"Aggregation complete. Shape: {unified_expr_df_grouped.shape}")
                    unified_expr_df = unified_expr_df_grouped.reindex(index=unique_std_ids, columns=final_sample_ids_mage, fill_value=0.0).astype(np.float32)
        except Exception as e_concat:
            logger.error(f"Error during optimized matrix creation for {dataset_name_for_log}: {e_concat}", exc_info=True)
            logger.warning("Falling back to previous matrix creation method for MAGE due to error.")
            # Fallback to original loop method if concat/groupby fails
            expr_matrix = np.full((num_genes_mage, num_samples_mage), np.nan, dtype=np.float32)
            sample_id_to_col_idx = {sid: i for i, sid in enumerate(final_sample_ids_mage)}
            std_to_originals_map_mage = {} # Re-create as it might not be in scope if optimized failed early
            for original_gid, mapped_std_id_val in gene_id_mapping.items():
                if mapped_std_id_val is not None:
                    std_to_originals_map_mage.setdefault(str(mapped_std_id_val), []).append(str(original_gid))

            for std_id, original_ids_for_this_std in std_to_originals_map_mage.items():
                if std_id not in std_id_to_idx:
                    logger.warning(f"Standardized ID '{std_id}' not in final index for MAGE (fallback). Skipping.")
                    continue
                gene_idx = std_id_to_idx[std_id]
                for sample_id, col_idx in sample_id_to_col_idx.items():
                    sample_df = sample_dfs.get(sample_id)
                    if sample_df is None: continue
                    if sample_id not in sample_df.columns:
                        logger.error(f"MAGE Aggregation Error (fallback): Column '{sample_id}' not found in sample_df. Columns: {sample_df.columns.tolist()}")
                        continue
                    relevant_orig_ids_in_sample = [orig_id for orig_id in original_ids_for_this_std if orig_id in sample_df.index]
                    if relevant_orig_ids_in_sample:
                        values = pd.to_numeric(sample_df.loc[relevant_orig_ids_in_sample, sample_id], errors='coerce')
                        if not values.empty and not values.isna().all():
                            expr_matrix[gene_idx, col_idx] = values.max()
            expr_matrix = np.nan_to_num(expr_matrix, nan=0.0)
            unified_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=final_sample_ids_mage)
    # End of Optimized Matrix Creation

    var_df_mage = pd.DataFrame(index=unique_std_ids)
    var_df_mage["gene_id"] = var_df_mage.index.astype(str)
    var_df_mage["original_ids"] = var_df_mage.index.map(lambda x: ";".join(std_to_originals_map_mage.get(x, [])))

    # gencode_map passed as parameter
    var_df_mage = add_gencode_annotations(var_df_mage, gencode_map)

    mappings = load_mappings()
    obs_df_mage = standardize_metadata(obs_df_mage, dataset_name_for_log, mappings)
    if obs_df_mage is None: # Check if standardize_metadata failed
        logger.error("MAGE: standardize_metadata returned None. Aborting MAGE processing.")
        return None


    adata = create_standard_anndata(unified_expr_df.T, obs_df_mage, var_df_mage, mage_specific_cfg, dataset_name_for_log)

    if mage_specific_cfg:
        adata = apply_dataset_specific_metadata(adata, mage_specific_cfg)

    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed MAGE data: {adata.n_obs} samples and {adata.n_vars} genes")
    else:
        logger.error(f"Failed to save AnnData for MAGE.")
        return None
    return adata


# ====================================================================================
# ADNI-specific Processing
# ====================================================================================

def preprocess_adni_file(file_path):
    """Preprocess ADNI file to fix escaped tab characters if present."""
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        if '\\t' in content:
            logger.info(f"Fixing escaped tabs '\\t' in {file_path}")
            fixed_content = content.replace('\\t', '\t')

            fixed_file_path = file_path + '.fixed_tabs'
            with open(fixed_file_path, 'w', encoding='utf-8') as f_out:
                f_out.write(fixed_content)
            return fixed_file_path
    except Exception as e:
        logger.warning(f"Error during ADNI file preprocessing (tab fixing) for {file_path}: {e}")

    return file_path

def process_adni_data(input_dir, output_file, adni_demographics_file=None, adni_diagnosis_file=None, dict_file=None, metadata_dir=None, gencode_map=None): # Added gencode_map
    logger.info(f"Processing ADNI data from {input_dir} to capture longitudinal info.")
    dataset_name_for_log = "ADNI"
    adni_specific_cfg = load_dataset_specific_metadata(metadata_dir, "adni") or {}

    sample_dfs = {}
    sample_metadata_list = []
    all_gene_ids = set()
    processed_expression_files_count = 0

    demographics_df = None
    if adni_demographics_file and os.path.exists(adni_demographics_file):
        logger.info(f"ADNI: Loading demographics metadata from {adni_demographics_file}")
        try:
            demographics_df = pd.read_csv(adni_demographics_file, low_memory=False)
            if 'RID' in demographics_df.columns:
                demographics_df['RID'] = pd.to_numeric(demographics_df['RID'], errors='coerce').astype('Int64')
            if 'PTID' in demographics_df.columns:
                demographics_df['PTID'] = demographics_df['PTID'].astype(str)
            demographics_df = demographics_df.drop_duplicates(subset=['RID' if 'RID' in demographics_df.columns else 'PTID'], keep='first')
        except Exception as e:
            logger.error(f"ADNI: Error loading demographics file {adni_demographics_file}: {e}", exc_info=True)
            demographics_df = None

    diagnosis_long_df = None
    all_subject_diagnoses_for_uns = {}

    if adni_diagnosis_file and os.path.exists(adni_diagnosis_file):
        logger.info(f"ADNI: Loading longitudinal diagnosis data from {adni_diagnosis_file}")
        try:
            dx_df_raw = pd.read_csv(adni_diagnosis_file, low_memory=False)
            if 'RID' not in dx_df_raw.columns:
                logger.error("ADNI Diagnosis file must contain 'RID' column.")
                all_subject_diagnoses_for_uns = {}
                diagnosis_long_df = None
            else:
                dx_df_raw['EXAMDATE'] = pd.to_datetime(dx_df_raw['EXAMDATE'], errors='coerce')
                diagnosis_long_df = dx_df_raw.copy()
                diagnosis_long_df['RID'] = pd.to_numeric(diagnosis_long_df['RID'], errors='coerce').astype('Int64')
                diagnosis_long_df['DIAGNOSIS'] = pd.to_numeric(diagnosis_long_df['DIAGNOSIS'], errors='coerce')
                if 'VISCODE2' in diagnosis_long_df.columns:
                    diagnosis_long_df['VISCODE2'] = diagnosis_long_df['VISCODE2'].astype(str).str.lower().str.strip()
                logger.info(f"ADNI: Loaded {len(diagnosis_long_df)} total diagnosis entries.")
                diag_map_codes = adni_specific_cfg.get("diagnosis_mapping", {
                    1: "Cognitively Normal", 2: "Mild Cognitive Impairment", 3: "Alzheimer's Disease",
                    -4: "Missing/Unknown"
                })
                dx_df_raw_cleaned = dx_df_raw.dropna(subset=['RID', 'EXAMDATE', 'DIAGNOSIS'])
                dx_df_raw_cleaned.loc[:, 'RID'] = pd.to_numeric(dx_df_raw_cleaned['RID'], errors='coerce').astype('Int64')
                for rid_val, group in dx_df_raw_cleaned.groupby('RID'):
                    subject_diagnoses = []
                    for _, row in group.sort_values(by='EXAMDATE').iterrows():
                        diag_code = row['DIAGNOSIS']
                        diag_label = diag_map_codes.get(diag_code, f"Unknown code: {diag_code}")
                        exam_date_val = row['EXAMDATE']
                        subject_diagnoses.append({
                            "visit_code": str(row.get('VISCODE2', 'N/A')),
                            "exam_date": str(exam_date_val.date()) if pd.notna(exam_date_val) and hasattr(exam_date_val, 'date') else 'N/A',
                            "diagnosis_code": int(diag_code) if pd.notna(diag_code) else None,
                            "diagnosis_label": diag_label
                        })
                    if subject_diagnoses:
                        all_subject_diagnoses_for_uns[str(int(rid_val))] = subject_diagnoses
                logger.info(f"ADNI: Processed full diagnostic history for {len(all_subject_diagnoses_for_uns)} subjects for adata.uns.")
        except Exception as e:
            logger.error(f"ADNI: Error loading or processing diagnosis file {adni_diagnosis_file}: {e}", exc_info=True)
            diagnosis_long_df = None
            all_subject_diagnoses_for_uns = {}

    subject_dirs = glob.glob(os.path.join(input_dir, "*_S_*"))
    expression_files = []
    for subject_dir_path in subject_dirs:
        expression_files.extend(glob.glob(os.path.join(subject_dir_path, "*.txt")))
        expression_files.extend(glob.glob(os.path.join(subject_dir_path, "*.csv")))
        expression_files.extend(glob.glob(os.path.join(subject_dir_path, "**", "*.txt"), recursive=True))
        expression_files.extend(glob.glob(os.path.join(subject_dir_path, "**", "*.csv"), recursive=True))

    expression_files.extend(glob.glob(os.path.join(input_dir, "*.txt")))
    expression_files.extend(glob.glob(os.path.join(input_dir, "*.csv")))
    expression_files = sorted(list(set(expression_files)))

    if not expression_files:
        logger.error(f"No expression files (.txt or .csv) found for ADNI in {input_dir} or its subdirectories.")
        return None

    for file_path_str in expression_files:
        if "ADNI_Gene_Expression_Profile.csv" in os.path.basename(file_path_str) or \
           "ADNI_Gene_Expression_Profile_DICT.csv" in os.path.basename(file_path_str):
            logger.info(f"ADNI: Skipping summary/dictionary file {file_path_str} in per-sample processing loop.")
            continue

        subject_id_from_path = "unknown_subject"
        path_obj = Path(file_path_str)
        parent_dir_name = path_obj.parent.name
        if re.match(r"^\d{3}_S_\d{4,}$", parent_dir_name):
            subject_id_from_path = parent_dir_name
        else:
            match_fn = re.search(r"(\d{3}_S_\d{4,})", path_obj.name)
            if match_fn:
                subject_id_from_path = match_fn.group(1)
            else:
                logger.warning(f"Could not determine ADNI subject ID for file: {file_path_str}. Using 'unknown_subject'.")

        file_name = os.path.basename(file_path_str)
        sample_id_suffix = Path(file_name).stem
        sample_id = f"{subject_id_from_path}_{sample_id_suffix}"

        processed_file_path = preprocess_adni_file(file_path_str)

        try:
            df = None
            file_basename_lower = os.path.basename(processed_file_path).lower()
            if "_gencode_v24_pruned.csv" in file_basename_lower:
                logger.debug(f"ADNI: Attempting tab delimiter for pruned CSV: {processed_file_path}")
                try: df = pd.read_csv(processed_file_path, sep='\t', comment='#', low_memory=False)
                except Exception:
                    logger.warning(f"Tab delimiter failed for pruned CSV {processed_file_path}. Trying comma.")
                    try: df = pd.read_csv(processed_file_path, sep=',', comment='#', low_memory=False)
                    except Exception as e_comma:
                        logger.error(f"Comma delimiter also failed for pruned CSV {processed_file_path}: {e_comma}. Skipping.")
                        if processed_file_path != file_path_str: os.remove(processed_file_path)
                        continue
            elif file_basename_lower.endswith('.csv'):
                try: df = pd.read_csv(processed_file_path, sep=',', comment='#', low_memory=False)
                except Exception:
                    logger.warning(f"Comma delimiter failed for CSV {processed_file_path}. Trying tab.")
                    try: df = pd.read_csv(processed_file_path, sep='\t', comment='#', low_memory=False)
                    except Exception as e_tab:
                        logger.error(f"Tab delimiter also failed for CSV {processed_file_path}: {e_tab}. Skipping.")
                        if processed_file_path != file_path_str: os.remove(processed_file_path)
                        continue
            else:
                df = pd.read_csv(processed_file_path, sep='\t', comment='#', low_memory=False)

            if df is None:
                if processed_file_path != file_path_str: os.remove(processed_file_path)
                continue
            if processed_file_path != file_path_str and os.path.exists(processed_file_path):
                os.remove(processed_file_path)

            if df.empty:
                logger.warning(f"Empty data after reading ADNI file: {file_path_str}")
                continue

            gene_id_col, expr_val_col = identify_columns(df)
            if gene_id_col is None or expr_val_col is None:
                logger.warning(f"Could not identify required columns in ADNI file {file_path_str}. Gene: {gene_id_col}, Expr: {expr_val_col}. Columns: {df.columns.tolist()}")
                continue

            df[gene_id_col] = df[gene_id_col].astype(str)
            df[expr_val_col] = pd.to_numeric(df[expr_val_col], errors='coerce')
            simple_df = df.set_index(gene_id_col)[[expr_val_col]]
            simple_df = simple_df.rename(columns={expr_val_col: sample_id})
            sample_dfs[sample_id] = simple_df
            all_gene_ids.update(simple_df.index)

            current_metadata = {
                "sample_id_from_file": sample_id, "subject_id": subject_id_from_path,
                "dataset": adni_specific_cfg.get("label", "ADNI"),
                "data_type": adni_specific_cfg.get("data_type", "Microarray"),
                "expression_unit": adni_specific_cfg.get("obs_columns", {}).get("expression_unit", "Normalized intensity"),
                "tissue": adni_specific_cfg.get("obs_columns", {}).get("tissue", "blood"),
                "platform": adni_specific_cfg.get("platform", "Affymetrix Human Genome U219 Array"),
                "expression_sample_visit_code": "unknown", "age_at_expression_sampling": "",
                "diagnosis_at_expression_sampling_code": pd.NA,
                "diagnosis_at_expression_sampling_label": "Unknown",
                "date_of_expression_sampling": ""
            }
            parsed_visit_code = "unknown"
            filename_lower = sample_id_suffix.lower()
            visit_pattern = r"_(bl|sc|scrn|screening|m\d{1,3}|v\d{1,3}|aim\d{1,3})[._]"
            visit_match = re.search(visit_pattern, filename_lower)
            if visit_match:
                parsed_visit_code = visit_match.group(1)
                if parsed_visit_code.startswith("aim"): parsed_visit_code = "m" + parsed_visit_code[3:].zfill(2)
                elif parsed_visit_code in ["scrn", "screening"]: parsed_visit_code = "sc"
            current_metadata["expression_sample_visit_code"] = parsed_visit_code
            subject_rid = None
            try:
                if "_S_" in subject_id_from_path: subject_rid = int(subject_id_from_path.split('_S_')[-1])
            except ValueError: logger.warning(f"Could not parse RID for {subject_id_from_path}")
            if demographics_df is not None and subject_rid is not None:
                demog_subject_row = demographics_df[demographics_df['RID'] == subject_rid]
                if not demog_subject_row.empty:
                    demog_data = demog_subject_row.iloc[0]
                    sex_col_demog = next((c for c in ['PTGENDER', 'SEX'] if c in demog_data.index), None)
                    if sex_col_demog and pd.notna(demog_data[sex_col_demog]):
                        sex_val = demog_data[sex_col_demog]
                        sex_map = {1:'male', 2:'female', '1':'male', '2':'female', 'Male':'male', 'Female':'female'}
                        current_metadata['sex'] = sex_map.get(str(sex_val).lower(), sex_map.get(int(sex_val) if isinstance(sex_val, (int,float)) or str(sex_val).isdigit() else None, 'unknown'))
                    current_metadata['PTDOBYY_temp'] = pd.to_numeric(demog_data.get('PTDOBYY'), errors='coerce')
                    if adni_specific_cfg.get('race_source_column') in demog_data:
                        race_val = demog_data[adni_specific_cfg['race_source_column']]
                        race_map = adni_specific_cfg.get('source_race_to_standard_label_map',{})
                        current_metadata['race'] = race_map.get(str(race_val), race_map.get(int(race_val) if pd.notna(race_val) and str(race_val).isdigit() else race_val, 'unknown or not reported'))
                    else: current_metadata['race'] = 'unknown or not reported'
                    if adni_specific_cfg.get('ethnicity_source_column') in demog_data:
                        eth_val = demog_data[adni_specific_cfg['ethnicity_source_column']]
                        eth_map = adni_specific_cfg.get('source_ethnicity_to_standard_label_map',{})
                        current_metadata['is_hispanic_or_latino'] = eth_map.get(str(eth_val), eth_map.get(int(eth_val) if pd.notna(eth_val) and str(eth_val).isdigit() else eth_val, 'unknown or not reported'))
                    else: current_metadata['is_hispanic_or_latino'] = 'unknown or not reported'
                else: logger.warning(f"ADNI: No demographics for RID {subject_rid}")
            exam_date_for_this_visit = None
            if diagnosis_long_df is not None and subject_rid is not None and current_metadata["expression_sample_visit_code"] != "unknown":
                visit_dx_entry = diagnosis_long_df[(diagnosis_long_df['RID'] == subject_rid) & (diagnosis_long_df['VISCODE2'] == current_metadata["expression_sample_visit_code"])]
                if not visit_dx_entry.empty:
                    diag_data = visit_dx_entry.iloc[0]
                    diag_code = diag_data.get('DIAGNOSIS'); exam_date_for_this_visit = diag_data.get('EXAMDATE')
                    diag_map = adni_specific_cfg.get("diagnosis_mapping", {1:"CN", 2:"MCI", 3:"AD", -4:"Unknown"})
                    current_metadata['diagnosis_at_expression_sampling_label'] = diag_map.get(diag_code, "Unknown")
                    current_metadata['diagnosis_at_expression_sampling_code'] = diag_code if pd.notna(diag_code) else pd.NA
                    current_metadata['date_of_expression_sampling'] = str(exam_date_for_this_visit.date()) if pd.notna(exam_date_for_this_visit) and hasattr(exam_date_for_this_visit, 'date') else ""
                else: logger.warning(f"ADNI: No DXSUM for RID {subject_rid}, VISCODE2 {current_metadata['expression_sample_visit_code']}")
            elif diagnosis_long_df is None: logger.debug("ADNI: DXSUM not loaded/RID missing.")
            if pd.notna(current_metadata.get('PTDOBYY_temp')) and pd.notna(exam_date_for_this_visit) and hasattr(exam_date_for_this_visit, 'year'):
                current_metadata['age_at_expression_sampling'] = str(exam_date_for_this_visit.year - int(current_metadata['PTDOBYY_temp']))
            elif demographics_df is not None and subject_rid is not None:
                 demog_row = demographics_df[demographics_df['RID'] == subject_rid]
                 if not demog_row.empty and 'AGE' in demog_row.columns and pd.notna(demog_row['AGE'].iloc[0]):
                    age_baseline = pd.to_numeric(demog_row['AGE'].iloc[0], errors='coerce')
                    if pd.notna(age_baseline): current_metadata['age_at_expression_sampling'] = str(int(age_baseline))
            current_metadata.pop('PTDOBYY_temp', None)
            sample_metadata_list.append(current_metadata)
            processed_expression_files_count += 1
        except Exception as e:
            logger.error(f"Error processing ADNI file {file_path_str} for sample {sample_id}: {e}", exc_info=True)
            if processed_file_path != file_path_str and os.path.exists(processed_file_path):
                os.remove(processed_file_path)

    if processed_expression_files_count == 0:
        logger.error(f"No valid ADNI expression data processed from {input_dir}")
        return None
    logger.info(f"Processed {processed_expression_files_count} ADNI expression files.")

    obs_df_adni = pd.DataFrame(sample_metadata_list)
    if 'sample_id_from_file' not in obs_df_adni.columns or obs_df_adni.empty:
        logger.error(f"{dataset_name_for_log}: 'sample_id_from_file' key missing or obs_df is empty.")
        return None
    obs_df_adni = obs_df_adni.set_index('sample_id_from_file')
    obs_df_adni.index.name = "sample_id"
    final_sample_ids_adni = obs_df_adni.index.tolist()
    if obs_df_adni.index.has_duplicates:
        logger.warning(f"{dataset_name_for_log}: Correcting duplicate sample IDs in obs_df index.")
        obs_df_adni = obs_df_adni[~obs_df_adni.index.duplicated(keep='first')]
        final_sample_ids_adni = obs_df_adni.index.tolist()
    sample_dfs = {sid: df_val for sid, df_val in sample_dfs.items() if sid in final_sample_ids_adni}

    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs for {dataset_name_for_log}")
    gene_id_mapping = {gid: standardize_ensembl_id(gid) for gid in all_gene_ids}
    valid_mapped_values_adni = [str(gid) for gid in gene_id_mapping.values() if gid is not None]
    if not valid_mapped_values_adni:
        logger.error(f"{dataset_name_for_log}: No gene IDs could be meaningfully standardized. Aborting.")
        return None
    unique_std_ids = sorted(list(set(valid_mapped_values_adni)))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs for {dataset_name_for_log}.")

    # Optimized Matrix Creation for ADNI
    logger.info(f"Creating unified expression matrix for {dataset_name_for_log} ({len(unique_std_ids)} genes x {len(final_sample_ids_adni)} samples) using optimized method.")
    all_series = [s_df[sample_id].rename(sample_id) for sample_id, s_df in sample_dfs.items()]
    if not all_series:
        logger.error(f"No data series collected for {dataset_name_for_log}. Aborting matrix creation.")
        unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_adni).astype(np.float32) # Create empty if no series
    else:
        try:
            concat_df = pd.concat(all_series, axis=1, join='outer').astype(np.float32)
            logger.info(f"Concatenated DataFrame shape: {concat_df.shape}")

            std_to_originals_map_adni = {}
            for original_gid, mapped_std_id_val in gene_id_mapping.items():
                if mapped_std_id_val is not None:
                    std_to_originals_map_adni.setdefault(str(mapped_std_id_val), []).append(str(original_gid))

            original_to_std_series_map = {orig_id: std_id_val for std_id_val, orig_ids_list in std_to_originals_map_adni.items() for orig_id in orig_ids_list}
            mappable_original_ids = [idx for idx in concat_df.index if idx in original_to_std_series_map]

            if not mappable_original_ids:
                logger.warning(f"No mappable original gene IDs in concat_df for {dataset_name_for_log}. Resulting matrix might be all zeros or have unexpected shape.")
                unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_adni).astype(np.float32)
            else:
                filtered_concat_df = concat_df.loc[mappable_original_ids]
                std_ids_for_aggregation = filtered_concat_df.index.map(original_to_std_series_map)
                if std_ids_for_aggregation.isna().all():
                    logger.error(f"All standardized IDs for aggregation are NaN for {dataset_name_for_log}. Check mapping.")
                    unified_expr_df = pd.DataFrame(0.0, index=unique_std_ids, columns=final_sample_ids_adni).astype(np.float32)
                else:
                    filtered_concat_df.index = std_ids_for_aggregation
                    unified_expr_df_grouped = filtered_concat_df.groupby(level=0).max()
                    unified_expr_df = unified_expr_df_grouped.reindex(index=unique_std_ids, columns=final_sample_ids_adni, fill_value=0.0).astype(np.float32)
        except Exception as e_concat:
            logger.error(f"Error during optimized matrix creation for {dataset_name_for_log}: {e_concat}", exc_info=True)
            logger.warning("Falling back to previous (slower) matrix creation method for ADNI due to error.")
            # Fallback to original loop method
            expr_matrix = np.full((len(unique_std_ids), len(final_sample_ids_adni)), np.nan, dtype=np.float32)
            sample_id_to_col_idx = {sid: i for i, sid in enumerate(final_sample_ids_adni)}
            std_id_to_idx = {gid: i for i, gid in enumerate(unique_std_ids)} # Re-define if not in scope
            std_to_originals_map_adni = {} # Re-create as it might not be in scope if optimized failed early
            for original_gid, mapped_std_id_val in gene_id_mapping.items():
                if mapped_std_id_val is not None: std_to_originals_map_adni.setdefault(str(mapped_std_id_val), []).append(str(original_gid))

            for std_id, original_ids_for_this_std in std_to_originals_map_adni.items():
                if std_id not in std_id_to_idx: continue
                gene_idx = std_id_to_idx[std_id]
                for sample_id, col_idx in sample_id_to_col_idx.items():
                    sample_df = sample_dfs.get(sample_id)
                    if sample_df is None: continue
                    if sample_id not in sample_df.columns: continue
                    relevant_orig_ids = [orig_id for orig_id in original_ids_for_this_std if orig_id in sample_df.index]
                    if relevant_orig_ids:
                        values = pd.to_numeric(sample_df.loc[relevant_orig_ids, sample_id], errors='coerce')
                        if not values.empty and not values.isna().all(): expr_matrix[gene_idx, col_idx] = values.max()
            expr_matrix = np.nan_to_num(expr_matrix, nan=0.0)
            unified_expr_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=final_sample_ids_adni)

    logger.info(f"Optimized ADNI matrix created. Shape: {unified_expr_df.shape}")
    # End of Optimized Matrix Creation

    var_df_adni = pd.DataFrame(index=unique_std_ids)
    var_df_adni["gene_id"] = var_df_adni.index.astype(str)
    var_df_adni["original_ids"] = var_df_adni.index.map(lambda x: ";".join(std_to_originals_map_adni.get(x, [])))

    # gencode_map passed as parameter
    var_df_adni = add_gencode_annotations(var_df_adni, gencode_map)

    mappings = load_mappings()
    obs_df_adni = standardize_metadata(obs_df_adni, dataset_name_for_log, mappings)
    if obs_df_adni is None:
        logger.error("ADNI: standardize_metadata returned None. Aborting ADNI processing.")
        return None

    obs_df_adni.rename(columns={
        'age_at_expression_sampling': 'age',
        'diagnosis_at_expression_sampling_code': 'diagnosis_code',
        'diagnosis_at_expression_sampling_label': 'diagnosis',
        'date_of_expression_sampling': 'visit_date'
    }, inplace=True)

    for col, dtype, default in [
        ('age', str, ''), ('sex', 'category', 'unknown'),
        ('race', 'category', 'unknown or not reported'),
        ('is_hispanic_or_latino', 'category', 'unknown or not reported'),
        ('self_reported_ethnicity', 'category', 'unknown or not reported'),
        ('self_reported_ethnicity_ontology_term_id', 'category', 'unknown'),
        ('diagnosis', 'category', 'Unknown'),
        ('diagnosis_code', 'Int64', pd.NA),
        ('visit_date', str, ''),
        ('expression_sample_visit_code', str, 'unknown')
    ]:
        if col not in obs_df_adni.columns:
            obs_df_adni[col] = default
        if dtype == 'category':
            current_categories = None
            if col == 'race': current_categories = NIH_RACE_CATEGORIES_LOWER
            elif col == 'is_hispanic_or_latino': current_categories = NIH_ETHNICITY_CATEGORIES_LOWER
            elif col == 'self_reported_ethnicity': current_categories = SRE_BASE_CATEGORIES_LOWER
            elif col == 'self_reported_ethnicity_ontology_term_id': pass # Handled by standardize_metadata

            if current_categories:
                 obs_df_adni[col] = pd.Categorical(
                     obs_df_adni[col].astype(str).fillna(str(default) if pd.notna(default) else ''),
                     categories=current_categories, ordered=False)
            elif col != 'self_reported_ethnicity_ontology_term_id':
                 obs_df_adni[col] = obs_df_adni[col].astype(str).fillna(str(default) if pd.notna(default) else '').astype('category')
        elif dtype == 'Int64':
            obs_df_adni[col] = pd.to_numeric(obs_df_adni[col], errors='coerce').astype('Int64')
        else:
            obs_df_adni[col] = obs_df_adni[col].astype(str).fillna(str(default) if pd.notna(default) else '')

    # ---- Moved Deduplication Logic Here (ADNI) ----
    if obs_df_adni.columns.has_duplicates:
        duplicated_cols_obs = obs_df_adni.columns[obs_df_adni.columns.duplicated(keep=False)].unique().tolist()
        logger.warning(f"ADNI: obs_df_adni has duplicate columns before creating AnnData. Duplicates: {duplicated_cols_obs}")
        obs_df_adni = obs_df_adni.loc[:, ~obs_df_adni.columns.duplicated(keep='first')]
        logger.info(f"ADNI: Deduplicated obs_df_adni columns. New shape: {obs_df_adni.shape}, Columns: {obs_df_adni.columns.tolist()}")

    if var_df_adni.columns.has_duplicates:
        duplicated_cols_var = var_df_adni.columns[var_df_adni.columns.duplicated(keep=False)].unique().tolist()
        logger.warning(f"ADNI: var_df_adni has duplicate columns before creating AnnData. Duplicates: {duplicated_cols_var}")
        var_df_adni = var_df_adni.loc[:, ~var_df_adni.columns.duplicated(keep='first')]
        logger.info(f"ADNI: Deduplicated var_df_adni columns. New shape: {var_df_adni.shape}, Columns: {var_df_adni.columns.tolist()}")
    # ---- End Moved Deduplication Logic ----

    logger.info(f"ADNI: Columns in obs_df_adni before AnnData creation: {obs_df_adni.columns.tolist()}")
    logger.info(f"ADNI: dtypes in obs_df_adni before AnnData creation:\n{obs_df_adni.dtypes}")
    for col_debug in obs_df_adni.columns:
        if isinstance(obs_df_adni[col_debug], pd.DataFrame):
            logger.error(f"ADNI ERROR (FINAL CHECK): obs_df_adni['{col_debug}'] is still a DataFrame!")
    for col_debug_var in var_df_adni.columns:
        if isinstance(var_df_adni[col_debug_var], pd.DataFrame):
            logger.error(f"ADNI ERROR (FINAL CHECK): var_df_adni['{col_debug_var}'] is still a DataFrame!")

    adata = create_standard_anndata(unified_expr_df.T, obs_df_adni, var_df_adni, adni_specific_cfg, dataset_name_for_log)

    if all_subject_diagnoses_for_uns:
        adata.uns['longitudinal_diagnoses_adni'] = all_subject_diagnoses_for_uns
        logger.info("ADNI: Added full longitudinal diagnosis history to adata.uns['longitudinal_diagnoses_adni']")

    if adni_specific_cfg:
        adata = apply_dataset_specific_metadata(adata, adni_specific_cfg)

    if save_anndata(adata, output_file):
        logger.info(f"Successfully processed ADNI data: {adata.n_obs} samples and {adata.n_vars} genes")
    else:
        logger.error(f"Failed to save AnnData for ADNI.")
        return None
    return adata


# ====================================================================================
# Main Processing Pipeline
# ====================================================================================
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Multi-Dataset Standardization Pipeline (Stage 1)")
    parser.add_argument("--encode-dir", help="Directory containing ENCODE cell line TPM files")
    parser.add_argument("--encode-entex-dir", help="Directory containing ENCODE ENTEx TPM files (can be same as encode-dir if co-located)")
    parser.add_argument("--gtex-file", help="Path to GTEx expression file (GCT format)")
    parser.add_argument("--mage-dir", help="Directory containing MAGE expression files")
    parser.add_argument(
        "--adni-dir", help="Directory containing ADNI sample directories with expression files"
    )
    parser.add_argument(
        "--metadata-dir",
        default=str(DEFAULT_METADATA_DIR),
        help=f"Directory containing metadata JSON files (default: {DEFAULT_METADATA_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"Output directory for standardized files (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument("--entex-metadata-file", help="Path to ENTEx metadata JSON file (can override default path if entex-dir is used)")
    parser.add_argument("--adni-demographics-file", help="Path to ADNI patient demographics CSV (e.g., PTDEMOG_*.csv)")
    parser.add_argument("--adni-diagnosis-file", help="Path to ADNI diagnosis summary CSV (e.g., DXSUM_*.csv)")
    return parser.parse_args()

def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    start_time = time.time()

    # Load GENCODE mapping once
    logger.info("Loading GENCODE mapping for the pipeline run...")
    gencode_map = load_gencode_mapping()
    if not gencode_map:
        logger.error("Failed to load GENCODE mapping. This is critical for gene annotations. Exiting.")
        sys.exit(1)
    logger.info("GENCODE mapping loaded successfully.")


    metadata_dir = Path(args.metadata_dir if args.metadata_dir else DEFAULT_METADATA_DIR)
    processed_datasets_summary = []

    # Process ENCODE (which can include ENTEx based on its args)
    if args.encode_dir or args.encode_entex_dir or args.entex_metadata_file:
        encode_output_file = Path(args.output_dir) / "encode_standardized_v1.h5ad"
        logger.info("Starting ENCODE/ENTEx processing...")
        encode_data = process_encode_data(
            args.encode_dir,
            args.encode_entex_dir,
            str(encode_output_file),
            args.entex_metadata_file,
            str(metadata_dir),
            gencode_map=gencode_map # Pass loaded map
        )
        if encode_data:
            processed_label = "ENCODE"
            if args.encode_entex_dir or args.entex_metadata_file:
                processed_label = "ENCODE_ENTEx_Combined" if args.encode_dir else "ENTEx"

            processed_datasets_summary.append((processed_label, encode_data.n_obs, encode_data.n_vars))
            logger.info(
                f"Successfully processed {processed_label} data: {encode_data.n_obs} samples and {encode_data.n_vars} genes"
            )
            if hasattr(encode_data, 'var') and 'mapping_source' in encode_data.var.columns:
                logger.info(f"{processed_label} gene statistics (mapping source):")
                gencode_stats = encode_data.var["mapping_source"].value_counts()
                for source, count in gencode_stats.items():
                    logger.info(f"  {source}: {count} genes ({count/encode_data.n_vars:.1%})")


    # Process GTEx
    if args.gtex_file:
        gtex_output_file = Path(args.output_dir) / "gtex_standardized_v1.h5ad"
        logger.info(f"Processing GTEx data from {args.gtex_file}")
        gtex_data = process_gtex_data(args.gtex_file, str(gtex_output_file), str(metadata_dir), gencode_map=gencode_map) # Pass loaded map
        if gtex_data:
            processed_datasets_summary.append(("GTEx", gtex_data.n_obs, gtex_data.n_vars))
            logger.info(
                f"Successfully processed GTEx data: {gtex_data.n_obs} samples, {gtex_data.n_vars} genes"
            )
            if hasattr(gtex_data, 'var') and 'mapping_source' in gtex_data.var.columns:
                logger.info("GTEx gene statistics (mapping source):")
                gencode_stats = gtex_data.var["mapping_source"].value_counts()
                for source, count in gencode_stats.items():
                    logger.info(f"  {source}: {count} genes ({count/gtex_data.n_vars:.1%})")

    # Process MAGE
    if args.mage_dir:
        mage_output_file = Path(args.output_dir) / "mage_standardized_v1.h5ad"
        logger.info(f"Processing MAGE data from {args.mage_dir}")
        mage_data = process_mage_data(args.mage_dir, str(mage_output_file), str(metadata_dir), gencode_map=gencode_map) # Pass loaded map
        if mage_data:
            processed_datasets_summary.append(("MAGE", mage_data.n_obs, mage_data.n_vars))
            logger.info(
                f"Successfully processed MAGE data: {mage_data.n_obs} samples, {mage_data.n_vars} genes"
            )

    # Process ADNI
    if args.adni_dir:
        adni_output_file = Path(args.output_dir) / "adni_standardized_v1.h5ad"
        logger.info(f"Processing ADNI data from directory: {args.adni_dir}")
        adni_data = process_adni_data(
            args.adni_dir,
            str(adni_output_file),
            adni_demographics_file=args.adni_demographics_file,
            adni_diagnosis_file=args.adni_diagnosis_file,
            dict_file=None,
            metadata_dir=str(metadata_dir),
            gencode_map=gencode_map # Pass loaded map
        )
        if adni_data:
            processed_datasets_summary.append(("ADNI", adni_data.n_obs, adni_data.n_vars))
            logger.info(
                f"Successfully processed ADNI data: {adni_data.n_obs} samples, {adni_data.n_vars} genes"
            )

    logger.info("=== Stage 1 Processing Summary ===")
    if processed_datasets_summary:
        for name, obs, var in processed_datasets_summary:
            logger.info(f"  {name}: {obs} samples, {var} genes.")
    else:
        logger.info("  No datasets were processed in this run based on provided arguments.")

    total_time = time.time() - start_time
    logger.info(f"Total Stage 1 processing time: {total_time:.2f} seconds")


if __name__ == "__main__":
    main()
````

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/standardize_metadata.py`

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

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2/validate_standardized_datasets.py`

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

