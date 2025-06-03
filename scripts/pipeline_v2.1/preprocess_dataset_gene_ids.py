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

        # 5. Handle complex versioned Ensembl IDs (e.g., "ENSG00000000003.14;ENSG00000000003.10")
        if not target_id.startswith(('ENSG', 'ENTREZ:', 'gSpikein')) and 'ENSG' in original_id_source:
            # Try to extract clean ENSG ID from complex format
            ensg_parts = []
            if ';' in original_id_source:
                # Split by semicolon and extract ENSG IDs
                parts = original_id_source.split(';')
                for part in parts:
                    if 'ENSG' in part:
                        clean_ensg = part.split('.')[0].strip()  # Remove version
                        if clean_ensg.startswith('ENSG'):
                            ensg_parts.append(clean_ensg)
            elif 'ENSG' in original_id_source:
                # Single ENSG, possibly with version
                clean_ensg = original_id_source.split('.')[0].strip()
                if clean_ensg.startswith('ENSG'):
                    ensg_parts.append(clean_ensg)
            
            # Use the first valid ENSG ID found
            if ensg_parts:
                candidate_ensg = ensg_parts[0]
                if candidate_ensg in ensembl_to_info:
                    target_id = candidate_ensg
                    mapping_source = 'ensg_extraction'
                    mapped_count += 1
                else:
                    # ENSG not in reference, but still use it
                    target_id = candidate_ensg
                    mapping_source = 'ensg_extraction_unref'
                    mapped_count += 1
            else:
                target_id = input_var_name
                mapping_source = 'unmapped'
                unmapped_count += 1
        
        # 6. If target_id is still empty or not standard, mark as unmapped
        elif not target_id.startswith(('ENSG', 'ENTREZ:', 'gSpikein')):
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
        else: # Entrez ID, extracted ENSG, or Unmapped
            var_columns['gene_name'].append('')
            var_columns['gene_type'].append('unknown')
            var_columns['chromosome'].append('')
            if target_id.startswith('ENTREZ:'):
                var_columns['mapping_confidence'].append('low')
            elif target_id.startswith('ENSG'):
                var_columns['mapping_confidence'].append('medium')
            else:
                var_columns['mapping_confidence'].append('none')

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