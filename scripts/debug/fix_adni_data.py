#!/usr/bin/env python3
"""
ADNI Data Processor

This script processes ADNI gene expression data files with proper handling
of escaped tab characters and creates a standardized AnnData object.
"""

import os
import sys
import glob
import logging
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('adni_processor')

def fix_file_delimiters(file_path):
    """Fix escaped tab delimiters in file."""
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

def parse_adni_file(file_path):
    """Parse ADNI expression file with proper handling of delimiters."""
    try:
        # Fix delimiters if needed
        fixed_path = fix_file_delimiters(file_path)
        
        # Try reading with tab delimiter
        try:
            df = pd.read_csv(fixed_path, sep='\t')
            
            # Check if we have at least the expected columns
            if len(df.columns) >= 2:
                # Identify gene_id and TPM columns
                gene_id_col = None
                tpm_col = None
                
                for col in df.columns:
                    if col.lower() == 'gene_id':
                        gene_id_col = col
                    elif col.lower() == 'tpm':
                        tpm_col = col
                
                if gene_id_col is not None and tpm_col is not None:
                    return df[[gene_id_col, tpm_col]].rename(
                        columns={gene_id_col: 'gene_id', tpm_col: 'TPM'})
            
            # If columns aren't recognized, check if headers contain tabs
            if len(df.columns) == 1 and '\t' in df.columns[0]:
                # Split the header
                headers = df.columns[0].split('\t')
                
                # Split the content of the first column
                new_df = df[df.columns[0]].str.split('\t', expand=True)
                new_df.columns = headers
                
                # Find gene_id and TPM columns
                gene_id_col = None
                tpm_col = None
                
                for col in new_df.columns:
                    if col.lower() == 'gene_id':
                        gene_id_col = col
                    elif col.lower() == 'tpm':
                        tpm_col = col
                
                if gene_id_col is not None and tpm_col is not None:
                    return new_df[[gene_id_col, tpm_col]].rename(
                        columns={gene_id_col: 'gene_id', tpm_col: 'TPM'})
        
        except Exception as e:
            logger.warning(f"Failed to parse with tab delimiter: {e}")
            
            # Try with comma delimiter
            df = pd.read_csv(fixed_path)
            
            # Try to find gene_id and TPM columns
            gene_id_col = None
            tpm_col = None
            
            for col in df.columns:
                if col.lower() == 'gene_id':
                    gene_id_col = col
                elif col.lower() == 'tpm':
                    tpm_col = col
            
            if gene_id_col is not None and tpm_col is not None:
                return df[[gene_id_col, tpm_col]].rename(
                    columns={gene_id_col: 'gene_id', tpm_col: 'TPM'})
        
        # Clean up temporary file if it was created
        if fixed_path != file_path and os.path.exists(fixed_path):
            os.remove(fixed_path)
        
        logger.warning(f"Could not identify gene_id and TPM columns in {file_path}")
        logger.warning(f"Available columns: {df.columns.tolist()}")
        return None
    
    except Exception as e:
        logger.error(f"Error processing {file_path}: {e}")
        return None

def process_adni_dataset(input_dir, output_file):
    """Process all ADNI samples and create standardized AnnData."""
    logger.info(f"Processing ADNI data from {input_dir}")
    
    # Find all subject directories
    subject_dirs = glob.glob(os.path.join(input_dir, "*_S_*"))
    
    if not subject_dirs:
        logger.error(f"No subject directories found in {input_dir}")
        return None
    
    logger.info(f"Found {len(subject_dirs)} subject directories")
    
    # Process each subject
    sample_dfs = {}
    sample_metadata = {}
    all_gene_ids = set()
    processed_count = 0
    
    for subject_dir in subject_dirs:
        subject_id = os.path.basename(subject_dir)
        logger.info(f"Processing subject {subject_id}")
        
        # Find all CSV files for this subject
        csv_files = glob.glob(os.path.join(subject_dir, "*.csv"))
        
        if not csv_files:
            logger.warning(f"No CSV files found for subject {subject_id}")
            continue
        
        # Prioritize gencode_v24 files if available
        gencode_files = [f for f in csv_files if "gencode_v24" in f]
        files_to_process = gencode_files if gencode_files else csv_files
        
        # Process files for this subject
        for file_path in files_to_process:
            file_name = os.path.basename(file_path)
            sample_id = f"{subject_id}_{os.path.splitext(file_name)[0]}"
            
            # Parse the file
            expr_df = parse_adni_file(file_path)
            
            if expr_df is not None:
                # Store the expression data
                sample_dfs[sample_id] = expr_df
                
                # Collect gene IDs
                all_gene_ids.update(expr_df['gene_id'])
                
                # Create metadata
                metadata = {
                    'sample_id': sample_id,
                    'subject_id': subject_id,
                    'dataset': 'ADNI',
                    'data_type': 'RNA-seq',
                    'expression_unit': 'TPM',
                    'tissue': 'blood',  # ADNI samples are from blood
                    'processing': 'gencode_v24' if 'gencode_v24' in file_name else 'unknown',
                }
                
                # Store metadata
                sample_metadata[sample_id] = metadata
                processed_count += 1
                
                # Only process one file per subject (the first successful one)
                break
    
    if not sample_dfs:
        logger.error("No valid expression data files found")
        return None
    
    logger.info(f"Processed {processed_count} samples from {len(subject_dirs)} subjects")
    
    # Standardize gene IDs by removing version numbers
    logger.info(f"Standardizing {len(all_gene_ids)} gene IDs")
    gene_id_mapping = {}
    for gene_id in all_gene_ids:
        # Remove version number if present
        if '.' in gene_id:
            std_id = gene_id.split('.')[0]
        else:
            std_id = gene_id
        gene_id_mapping[gene_id] = std_id
    
    # Create unified expression matrix
    unique_std_ids = sorted(set(gene_id_mapping.values()))
    logger.info(f"Found {len(unique_std_ids)} unique standardized gene IDs")
    
    # Create mappings for faster lookups
    std_id_to_idx = {id: idx for idx, id in enumerate(unique_std_ids)}
    sample_ids = list(sample_dfs.keys())
    
    # Create matrix
    expr_matrix = np.zeros((len(unique_std_ids), len(sample_ids)), dtype=np.float32)
    
    # Fill matrix
    for sample_idx, sample_id in enumerate(sample_ids):
        expr_df = sample_dfs[sample_id]
        
        for _, row in expr_df.iterrows():
            gene_id = row['gene_id']
            expr_val = row['TPM']
            
            if gene_id in gene_id_mapping:
                std_id = gene_id_mapping[gene_id]
                gene_idx = std_id_to_idx[std_id]
                expr_matrix[gene_idx, sample_idx] = expr_val
    
    # Create DataFrame
    data_df = pd.DataFrame(expr_matrix, index=unique_std_ids, columns=sample_ids)
    
    # Create observation metadata
    obs_df = pd.DataFrame([sample_metadata[sample_id] for sample_id in sample_ids])
    obs_df.index = sample_ids
    
    # Create variable metadata
    var_df = pd.DataFrame(index=unique_std_ids)
    var_df['gene_id'] = var_df.index
    var_df['original_gene_id'] = [
        next((g for g in all_gene_ids if gene_id_mapping[g] == std_id), std_id)
        for std_id in unique_std_ids
    ]
    
    # Create AnnData object
    adata = ad.AnnData(X=data_df.T, obs=obs_df, var=var_df)
    
    # Add dataset info
    adata.uns['dataset_info'] = {
        'source': 'ADNI',
        'gencode_version': '24',
        'data_type': 'RNA-seq',
        'expression_unit': 'TPM',
        'samples': len(obs_df),
        'genes': len(var_df),
        'subject_count': len(set(obs_df['subject_id'])),
    }
    
    # Add reference genome info
    adata.uns['reference_genome'] = 'hg38'
    adata.uns['gencode_version'] = '24'
    
    # Save the standardized AnnData
    try:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Save the AnnData
        adata.write(output_file)
        logger.info(f"Saved standardized ADNI dataset to {output_file}")
        
        return adata
    
    except Exception as e:
        logger.error(f"Error saving AnnData: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Process ADNI gene expression data")
    parser.add_argument("--input-dir", required=True, help="Directory containing ADNI data")
    parser.add_argument("--output-file", required=True, help="Output file path for standardized AnnData")
    
    args = parser.parse_args()
    
    process_adni_dataset(args.input_dir, args.output_file)

if __name__ == "__main__":
    main()