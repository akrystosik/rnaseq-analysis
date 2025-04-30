#!/usr/bin/env python3
"""
Multi-Dataset Expression Query Tool

A tool to query gene expression across standardized datasets (ENCODE, GTEx, ENTEx, etc.).
Works with both gene symbols and Ensembl IDs.

Usage:
  # Query by gene symbol across all datasets
  python expression-query.py --symbol TP53
  
  # Query by Ensembl ID
  python expression-query.py --gene ENSG00000141510
  
  # Query a specific dataset
  python expression-query.py --symbol TP53 --dataset encode
  
  # Query ENTEx data
  python expression-query.py --symbol TP53 --dataset entex
  
  # Query a specific tissue/cell line
  python expression-query.py --symbol TP53 --tissue "lung"
  
  # List available tissues
  python expression-query.py --list-tissues
  
  # List available datasets
  python expression-query.py --list-datasets
  
  # Search for genes matching a pattern
  python expression-query.py --search-gene "TP53"
"""

import argparse
import os
import sys
import time
import pandas as pd
import numpy as np
import anndata as ad
import glob
from pathlib import Path

# Define paths
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data")

DEFAULT_DATASETS = {
    'encode': BASE_DIR / 'encode_standardized.h5ad',
    'gtex': BASE_DIR / 'gtex_standardized.h5ad',
    'mage': BASE_DIR / 'mage_standardized.h5ad',
    'adni': BASE_DIR / 'adni_standardized.h5ad',
    'entex': BASE_DIR / 'entex_standardized.h5ad'  # Added ENTEx dataset
}

def find_datasets(base_dir=BASE_DIR):
    """Find all h5ad files in the base directory"""
    h5ad_files = glob.glob(str(base_dir / "*.h5ad"))
    
    datasets = {}
    for file_path in h5ad_files:
        dataset_name = os.path.basename(file_path).replace('_standardized.h5ad', '')
        datasets[dataset_name] = file_path
    
    return datasets

def load_dataset(file_path, dataset_name=None):
    """Load a dataset from an h5ad file"""
    if not os.path.exists(file_path):
        print(f"Error: Dataset file not found: {file_path}")
        return None
        
    name = dataset_name or os.path.basename(file_path).replace('_standardized.h5ad', '')
    print(f"Loading {name} dataset from {file_path}...")
    start_time = time.time()
    
    try:
        adata = ad.read_h5ad(file_path)
        load_time = time.time() - start_time
        print(f"Loaded {name} dataset with {adata.shape[0]} samples and {adata.shape[1]} genes in {load_time:.2f} seconds")
        
        # Check if gene_name column exists
        if 'gene_name' not in adata.var.columns:
            print(f"Warning: Gene name mapping not found in {name} dataset")
        
        return adata
    
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def load_all_datasets(datasets_dict):
    """Load all available datasets"""
    loaded_datasets = {}
    
    for name, path in datasets_dict.items():
        adata = load_dataset(path, name)
        if adata is not None:
            loaded_datasets[name] = adata
    
    if not loaded_datasets:
        print("Error: No datasets could be loaded")
        sys.exit(1)
    
    return loaded_datasets

def find_gene_by_symbol(datasets, symbol):
    """Find a gene by its symbol across all datasets"""
    matches = {}
    
    for dataset_name, adata in datasets.items():
        if 'gene_name' not in adata.var.columns:
            continue
            
        # Look for exact match (case-sensitive)
        exact_matches = adata.var[adata.var['gene_name'] == symbol]
        
        if len(exact_matches) > 0:
            # Use the first match if multiple exist
            matches[dataset_name] = (exact_matches.index[0], 'exact')
            continue
        
        # Try case-insensitive match
        case_insensitive = adata.var[adata.var['gene_name'].str.lower() == symbol.lower()]
        
        if len(case_insensitive) > 0:
            matches[dataset_name] = (case_insensitive.index[0], 'case_insensitive')
    
    if not matches:
        print(f"Error: Gene symbol '{symbol}' not found in any dataset")
        return {}
    
    # Report matches
    print(f"Found '{symbol}' in {len(matches)} datasets:")
    for dataset, (gene_id, match_type) in matches.items():
        match_desc = "exact match" if match_type == 'exact' else "case-insensitive match"
        print(f"  {dataset}: {gene_id} ({match_desc})")
    
    return matches

def find_gene_by_ensembl(datasets, ensembl_id):
    """Find a gene by its Ensembl ID across all datasets"""
    matches = {}
    
    # Extract base ID (without version)
    if '.' in ensembl_id:
        base_id = ensembl_id.split('.')[0]
    else:
        base_id = ensembl_id
    
    for dataset_name, adata in datasets.items():
        # Try exact match
        if ensembl_id in adata.var_names:
            matches[dataset_name] = (ensembl_id, 'exact')
            continue
        
        # Try match with base ID (without version)
        base_matches = [g for g in adata.var_names if g == base_id or g.startswith(base_id + '.')]
        
        if base_matches:
            matches[dataset_name] = (base_matches[0], 'base_id')
    
    if not matches:
        print(f"Error: Gene ID '{ensembl_id}' not found in any dataset")
        return {}
    
    # Report matches
    print(f"Found gene ID '{ensembl_id}' in {len(matches)} datasets:")
    for dataset, (gene_id, match_type) in matches.items():
        match_desc = "exact match" if match_type == 'exact' else "version match"
        print(f"  {dataset}: {gene_id} ({match_desc})")
    
    return matches

def find_tissue_column(adata):
    """Find the tissue column in the dataset"""
    tissue_columns = ['tissue', 'tissue_type', 'tissue_of_origin', 'SMTSD', 'cell_line']
    
    for col in tissue_columns:
        if col in adata.obs.columns:
            return col
    
    return None

def find_subject_column(adata):
    """Find the subject column in the dataset"""
    subject_columns = ['subject_id', 'donor_id', 'donor', 'SUBJID']  # Added 'donor' for ENTEx compatibility
    
    for col in subject_columns:
        if col in adata.obs.columns:
            return col
    
    return None

def query_gene_expression(datasets, gene_matches, tissue=None, subject=None):
    """Query expression for a gene across datasets, optionally filtered by tissue and subject"""
    results = []
    
    for dataset_name, (gene_id, _) in gene_matches.items():
        adata = datasets[dataset_name]
        
        # Get tissue and subject columns
        tissue_col = find_tissue_column(adata)
        subject_col = find_subject_column(adata)
        
        # Create filter mask
        mask = pd.Series(True, index=adata.obs.index)
        
        if tissue and tissue_col:
            # Case-insensitive match on tissue
            tissue_mask = adata.obs[tissue_col].str.lower().str.contains(tissue.lower(), na=False)
            mask = mask & tissue_mask
            
            # If no matches, report and skip
            if not any(tissue_mask):
                print(f"Warning: No samples in '{dataset_name}' match tissue '{tissue}'")
                continue
        
        if subject and subject_col:
            # Match on subject
            subject_mask = adata.obs[subject_col] == subject
            mask = mask & subject_mask
            
            # If no matches, report and skip
            if not any(subject_mask):
                print(f"Warning: No samples in '{dataset_name}' match subject '{subject}'")
                continue
        
        # Apply filter
        filtered_data = adata[mask, gene_id]
        
        if filtered_data.shape[0] == 0:
            print(f"Warning: No matching samples in '{dataset_name}' after filtering")
            continue
        
        # Get expression values
        if isinstance(filtered_data.X, np.ndarray):
            expression_values = filtered_data.X
        else:
            expression_values = filtered_data.X.toarray()
        
        # Get gene information
        gene_name = adata.var.loc[gene_id, 'gene_name'] if 'gene_name' in adata.var.columns else None
        gene_type = adata.var.loc[gene_id, 'gene_type'] if 'gene_type' in adata.var.columns else None
        
        # Get sample information
        if tissue_col:
            tissues = filtered_data.obs[tissue_col].value_counts().to_dict()
        else:
            tissues = {"unknown": filtered_data.shape[0]}
        
        if subject_col:
            subjects = filtered_data.obs[subject_col].value_counts().to_dict()
        else:
            subjects = {"unknown": filtered_data.shape[0]}
        
        # Calculate summary statistics
        result = {
            'dataset': dataset_name,
            'gene_id': gene_id,
            'gene_name': gene_name,
            'gene_type': gene_type,
            'sample_count': int(filtered_data.shape[0]),
            'mean_expression': float(np.mean(expression_values)),
            'median_expression': float(np.median(expression_values)),
            'min_expression': float(np.min(expression_values)),
            'max_expression': float(np.max(expression_values)),
            'tissues': tissues,
            'subjects': subjects
        }
        
        results.append(result)
    
    return results

def print_expression_results(results):
    """Print expression query results"""
    if not results:
        print("No expression results found matching your criteria")
        return
    
    print("\nExpression Results:")
    
    for result in results:
        print(f"\n===== {result['dataset'].upper()} Dataset =====")
        print(f"Gene ID: {result['gene_id']}")
        
        if result['gene_name'] is not None:
            print(f"Gene Symbol: {result['gene_name']}")
        
        if result['gene_type'] is not None:
            print(f"Gene Type: {result['gene_type']}")
        
        print(f"Sample Count: {result['sample_count']}")
        print(f"Expression (TPM):")
        print(f"  Mean: {result['mean_expression']:.4f}")
        print(f"  Median: {result['median_expression']:.4f}")
        print(f"  Min: {result['min_expression']:.4f}")
        print(f"  Max: {result['max_expression']:.4f}")
        
        # Print tissue breakdown
        if len(result['tissues']) > 1:
            print("\nTissue Breakdown:")
            for tissue, count in sorted(result['tissues'].items(), key=lambda x: x[1], reverse=True):
                print(f"  {tissue}: {count} samples")
        
        # Print subject breakdown
        if len(result['subjects']) > 1 and len(result['subjects']) <= 10:
            print("\nSubject Breakdown:")
            for subject, count in sorted(result['subjects'].items(), key=lambda x: x[1], reverse=True):
                print(f"  {subject}: {count} samples")
        elif len(result['subjects']) > 10:
            print(f"\nSubjects: {len(result['subjects'])} unique subjects")

def list_tissues(datasets):
    """List all tissues across datasets"""
    all_tissues = {}
    
    for dataset_name, adata in datasets.items():
        tissue_col = find_tissue_column(adata)
        
        if not tissue_col:
            print(f"Warning: No tissue information found in {dataset_name} dataset")
            continue
        
        # Get tissue counts
        tissues = adata.obs[tissue_col].value_counts().to_dict()
        all_tissues[dataset_name] = tissues
    
    if not all_tissues:
        print("No tissue information found in any dataset")
        return
    
    # Print tissues by dataset
    for dataset, tissues in all_tissues.items():
        print(f"\n===== {dataset.upper()} Dataset =====")
        print(f"Found {len(tissues)} tissue types with {sum(tissues.values())} total samples:")
        
        for tissue, count in sorted(tissues.items(), key=lambda x: x[1], reverse=True):
            print(f"  {tissue}: {count} samples")

def list_datasets(datasets):
    """List all available datasets"""
    print(f"\nFound {len(datasets)} datasets:")
    
    for dataset_name, adata in datasets.items():
        # Get metadata
        sample_count = adata.shape[0]
        gene_count = adata.shape[1]
        
        # Get tissue information
        tissue_col = find_tissue_column(adata)
        if tissue_col:
            tissue_count = len(adata.obs[tissue_col].unique())
        else:
            tissue_count = "unknown"
        
        # Get subject information
        subject_col = find_subject_column(adata)
        if subject_col:
            subject_count = len(adata.obs[subject_col].unique())
        else:
            subject_count = "unknown"
        
        # Print summary
        print(f"\n===== {dataset_name.upper()} Dataset =====")
        print(f"Samples: {sample_count}")
        print(f"Genes: {gene_count}")
        print(f"Tissues/Cell Types: {tissue_count}")
        print(f"Subjects/Donors: {subject_count}")
        
        # Dataset-specific information
        if 'dataset_info' in adata.uns:
            info = adata.uns['dataset_info']
            for key, value in info.items():
                if isinstance(value, (str, int, float, bool)):
                    print(f"{key}: {value}")

def search_genes(datasets, query, max_results=50):
    """Search for genes by ID or symbol across all datasets"""
    all_matches = {}
    
    for dataset_name, adata in datasets.items():
        dataset_matches = []
        
        # Search in gene IDs
        id_matches = [g for g in adata.var_names if query.lower() in g.lower()]
        if id_matches:
            dataset_matches.extend([(g, 'ensembl_id') for g in id_matches[:max_results]])
        
        # Search in gene symbols
        if 'gene_name' in adata.var.columns:
            # Handle possible NaNs in gene_name column
            symbol_mask = adata.var['gene_name'].notna() & adata.var['gene_name'].str.contains(query, case=False, na=False)
            symbol_matches = adata.var[symbol_mask]
            
            if not symbol_matches.empty:
                dataset_matches.extend([
                    (idx, 'gene_symbol') 
                    for idx in symbol_matches.index[:max_results]
                ])
        
        # Store unique matches
        unique_matches = []
        seen_genes = set()
        
        for gene_id, match_type in dataset_matches:
            if gene_id not in seen_genes:
                seen_genes.add(gene_id)
                unique_matches.append((gene_id, match_type))
        
        all_matches[dataset_name] = unique_matches
    
    # Print results
    total_matches = sum(len(matches) for matches in all_matches.values())
    print(f"\nFound {total_matches} matches for '{query}' across {len(all_matches)} datasets.")
    
    for dataset_name, matches in all_matches.items():
        if not matches:
            continue
            
        print(f"\n===== {dataset_name.upper()} Dataset ({len(matches)} matches) =====")
        
        for i, (gene_id, match_type) in enumerate(matches[:20]):
            if i >= 20 and len(matches) > 20:
                print(f"... and {len(matches) - 20} more matches")
                break
                
            adata = datasets[dataset_name]
            gene_name = adata.var.loc[gene_id, 'gene_name'] if 'gene_name' in adata.var.columns and gene_id in adata.var.index else 'N/A'
            gene_type = adata.var.loc[gene_id, 'gene_type'] if 'gene_type' in adata.var.columns and gene_id in adata.var.index else 'N/A'
            
            print(f"- {gene_id} (Symbol: {gene_name}, Type: {gene_type})")

def compare_expression_across_datasets(datasets, gene_matches, tissue=None):
    """
    Compare expression of a gene across different datasets, 
    optionally filtered by tissue type.
    """
    results = []
    
    for dataset_name, (gene_id, _) in gene_matches.items():
        adata = datasets[dataset_name]
        
        # Get tissue column
        tissue_col = find_tissue_column(adata)
        
        # Create filter mask
        mask = pd.Series(True, index=adata.obs.index)
        
        if tissue and tissue_col:
            # Case-insensitive match on tissue
            tissue_mask = adata.obs[tissue_col].str.lower().str.contains(tissue.lower(), na=False)
            mask = mask & tissue_mask
            
            # If no matches, skip
            if not any(tissue_mask):
                continue
        
        # Apply filter
        filtered_data = adata[mask, gene_id]
        
        if filtered_data.shape[0] == 0:
            continue
        
        # Get expression values
        if isinstance(filtered_data.X, np.ndarray):
            expression_values = filtered_data.X
        else:
            expression_values = filtered_data.X.toarray()
        
        # Calculate summary statistics
        result = {
            'dataset': dataset_name,
            'gene_id': gene_id,
            'sample_count': int(filtered_data.shape[0]),
            'mean_expression': float(np.mean(expression_values)),
            'median_expression': float(np.median(expression_values)),
            'std_expression': float(np.std(expression_values))
        }
        
        # Add tissue information if available
        if tissue_col:
            result['tissues'] = filtered_data.obs[tissue_col].value_counts().to_dict()
        
        results.append(result)
    
    return results

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Multi-Dataset Expression Query Tool')
    
    # Query options
    query_group = parser.add_argument_group('Query Options')
    query_group.add_argument('--symbol', help='Query by gene symbol (e.g., TP53)')
    query_group.add_argument('--gene', help='Query by Ensembl gene ID (e.g., ENSG00000141510)')
    query_group.add_argument('--dataset', help='Dataset to query (e.g., encode, gtex, entex, all)')
    query_group.add_argument('--tissue', help='Tissue or cell type (e.g., "lung", "blood")')
    query_group.add_argument('--subject', help='Subject or donor ID')
    
    # List options
    list_group = parser.add_argument_group('List Options')
    list_group.add_argument('--list-tissues', action='store_true', help='List all tissues in the datasets')
    list_group.add_argument('--list-datasets', action='store_true', help='List all available datasets')
    list_group.add_argument('--search-gene', help='Search for genes matching a pattern')
    
    # Compare options
    compare_group = parser.add_argument_group('Compare Options')
    compare_group.add_argument('--compare', action='store_true', 
                               help='Compare gene expression across datasets')
    
    # General options
    parser.add_argument('--data-dir', default=BASE_DIR, help='Directory containing standardized datasets')
    
    args = parser.parse_args()
    
    # Find available datasets
    available_datasets = find_datasets(args.data_dir)
    
    if not available_datasets:
        print(f"Error: No dataset files found in {args.data_dir}")
        sys.exit(1)
    
    print(f"Found {len(available_datasets)} available datasets: {', '.join(available_datasets.keys())}")
    
    # Determine which datasets to load
    if args.dataset and args.dataset.lower() != 'all':
        if args.dataset.lower() not in available_datasets:
            print(f"Error: Dataset '{args.dataset}' not found. Available: {', '.join(available_datasets.keys())}")
            sys.exit(1)
        
        datasets_to_load = {args.dataset.lower(): available_datasets[args.dataset.lower()]}
    else:
        datasets_to_load = available_datasets
    
    # Load datasets
    datasets = load_all_datasets(datasets_to_load)
    
    # Handle list options
    if args.list_tissues:
        list_tissues(datasets)
        return
    
    if args.list_datasets:
        list_datasets(datasets)
        return
    
    if args.search_gene:
        search_genes(datasets, args.search_gene)
        return
    
    # Handle query options
    if args.symbol or args.gene:
        # Find the gene
        if args.symbol:
            gene_matches = find_gene_by_symbol(datasets, args.symbol)
        else:
            gene_matches = find_gene_by_ensembl(datasets, args.gene)
        
        if not gene_matches:
            return
        
        # If comparing across datasets
        if args.compare:
            results = compare_expression_across_datasets(datasets, gene_matches, args.tissue)
            
            if not results:
                print("No comparable expression data found across datasets")
                return
                
            print("\nExpression Comparison Across Datasets:")
            for result in sorted(results, key=lambda x: x['mean_expression'], reverse=True):
                print(f"{result['dataset']} ({result['sample_count']} samples): Mean TPM {result['mean_expression']:.4f} ± {result['std_expression']:.4f}")
            
            return
        
        # Query expression
        start_time = time.time()
        results = query_gene_expression(datasets, gene_matches, args.tissue, args.subject)
        query_time = time.time() - start_time
        
        print_expression_results(results)
        print(f"\nQuery completed in {query_time:.4f} seconds")
        
        return
    
    # If we get here, no valid options were provided
    parser.print_help()

if __name__ == "__main__":
    main()