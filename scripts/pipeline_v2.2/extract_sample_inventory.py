#!/usr/bin/env python3
"""
Extract comprehensive sample and subject ID inventory from RNA-seq pipeline datasets.
"""

import pandas as pd
import anndata as ad
import numpy as np
import os
import re
from collections import defaultdict
import json

def extract_subject_from_sample(sample_id, dataset):
    """
    Extract subject ID from sample ID based on dataset-specific patterns.
    """
    if dataset == 'GTEx':
        # GTEx pattern: GTEX-XXXXX-YYYY-ZZ-AA-BB where GTEX-XXXXX is subject
        match = re.match(r'(GTEX-[A-Z0-9]+)', str(sample_id))
        return match.group(1) if match else sample_id
    
    elif dataset == 'ENCODE':
        # ENCODE typically uses cell line names or donor IDs
        # Often sample == subject for cell lines
        return sample_id
    
    elif dataset == 'MAGE':
        # 1000 Genomes: HG/NA prefix patterns
        # Sample IDs like HG00096, NA18486 - subject is the same
        return sample_id
    
    elif dataset == 'ADNI':
        # ADNI subject IDs are typically consistent
        return sample_id
    
    else:
        return sample_id

def analyze_h5ad_file(file_path, dataset_name):
    """
    Extract sample and subject information from H5AD file.
    """
    print(f"Analyzing {dataset_name}: {file_path}")
    
    try:
        adata = ad.read_h5ad(file_path)
        
        # Get basic info
        n_samples = adata.n_obs
        n_genes = adata.n_vars
        
        # Extract sample IDs (observation names)
        sample_ids = adata.obs_names.tolist()
        
        # Try to get subject IDs from metadata
        subject_ids = []
        if 'subject_id' in adata.obs.columns:
            subject_ids = adata.obs['subject_id'].tolist()
        elif 'donor_id' in adata.obs.columns:
            subject_ids = adata.obs['donor_id'].tolist()
        elif 'individual_id' in adata.obs.columns:
            subject_ids = adata.obs['individual_id'].tolist()
        else:
            # Extract from sample IDs using pattern matching
            subject_ids = [extract_subject_from_sample(sid, dataset_name) for sid in sample_ids]
        
        # Get data type if available
        data_types = []
        if 'data_type' in adata.obs.columns:
            data_types = adata.obs['data_type'].tolist()
        else:
            data_types = ['RNA-seq'] * len(sample_ids)  # Default assumption
        
        # Get other relevant metadata columns
        metadata_cols = list(adata.obs.columns)
        
        results = {
            'dataset': dataset_name,
            'file_path': file_path,
            'n_samples': n_samples,
            'n_genes': n_genes,
            'sample_ids': sample_ids,
            'subject_ids': subject_ids,
            'data_types': data_types,
            'metadata_columns': metadata_cols,
            'unique_subjects': len(set(subject_ids)),
        }
        
        # Analyze ID patterns
        sample_patterns = analyze_id_patterns(sample_ids, f"{dataset_name}_samples")
        subject_patterns = analyze_id_patterns(subject_ids, f"{dataset_name}_subjects")
        
        results['sample_patterns'] = sample_patterns
        results['subject_patterns'] = subject_patterns
        
        return results
        
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def analyze_id_patterns(id_list, id_type):
    """
    Analyze patterns in ID lists.
    """
    patterns = defaultdict(int)
    
    for id_val in id_list:
        id_str = str(id_val)
        
        # Common patterns
        if id_str.startswith('GTEX-'):
            patterns['GTEX_format'] += 1
        elif id_str.startswith('HG'):
            patterns['HG_format'] += 1
        elif id_str.startswith('NA'):
            patterns['NA_format'] += 1
        elif re.match(r'^[A-Z0-9]+-[A-Z0-9]+', id_str):
            patterns['hyphenated'] += 1
        elif re.match(r'^\d+$', id_str):
            patterns['numeric_only'] += 1
        elif re.match(r'^[A-Za-z]+\d+', id_str):
            patterns['alpha_numeric'] += 1
        else:
            patterns['other'] += 1
    
    return dict(patterns)

def main():
    """
    Main function to extract sample inventory from pipeline datasets.
    """
    # Path to the most recent run
    base_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250528_235853"
    
    # Dataset files to analyze
    dataset_files = {
        'ENCODE': f"{base_path}/encode_standardized_v2.h5ad",
        'GTEx': f"{base_path}/gtex_standardized_v2.h5ad", 
        'MAGE': f"{base_path}/mage_standardized_v2.h5ad",
        'ADNI': f"{base_path}/adni_standardized_v2.h5ad"
    }
    
    all_results = []
    inventory_data = []
    
    print("Starting sample inventory extraction...")
    
    for dataset, file_path in dataset_files.items():
        if os.path.exists(file_path):
            result = analyze_h5ad_file(file_path, dataset)
            if result:
                all_results.append(result)
                
                # Add to inventory
                for i, sample_id in enumerate(result['sample_ids']):
                    inventory_data.append({
                        'sample_id': sample_id,
                        'subject_id': result['subject_ids'][i],
                        'dataset': dataset,
                        'data_type': result['data_types'][i]
                    })
        else:
            print(f"File not found: {file_path}")
    
    # Create summary statistics
    summary = {}
    for result in all_results:
        dataset = result['dataset']
        summary[dataset] = {
            'total_samples': result['n_samples'],
            'unique_subjects': result['unique_subjects'],
            'n_genes': result['n_genes'],
            'sample_patterns': result['sample_patterns'],
            'subject_patterns': result['subject_patterns'],
            'metadata_columns': result['metadata_columns']
        }
    
    # Save inventory CSV
    df_inventory = pd.DataFrame(inventory_data)
    output_csv = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/sample_subject_inventory.csv"
    df_inventory.to_csv(output_csv, index=False)
    print(f"Saved inventory to: {output_csv}")
    
    # Save summary JSON
    output_json = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/sample_inventory_summary.json"
    with open(output_json, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Saved summary to: {output_json}")
    
    # Print summary
    print("\n" + "="*60)
    print("SAMPLE INVENTORY SUMMARY")
    print("="*60)
    
    total_samples = 0
    total_subjects = 0
    
    for dataset, stats in summary.items():
        print(f"\n{dataset}:")
        print(f"  Total samples: {stats['total_samples']:,}")
        print(f"  Unique subjects: {stats['unique_subjects']:,}")
        print(f"  Genes: {stats['n_genes']:,}")
        print(f"  Sample ID patterns: {stats['sample_patterns']}")
        print(f"  Subject ID patterns: {stats['subject_patterns']}")
        
        total_samples += stats['total_samples']
        total_subjects += stats['unique_subjects']
    
    print(f"\nOVERALL TOTALS:")
    print(f"  Total samples across all datasets: {total_samples:,}")
    print(f"  Total unique subjects: {total_subjects:,}")
    
    # Show sample data examples
    print(f"\nSAMPLE DATA PREVIEW:")
    print(df_inventory.head(20).to_string(index=False))
    
    print(f"\nDATA TYPE DISTRIBUTION:")
    print(df_inventory['data_type'].value_counts())
    
    print(f"\nSAMPLES PER DATASET:")
    print(df_inventory['dataset'].value_counts())

if __name__ == "__main__":
    main()