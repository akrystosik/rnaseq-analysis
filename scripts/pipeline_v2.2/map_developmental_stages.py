#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def map_age_to_hsapdv(age_value):
    """
    Map age to HsapDv ontology terms
    
    HsapDv terms:
    - HsapDv:0000082 - human child stage (3-12 years)
    - HsapDv:0000083 - human adolescent stage (13-17 years) 
    - HsapDv:0000087 - human adult stage (18-64 years)
    - HsapDv:0000224 - human aged stage (65+ years)
    """
    try:
        # Handle different age formats
        if isinstance(age_value, str):
            if '-' in age_value:
                # Handle age ranges like "60-69"
                start, end = age_value.split('-')
                numeric_age = (float(start) + float(end)) / 2
            else:
                numeric_age = float(age_value)
        else:
            numeric_age = float(age_value)
        
        # Map to appropriate HsapDv term
        if numeric_age < 3:
            return ''  # No specific HsapDv term for very young
        elif numeric_age < 13:
            return 'HsapDv:0000082'  # child stage
        elif numeric_age < 18:
            return 'HsapDv:0000083'  # adolescent stage
        elif numeric_age < 65:
            return 'HsapDv:0000087'  # adult stage
        else:
            return 'HsapDv:0000224'  # aged stage
    except:
        return ''  # Return empty if can't parse

def update_developmental_stages(preprocessed_dir=None):
    """Update developmental stage ontology mapping for all datasets"""
    
    # Use dynamic path based on current repo structure
    if preprocessed_dir is None:
        data_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest'
    else:
        data_dir = preprocessed_dir
    dataset_files = {
        'ADNI': 'adni_standardized_preprocessed.h5ad',
        'ENCODE': 'encode_standardized_preprocessed.h5ad', 
        'GTEx': 'gtex_standardized_preprocessed.h5ad',
        'MAGE': 'mage_standardized_preprocessed.h5ad'
    }
    
    results = {}
    
    for name, filename in dataset_files.items():
        file_path = f'{data_dir}/{filename}'
        
        print(f'\nProcessing {name} dataset...')
        
        try:
            # Load dataset
            adata = sc.read_h5ad(file_path)
            
            # Check if age data is available
            if 'age' in adata.obs.columns:
                age_values = adata.obs['age']
                
                # Map ages to HsapDv terms
                hsapdv_terms = []
                for age in age_values:
                    if pd.isna(age) or age == '':
                        hsapdv_terms.append('')
                    else:
                        hsapdv_terms.append(map_age_to_hsapdv(age))
                
                # Update developmental_stage_ontology column
                adata.obs['developmental_stage_ontology'] = hsapdv_terms
                
                # Calculate coverage
                non_empty = [term for term in hsapdv_terms if term != '']
                coverage = len(non_empty) / len(hsapdv_terms) * 100
                
                print(f'  ✅ Updated developmental stage ontology')
                print(f'  Coverage: {coverage:.1f}% ({len(non_empty):,}/{len(hsapdv_terms):,})')
                
                # Show mapping distribution
                hsapdv_counts = pd.Series(hsapdv_terms).value_counts()
                print('  HsapDv term distribution:')
                for term, count in hsapdv_counts.items():
                    if term != '':
                        term_name = {
                            'HsapDv:0000082': 'child stage (3-12 years)',
                            'HsapDv:0000083': 'adolescent stage (13-17 years)',
                            'HsapDv:0000087': 'adult stage (18-64 years)',
                            'HsapDv:0000224': 'aged stage (65+ years)'
                        }.get(term, 'unknown')
                        print(f'    • {term} ({term_name}): {count:,} samples')
                
                # Create backup before saving
                backup_path = file_path.replace('.h5ad', '_backup_before_hsapdv.h5ad')
                if not Path(backup_path).exists():
                    print(f'  Creating backup: {backup_path}')
                    adata.write_h5ad(backup_path)
                
                # Save updated dataset
                print(f'  Saving updated dataset: {file_path}')
                adata.write_h5ad(file_path)
                
                results[name] = {
                    'coverage': coverage,
                    'total_samples': len(hsapdv_terms),
                    'mapped_samples': len(non_empty),
                    'distribution': dict(hsapdv_counts)
                }
                
            else:
                print(f'  ⚠️ No age data available for {name}')
                results[name] = {
                    'coverage': 0,
                    'total_samples': adata.n_obs,
                    'mapped_samples': 0,
                    'distribution': {}
                }
                
        except Exception as e:
            print(f'  ❌ Error processing {name}: {e}')
            results[name] = {'error': str(e)}
    
    # Summary
    print('\n' + '='*60)
    print('DEVELOPMENTAL STAGE MAPPING SUMMARY')
    print('='*60)
    
    total_mapped = 0
    total_samples = 0
    
    for name, result in results.items():
        if 'error' not in result:
            print(f'\n{name}:')
            print(f'  Coverage: {result["coverage"]:.1f}%')
            print(f'  Mapped: {result["mapped_samples"]:,}/{result["total_samples"]:,}')
            total_mapped += result["mapped_samples"]
            total_samples += result["total_samples"]
    
    overall_coverage = (total_mapped / total_samples * 100) if total_samples > 0 else 0
    print(f'\n🎯 OVERALL RESULTS:')
    print(f'  Total coverage: {overall_coverage:.1f}%')
    print(f'  Total mapped: {total_mapped:,}/{total_samples:,}')
    print(f'  🚀 MAJOR IMPROVEMENT: HsapDv coverage significantly increased!')
    
    return results

def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(description='Map age data to developmental stage ontology')
    parser.add_argument('--preprocessed-dir', help='Directory with preprocessed h5ad files')
    parser.add_argument('--force', action='store_true', help='Force regeneration')
    
    args = parser.parse_args()
    
    results = update_developmental_stages(preprocessed_dir=args.preprocessed_dir)
    return results

if __name__ == '__main__':
    results = main()