#!/usr/bin/env python3
"""
Check current status of age and sex fields in datasets.
"""

import scanpy as sc
import pandas as pd

def check_dataset_fields():
    data_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137'

    dataset_files = {
        'ADNI': 'adni_standardized_preprocessed.h5ad',
        'ENCODE': 'encode_standardized_preprocessed.h5ad',
        'GTEx': 'gtex_standardized_preprocessed.h5ad',
        'MAGE': 'mage_standardized_preprocessed.h5ad'
    }

    for name, filename in dataset_files.items():
        file_path = f'{data_dir}/{filename}'
        adata = sc.read_h5ad(file_path)
        
        print(f'--- {name} Dataset Check ---')
        print(f'Shape: {adata.shape}')
        
        # Check age field
        if 'age' in adata.obs.columns:
            age_values = adata.obs['age']
            non_empty_age = age_values[age_values != ''].dropna()
            print(f'Age field: {len(non_empty_age):,}/{len(adata.obs):,} with valid data')
            if len(non_empty_age) > 0:
                print(f'  Examples: {list(non_empty_age.unique())[:3]}')
            else:
                print('  All age values are empty or NaN')
        else:
            print('No age field found')
        
        # Check sex field  
        if 'sex' in adata.obs.columns:
            sex_counts = adata.obs['sex'].value_counts()
            print(f'Sex distribution: {dict(sex_counts)}')
        else:
            print('No sex field found')
        
        print()

if __name__ == '__main__':
    check_dataset_fields()