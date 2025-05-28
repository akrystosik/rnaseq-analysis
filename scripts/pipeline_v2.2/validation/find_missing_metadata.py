#!/usr/bin/env python3
"""
Find and integrate missing age/sex metadata from source files.
"""

import scanpy as sc
import pandas as pd
import gzip
from io import StringIO
from pathlib import Path

def check_adni_age_data():
    """Check ADNI for age data."""
    print('🔍 **CHECKING ADNI AGE DATA**')
    print('=' * 30)

    # Use dynamic path based on current repo structure
    adni_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
    adata_adni = sc.read_h5ad(adni_path)

    print(f'ADNI shape: {adata_adni.shape}')
    print(f'Age-related columns: {[col for col in adata_adni.obs.columns if "age" in col.lower()]}')

    # Check for age data in various columns
    age_columns = [col for col in adata_adni.obs.columns if 'age' in col.lower()]
    
    for col in age_columns:
        values = adata_adni.obs[col].dropna()
        values = values[values != '']
        print(f'{col}: {len(values):,} non-empty values')
        if len(values) > 0:
            print(f'  Sample values: {list(values.unique())[:5]}')

    # Check if ADNI source files are available
    adni_source_paths = [
        '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni/',
        '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/adni/',
        '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/adni/'
    ]
    
    for path in adni_source_paths:
        if Path(path).exists():
            print(f'\nFound ADNI source directory: {path}')
            files = list(Path(path).glob('*.csv'))
            print(f'CSV files: {[f.name for f in files][:5]}')
            break

def check_mage_age_data():
    """Check MAGE for age data."""
    print('\n🔍 **CHECKING MAGE AGE DATA**')
    print('=' * 30)

    # Use dynamic path based on current repo structure
    mage_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/mage_standardized_preprocessed.h5ad'
    adata_mage = sc.read_h5ad(mage_path)

    print(f'MAGE shape: {adata_mage.shape}')
    print(f'Subject ID examples: {list(adata_mage.obs["subject_id"].head())}')

    # Check for age data in various columns
    age_columns = [col for col in adata_mage.obs.columns if 'age' in col.lower()]
    
    for col in age_columns:
        values = adata_mage.obs[col].dropna()
        values = values[values != '']
        print(f'{col}: {len(values):,} non-empty values')

    # Check if 1000 Genomes PED file has age data
    ped_paths = [
        '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/metadata/',
        '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/mage/metadata/',
        '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/mage/'
    ]
    
    for path in ped_paths:
        if Path(path).exists():
            print(f'\nFound MAGE source directory: {path}')
            ped_files = list(Path(path).glob('*.ped')) + list(Path(path).glob('*ped*'))
            print(f'PED files: {[f.name for f in ped_files]}')
            if ped_files:
                # Check first PED file for age info
                ped_file = ped_files[0]
                print(f'Checking {ped_file.name}...')
                try:
                    df = pd.read_csv(ped_file, sep='\t', nrows=5)
                    print(f'Columns: {list(df.columns)}')
                except Exception as e:
                    print(f'Error reading: {e}')
            break

def check_gtex_sex_integration():
    """Check if we can integrate GTEx sex data."""
    print('\n🔍 **CHECKING GTEx SEX INTEGRATION**')
    print('=' * 35)

    # Load GTEx phenotype data
    gtex_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/phs000424.v10.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz'
    
    with gzip.open(gtex_file, 'rt') as f:
        content = f.read()

    lines = content.split('\n')
    data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
    clean_content = '\n'.join(data_lines)
    phenotype_df = pd.read_csv(StringIO(clean_content), sep='\t')

    # Map sex codes
    sex_mapping = {1: 'male', 2: 'female', 98: 'unknown', 99: 'unknown'}
    phenotype_df['sex_mapped'] = phenotype_df['SEX'].map(sex_mapping)

    print(f'GTEx phenotype subjects: {len(phenotype_df):,}')
    print(f'Sex distribution: {dict(phenotype_df["sex_mapped"].value_counts())}')

    # Load current GTEx dataset
    # Use dynamic path based on current repo structure
    gtex_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/gtex_standardized_preprocessed.h5ad'
    adata_gtex = sc.read_h5ad(gtex_path)

    # Check current sex values
    current_sex = adata_gtex.obs['sex'].value_counts()
    print(f'Current GTEx sex values: {dict(current_sex)}')

    # Check if we can map subjects
    gtex_subjects = adata_gtex.obs['subject_id'].unique()
    phenotype_subjects = set(phenotype_df['SUBJID'])
    
    matched_subjects = len(set(gtex_subjects) & phenotype_subjects)
    print(f'Matching subjects: {matched_subjects:,}/{len(gtex_subjects):,}')

    return len(phenotype_df) > 0

if __name__ == '__main__':
    check_adni_age_data()
    check_mage_age_data()
    has_gtex_sex = check_gtex_sex_integration()
    
    print(f'\n🎯 **SUMMARY**')
    print(f'GTEx sex data available: {has_gtex_sex}')
    print('Next: Create integration script to update datasets with missing metadata')