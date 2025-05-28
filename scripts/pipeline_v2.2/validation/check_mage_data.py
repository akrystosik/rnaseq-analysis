#!/usr/bin/env python3

import scanpy as sc
import pandas as pd

# Load MAGE dataset
adata = sc.read_h5ad('/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137/mage_standardized_preprocessed.h5ad')

print('MAGE Age Data Analysis:')
print('=' * 30)

# Check age field
if 'age' in adata.obs.columns:
    age_values = adata.obs['age'].dropna()
    non_empty_age = age_values[age_values != '']
    print(f'Age field present: Yes')
    print(f'Non-empty age values: {len(non_empty_age)} out of {len(adata.obs)} ({len(non_empty_age)/len(adata.obs)*100:.1f}%)')
    if len(non_empty_age) > 0:
        print(f'Sample age values: {list(non_empty_age.unique())[:10]}')
    else:
        print('All age values are empty or NaN')

print('\nMAGE Ethnicity vs Age Comparison:')
print('=' * 40)

# Check ethnicity coverage
eth_values = adata.obs['self_reported_ethnicity'].dropna()
non_empty_eth = eth_values[eth_values != '']
print(f'Ethnicity coverage: {len(non_empty_eth)} ({len(non_empty_eth)/len(adata.obs)*100:.1f}%)')
print(f'Unique ethnicities: {list(non_empty_eth.unique())}')

print(f'\nTotal MAGE samples: {len(adata.obs)}')

# Check source metadata files for age availability
print('\n1000 Genomes Data Analysis:')
print('=' * 35)
print('MAGE represents 1000 Genomes Project samples (lymphoblastoid cell lines)')
print('1000 Genomes typically provides:')
print('  ✅ Population/ethnicity: Available (from population codes)')  
print('  ✅ Sex: Available (from PED files)')
print('  ❌ Age: NOT available (privacy protection - samples are immortalized cell lines)')
print('  ✅ Technical data: Available (RNA quality, sequencing metrics)')
print('\nReason: 1000 Genomes focuses on genetic variation, not age-related studies.')
print('Cell lines are immortalized, making original donor age less relevant.')