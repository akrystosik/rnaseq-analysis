#!/usr/bin/env python3

import scanpy as sc
import pandas as pd

# Load all datasets and check developmental stage ontology coverage
data_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137'
dataset_files = {
    'ADNI': 'adni_standardized_preprocessed.h5ad',
    'ENCODE': 'encode_standardized_preprocessed.h5ad', 
    'GTEx': 'gtex_standardized_preprocessed.h5ad',
    'MAGE': 'mage_standardized_preprocessed.h5ad'
}

print('Developmental Stage Ontology Analysis:')
print('=' * 45)

for name, filename in dataset_files.items():
    file_path = f'{data_dir}/{filename}'
    try:
        adata = sc.read_h5ad(file_path)
        
        print(f'\n{name} Dataset:')
        print(f'  Total samples: {adata.n_obs:,}')
        
        # Check if developmental_stage_ontology exists
        if 'developmental_stage_ontology' in adata.obs.columns:
            dev_stage_values = adata.obs['developmental_stage_ontology']
            
            # Count non-empty values
            non_empty = (dev_stage_values != '') & (dev_stage_values.notna())
            coverage = non_empty.sum() / adata.n_obs * 100
            
            print(f'  Developmental stage ontology coverage: {coverage:.1f}%')
            
            # Show unique values
            unique_values = dev_stage_values.unique()
            print(f'  Unique values: {list(unique_values)}')
            
            # Check what's missing
            empty_count = (dev_stage_values == '').sum()
            na_count = dev_stage_values.isna().sum()
            print(f'  Empty strings: {empty_count}')
            print(f'  NaN values: {na_count}')
            
        else:
            print('  ❌ No developmental_stage_ontology column')
        
        # Check for age data to map
        if 'age' in adata.obs.columns:
            age_values = adata.obs['age'].dropna()
            age_values = age_values[age_values != '']
            print(f'  Available age data: {len(age_values)} samples')
            if len(age_values) > 0:
                print(f'  Age value examples: {list(age_values.unique())[:5]}')
                
                # Try to categorize ages for HsapDv mapping
                try:
                    # Convert to numeric where possible
                    numeric_ages = []
                    for age in age_values:
                        try:
                            if isinstance(age, str) and '-' in age:
                                # Handle age ranges like "60-69"
                                start, end = age.split('-')
                                numeric_ages.append((float(start) + float(end)) / 2)
                            else:
                                numeric_ages.append(float(age))
                        except:
                            continue
                    
                    if numeric_ages:
                        print(f'  Numeric age range: {min(numeric_ages):.1f} - {max(numeric_ages):.1f}')
                        print(f'  Mean age: {sum(numeric_ages)/len(numeric_ages):.1f}')
                        
                        # Categorize for HsapDv mapping
                        categories = {}
                        for age in numeric_ages:
                            if age < 18:
                                categories.setdefault('child/adolescent', 0)
                                categories['child/adolescent'] += 1
                            elif age < 65:
                                categories.setdefault('adult', 0)
                                categories['adult'] += 1
                            else:
                                categories.setdefault('elderly', 0)
                                categories['elderly'] += 1
                        
                        print(f'  Age categories: {categories}')
                except Exception as e:
                    print(f'  Age parsing error: {e}')
                
    except Exception as e:
        print(f'  ❌ Error loading {name}: {e}')

print('\nHsapDv Ontology Terms for Age Mapping:')
print('=' * 45)
print('Common HsapDv terms for human development:')
print('  • HsapDv:0000087 - human adult stage (18+ years)')
print('  • HsapDv:0000224 - human aged stage (65+ years)')  
print('  • HsapDv:0000083 - human adolescent stage (13-17 years)')
print('  • HsapDv:0000082 - human child stage (3-12 years)')
print('  • HsapDv:0000003 - human adult stage, immature')
print('\nRecommendation: Map based on available age data to appropriate HsapDv terms')