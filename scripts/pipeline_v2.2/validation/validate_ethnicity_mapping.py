#!/usr/bin/env python3
"""
Validate the ethnicity mapping files.
"""

import pandas as pd

def validate_mappings():
    # Load both mappings
    sample_df = pd.read_csv('sample_ethnicity_mapping_with_ontology.csv')
    subject_df = pd.read_csv('subject_ethnicity_mapping_with_ontology.csv')

    print('=== VALIDATION SUMMARY ===')

    # Check GTEx subject ID extraction
    gtex_samples = sample_df[sample_df['dataset'] == 'GTEx']['sample_id'].unique()
    gtex_subjects = subject_df[subject_df['dataset'] == 'GTEx']['subject_id'].unique()

    print(f'GTEx samples: {len(gtex_samples):,}')
    print(f'GTEx subjects: {len(gtex_subjects):,}')

    # Show subject ID format examples
    print('\nGTEx subject ID examples:')
    for subject in gtex_subjects[:5]:
        print(f'  {subject}')

    # Final dataset summary
    print('\n=== FINAL DATASET SUMMARY ===')
    print('Subject-level mapping (recommended for partner):')
    for dataset in ['ADNI', 'ENCODE', 'GTEx', 'MAGE']:
        count = len(subject_df[subject_df['dataset'] == dataset])
        print(f'  {dataset}: {count:,} subjects')

    print(f'\nTotal unique subjects: {len(subject_df):,}')

    # Coverage calculation
    known_subjects = len(subject_df[subject_df['hancestro_term'] != 'unknown'])
    total_subjects = len(subject_df)
    coverage_pct = known_subjects / total_subjects * 100

    print(f'Coverage with HANCESTRO terms: {known_subjects:,}/{total_subjects:,} ({coverage_pct:.1f}%)')

if __name__ == '__main__':
    validate_mappings()