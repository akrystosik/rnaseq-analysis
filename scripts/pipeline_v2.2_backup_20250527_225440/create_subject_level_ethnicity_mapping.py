#!/usr/bin/env python3
"""
Create subject-level ethnicity mapping (deduplicated) for partner ancestry report.
"""

import pandas as pd

def create_subject_level_mapping():
    """Create deduplicated subject-level ethnicity mapping."""
    
    # Load the sample-level mapping
    df = pd.read_csv('sample_ethnicity_mapping_with_ontology.csv')
    
    print(f"📊 Original sample-level data: {len(df):,} entries")
    
    # For GTEx, extract subject ID and deduplicate
    # For other datasets, use sample_id as subject_id
    def get_subject_id(row):
        if row['dataset'] == 'GTEx':
            # GTEx sample format: GTEX-XXXXX-YYYY-ZZ -> subject is GTEX-XXXXX
            sample_id = row['sample_id']
            if isinstance(sample_id, str) and sample_id.startswith('GTEX-'):
                parts = sample_id.split('-')
                if len(parts) >= 2:
                    return f"{parts[0]}-{parts[1]}"  # GTEX-XXXXX
                else:
                    return sample_id  # fallback
            else:
                return sample_id
        else:
            # For ADNI, ENCODE, MAGE - sample_id is subject_id
            return row['sample_id']
    
    # Add subject_id column
    df['subject_id'] = df.apply(get_subject_id, axis=1)
    
    # Deduplicate by subject_id and dataset
    # Keep first occurrence (they should all have same ethnicity anyway)
    subject_df = df.drop_duplicates(subset=['subject_id', 'dataset'], keep='first')
    
    # Reorder columns for clarity
    subject_df = subject_df[['subject_id', 'dataset', 'ethnicity', 'hancestro_term']]
    subject_df = subject_df.sort_values(['dataset', 'subject_id'])
    
    # Save subject-level mapping
    output_file = 'subject_ethnicity_mapping_with_ontology.csv'
    subject_df.to_csv(output_file, index=False)
    
    print(f"✅ Subject-level ethnicity mapping saved to: {output_file}")
    print(f"📊 Deduplicated to: {len(subject_df):,} unique subjects")
    
    # Show before/after comparison
    print("\n📋 Before/After comparison:")
    print("Sample-level (with tissue replicates):")
    sample_counts = df['dataset'].value_counts().sort_index()
    for dataset, count in sample_counts.items():
        print(f"  {dataset}: {count:,}")
    
    print("\nSubject-level (deduplicated):")
    subject_counts = subject_df['dataset'].value_counts().sort_index()
    for dataset, count in subject_counts.items():
        print(f"  {dataset}: {count:,}")
    
    # Show ethnicity distribution for subjects
    print("\n🧬 Subject-level HANCESTRO ontology distribution:")
    hancestro_counts = subject_df['hancestro_term'].value_counts()
    for term, count in hancestro_counts.head(10).items():
        print(f"  {term}: {count:,}")
    
    # Validation check
    print(f"\n✅ Validation:")
    print(f"  Unique subject_ids: {subject_df['subject_id'].nunique():,}")
    print(f"  Total rows: {len(subject_df):,}")
    print(f"  No duplicates: {subject_df['subject_id'].nunique() == len(subject_df)}")
    
    return output_file

if __name__ == '__main__':
    create_subject_level_mapping()