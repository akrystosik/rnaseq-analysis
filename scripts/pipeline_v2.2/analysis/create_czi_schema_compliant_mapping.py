#!/usr/bin/env python3
"""
Create CZI single-cell curation schema compliant ethnicity mapping.
Field name: self_reported_ethnicity_ontology_term_id
"""

import pandas as pd
import argparse
from pathlib import Path

def create_czi_compliant_mapping(input_file=None, output_dir=None):
    """Create ethnicity mapping compliant with CZI schema v3.0.0."""
    
    # Use sample-level mapping as input (no intermediate files needed)
    if input_file is None:
        input_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/sample_ethnicity_mapping_with_ontology.csv'
    
    # Load sample-level mapping
    df = pd.read_csv(input_file)
    print(f"📊 Original sample-level data: {len(df):,} entries")
    
    # Create subject-level deduplication first
    def get_subject_id(row):
        if row['dataset'] == 'GTEx':
            # GTEx sample format: GTEX-XXXXX-YYYY-ZZ -> subject is GTEX-XXXXX
            sample_id = row['sample_id']
            if isinstance(sample_id, str) and sample_id.startswith('GTEX-'):
                parts = sample_id.split('-')
                if len(parts) >= 2:
                    return f"{parts[0]}-{parts[1]}"  # GTEX-XXXXX
            return sample_id
        else:
            # For other datasets, sample_id is the subject_id
            return row['sample_id']
    
    # Add subject_id column
    df['subject_id'] = df.apply(get_subject_id, axis=1)
    
    # Deduplicate to subject level (keep first occurrence for each subject)
    subject_df = df.drop_duplicates(subset=['subject_id'], keep='first')
    print(f"📊 Deduplicated to subject-level: {len(subject_df):,} unique subjects")
    
    # Now create CZI compliant version from subject-level data
    czi_df = subject_df.copy()
    czi_df = czi_df.rename(columns={
        'hancestro_term': 'self_reported_ethnicity_ontology_term_id'
    })
    
    # Keep ethnicity for human readability but make it clear it's for reference
    czi_df = czi_df.rename(columns={
        'ethnicity': 'ethnicity_label_for_reference'
    })
    
    # Reorder columns to match expected schema format
    czi_df = czi_df[['subject_id', 'dataset', 'self_reported_ethnicity_ontology_term_id', 'ethnicity_label_for_reference']]
    
    # Save CZI compliant version
    if output_dir:
        output_file = Path(output_dir) / 'subject_ethnicity_mapping_czi_compliant.csv'
    else:
        output_file = 'subject_ethnicity_mapping_czi_compliant.csv'
    czi_df.to_csv(output_file, index=False)
    
    print(f"✅ CZI schema compliant mapping saved to: {output_file}")
    print(f"📊 Total subjects: {len(czi_df):,}")
    
    # Validate against schema requirements
    print("\n=== CZI SCHEMA VALIDATION ===")
    
    ontology_terms = czi_df['self_reported_ethnicity_ontology_term_id'].value_counts()
    print("self_reported_ethnicity_ontology_term_id distribution:")
    for term, count in ontology_terms.items():
        print(f"  {term}: {count:,}")
    
    # Check compliance
    valid_terms = 0
    total_terms = len(ontology_terms)
    
    for term in ontology_terms.index:
        if (term.startswith('HANCESTRO:') and term.split(':')[1].isdigit()) or term in ['multiethnic', 'unknown']:
            valid_terms += 1
        else:
            print(f"❌ Invalid term: {term}")
    
    print(f"\n✅ Schema compliance: {valid_terms}/{total_terms} terms valid")
    
    # Show sample data
    print("\n📋 Sample CZI compliant data:")
    print(czi_df.head(10).to_string(index=False))
    
    return output_file

def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(description='Create CZI schema compliant ethnicity mapping')
    parser.add_argument('--preprocessed-dir', help='Preprocessed directory (not used)')
    parser.add_argument('--output-dir', help='Output directory for results')
    parser.add_argument('--force', action='store_true', help='Force regeneration')
    
    args = parser.parse_args()
    
    create_czi_compliant_mapping(output_dir=args.output_dir)

if __name__ == '__main__':
    main()