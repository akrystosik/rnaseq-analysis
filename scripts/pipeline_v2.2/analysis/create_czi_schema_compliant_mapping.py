#!/usr/bin/env python3
"""
Create CZI single-cell curation schema compliant ethnicity mapping.
Field name: self_reported_ethnicity_ontology_term_id
"""

import pandas as pd

def create_czi_compliant_mapping():
    """Create ethnicity mapping compliant with CZI schema v3.0.0."""
    
    # Load our current mapping
    df = pd.read_csv('subject_ethnicity_mapping_with_ontology.csv')
    
    # Rename columns to match CZI schema
    czi_df = df.copy()
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

if __name__ == '__main__':
    create_czi_compliant_mapping()