#!/usr/bin/env python3
"""
Enhance the complete ethnicity mapping with HANCESTRO ontology terms.
"""

import pandas as pd
import json

def load_hancestro_mapping():
    """Load the ethnicity to HANCESTRO ontology mapping."""
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/ethnicity_to_hancestro.json', 'r') as f:
        return json.load(f)

def create_ethnicity_mapping_with_ontology():
    """Add HANCESTRO ontology terms to the complete ethnicity mapping."""
    
    # Load the complete ethnicity mapping
    df = pd.read_csv('sample_ethnicity_mapping_complete.csv')
    
    # Load HANCESTRO mapping
    hancestro_map = load_hancestro_mapping()
    
    # Add ontology terms
    df['hancestro_term'] = df['ethnicity'].map(hancestro_map)
    
    # Handle any unmapped ethnicities
    unmapped = df[df['hancestro_term'].isna()]
    if len(unmapped) > 0:
        print(f"⚠️  Found {len(unmapped)} unmapped ethnicities:")
        for ethnicity in unmapped['ethnicity'].unique():
            print(f"  - {ethnicity}")
        # Set unmapped to 'unknown'
        df['hancestro_term'] = df['hancestro_term'].fillna('unknown')
    
    # Reorder columns for better readability
    df = df[['sample_id', 'dataset', 'ethnicity', 'hancestro_term']]
    
    # Save enhanced mapping
    output_file = 'sample_ethnicity_mapping_with_ontology.csv'
    df.to_csv(output_file, index=False)
    
    print(f"✅ Enhanced ethnicity mapping with ontology terms saved to: {output_file}")
    print(f"📊 Total samples: {len(df):,}")
    
    # Show ontology term distribution
    print("\n🧬 HANCESTRO ontology term distribution:")
    hancestro_counts = df['hancestro_term'].value_counts()
    for term, count in hancestro_counts.head(15).items():
        print(f"  {term}: {count:,}")
    
    # Show dataset summary with ontology coverage
    print("\n📋 Dataset summary with ontology coverage:")
    summary = df.groupby('dataset').agg({
        'sample_id': 'count',
        'hancestro_term': lambda x: (x != 'unknown').sum()
    }).rename(columns={'sample_id': 'total_samples', 'hancestro_term': 'with_hancestro'})
    summary['hancestro_coverage'] = (summary['with_hancestro'] / summary['total_samples'] * 100).round(1)
    print(summary.to_string())
    
    return output_file

if __name__ == '__main__':
    create_ethnicity_mapping_with_ontology()