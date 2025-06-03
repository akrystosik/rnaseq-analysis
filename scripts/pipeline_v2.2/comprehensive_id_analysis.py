#!/usr/bin/env python3
"""
Comprehensive analysis of sample and subject ID patterns across datasets.
"""

import pandas as pd
import numpy as np
import json
import re
from collections import defaultdict, Counter

def analyze_id_structure(df, dataset_name):
    """
    Analyze the structure and patterns of IDs for a specific dataset.
    """
    dataset_df = df[df['dataset'] == dataset_name]
    
    if len(dataset_df) == 0:
        return {}
    
    analysis = {
        'total_samples': len(dataset_df),
        'unique_subjects': dataset_df['subject_id'].nunique(),
        'samples_per_subject': len(dataset_df) / dataset_df['subject_id'].nunique(),
        'data_types': dataset_df['data_type'].value_counts().to_dict()
    }
    
    # Analyze sample ID patterns
    sample_patterns = analyze_patterns(dataset_df['sample_id'].tolist(), f"{dataset_name}_samples")
    subject_patterns = analyze_patterns(dataset_df['subject_id'].tolist(), f"{dataset_name}_subjects")
    
    analysis['sample_id_patterns'] = sample_patterns
    analysis['subject_id_patterns'] = subject_patterns
    
    # Get examples
    analysis['sample_id_examples'] = dataset_df['sample_id'].head(10).tolist()
    analysis['subject_id_examples'] = dataset_df['subject_id'].unique()[:10].tolist()
    
    return analysis

def analyze_patterns(id_list, context):
    """
    Analyze patterns in a list of IDs.
    """
    patterns = {
        'total_count': len(id_list),
        'unique_count': len(set(id_list)),
        'length_distribution': {},
        'prefix_patterns': defaultdict(int),
        'character_patterns': defaultdict(int),
        'format_patterns': defaultdict(int)
    }
    
    # Length distribution
    lengths = [len(str(id_val)) for id_val in id_list]
    patterns['length_distribution'] = dict(Counter(lengths))
    patterns['avg_length'] = np.mean(lengths)
    
    # Pattern analysis
    for id_val in id_list:
        id_str = str(id_val)
        
        # Prefix patterns (first 3-5 characters)
        if len(id_str) >= 3:
            patterns['prefix_patterns'][id_str[:3]] += 1
        if len(id_str) >= 5:
            patterns['prefix_patterns'][id_str[:5]] += 1
        
        # Character composition
        if re.match(r'^[A-Z]+\d+$', id_str):
            patterns['character_patterns']['alpha_prefix_numeric_suffix'] += 1
        elif re.match(r'^\d+$', id_str):
            patterns['character_patterns']['numeric_only'] += 1
        elif re.match(r'^[A-Za-z]+$', id_str):
            patterns['character_patterns']['alpha_only'] += 1
        elif re.match(r'^[A-Z]+[A-Z0-9_-]+$', id_str):
            patterns['character_patterns']['complex_alphanumeric'] += 1
        else:
            patterns['character_patterns']['other'] += 1
        
        # Format-specific patterns
        if id_str.startswith('GTEX-'):
            if '-' in id_str:
                parts = id_str.split('-')
                if len(parts) >= 2:
                    patterns['format_patterns']['GTEX_subject_format'] += 1
                if len(parts) >= 4:
                    patterns['format_patterns']['GTEX_sample_format'] += 1
        elif id_str.startswith('ENCFF'):
            patterns['format_patterns']['ENCODE_file_format'] += 1
        elif id_str.startswith('ENCDO'):
            patterns['format_patterns']['ENCODE_donor_format'] += 1
        elif id_str.startswith('HG'):
            patterns['format_patterns']['1000G_HG_format'] += 1
        elif id_str.startswith('NA'):
            patterns['format_patterns']['1000G_NA_format'] += 1
        elif re.match(r'^\d{3}_S_\d{4}$', id_str):
            patterns['format_patterns']['ADNI_format'] += 1
    
    # Convert defaultdicts to regular dicts for JSON serialization
    patterns['prefix_patterns'] = dict(patterns['prefix_patterns'])
    patterns['character_patterns'] = dict(patterns['character_patterns'])
    patterns['format_patterns'] = dict(patterns['format_patterns'])
    
    return patterns

def main():
    """
    Main function to perform comprehensive ID analysis.
    """
    # Load the inventory data
    df = pd.read_csv('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/sample_subject_inventory.csv')
    
    print("COMPREHENSIVE RNA-SEQ PIPELINE SAMPLE/SUBJECT ID ANALYSIS")
    print("=" * 65)
    
    # Overall statistics
    total_samples = len(df)
    total_subjects = df['subject_id'].nunique()
    datasets = df['dataset'].unique()
    
    print(f"\nOVERALL STATISTICS:")
    print(f"  Total samples: {total_samples:,}")
    print(f"  Total unique subjects: {total_subjects:,}")
    print(f"  Number of datasets: {len(datasets)}")
    print(f"  Datasets: {', '.join(datasets)}")
    
    # Dataset-specific analysis
    comprehensive_analysis = {}
    
    for dataset in sorted(datasets):
        print(f"\n{dataset.upper()} DATASET ANALYSIS:")
        print("-" * 40)
        
        analysis = analyze_id_structure(df, dataset)
        comprehensive_analysis[dataset] = analysis
        
        print(f"  Total samples: {analysis['total_samples']:,}")
        print(f"  Unique subjects: {analysis['unique_subjects']:,}")
        print(f"  Samples per subject: {analysis['samples_per_subject']:.1f}")
        print(f"  Data types: {analysis['data_types']}")
        
        print(f"\n  Sample ID Structure:")
        sample_patterns = analysis['sample_id_patterns']
        print(f"    Average length: {sample_patterns['avg_length']:.1f} characters")
        print(f"    Length distribution: {sample_patterns['length_distribution']}")
        print(f"    Format patterns: {sample_patterns['format_patterns']}")
        
        print(f"\n  Subject ID Structure:")
        subject_patterns = analysis['subject_id_patterns']
        print(f"    Average length: {subject_patterns['avg_length']:.1f} characters")
        print(f"    Length distribution: {subject_patterns['length_distribution']}")
        print(f"    Format patterns: {subject_patterns['format_patterns']}")
        
        print(f"\n  Sample ID Examples: {analysis['sample_id_examples'][:5]}")
        print(f"  Subject ID Examples: {analysis['subject_id_examples'][:5]}")
    
    # Identifier pattern summary
    print(f"\n\nIDENTIFIER PATTERN SUMMARY:")
    print("=" * 40)
    
    print(f"\nENCODE:")
    print(f"  - Sample IDs: ENCFF###### format (ENCODE File Format)")
    print(f"  - Subject IDs: ENCDO###### format (ENCODE Donor Format)")
    print(f"  - Pattern: Cell line/donor-based identifiers")
    
    print(f"\nGTEx:")
    print(f"  - Sample IDs: GTEX-XXXXX-####-##-#### format")
    print(f"  - Subject IDs: GTEX-XXXXX format (subject prefix)")
    print(f"  - Pattern: Hierarchical tissue-sample-subject structure")
    
    print(f"\nMAGE (1000 Genomes):")
    print(f"  - Sample IDs: HG##### or NA##### format")
    print(f"  - Subject IDs: Same as sample IDs (individual-level data)")
    print(f"  - Pattern: Population genetics identifiers")
    
    print(f"\nADNI:")
    print(f"  - Sample IDs: ###_S_#### format")
    print(f"  - Subject IDs: Same as sample IDs (subject-level data)")
    print(f"  - Pattern: Clinical study participant identifiers")
    
    # Data type distribution
    print(f"\n\nDATA TYPE DISTRIBUTION:")
    print("=" * 30)
    data_type_dist = df['data_type'].value_counts()
    for data_type, count in data_type_dist.items():
        print(f"  {data_type}: {count:,} samples")
    
    # Save comprehensive analysis
    output_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/comprehensive_id_analysis.json'
    with open(output_file, 'w') as f:
        json.dump(comprehensive_analysis, f, indent=2, default=str)
    
    print(f"\n\nDetailed analysis saved to: {output_file}")
    print(f"Sample inventory CSV: /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/sample_subject_inventory.csv")

if __name__ == "__main__":
    main()