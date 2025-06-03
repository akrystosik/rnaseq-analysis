#!/usr/bin/env python3
"""
Export ADNI diagnosis data for partner sharing.

This script extracts donor ID, most severe diagnosis, and diagnosis date
from the processed ADNI dataset.
"""

import scanpy as sc
import pandas as pd
import argparse
from pathlib import Path

def export_adni_diagnosis(input_file=None, output_dir=None):
    """Export ADNI diagnosis data to CSV from H5AD file.
    
    Args:
        input_file: Path to input H5AD file
        output_dir: Output directory for CSV export
    """
    
    # Use H5AD as the single source of truth
    if input_file is None:
        input_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
    
    # Load ADNI dataset
    print(f"Loading ADNI data from: {input_file}")
    adata = sc.read_h5ad(input_file)
    
    # Extract diagnosis data for each unique subject
    diagnosis_columns = [
        'subject_id',
        'worst_diagnosis_code', 
        'worst_diagnosis_label',
        'worst_diagnosis_date',
        'worst_diagnosis_visit'
    ]
    
    # Get unique subjects (deduplicate by subject_id)
    subject_diagnosis = adata.obs[diagnosis_columns].drop_duplicates(subset=['subject_id']).copy()
    
    # Rename columns for partner clarity
    subject_diagnosis = subject_diagnosis.rename(columns={
        'subject_id': 'donor_id',
        'worst_diagnosis_code': 'most_severe_diagnosis_code',
        'worst_diagnosis_label': 'most_severe_diagnosis',
        'worst_diagnosis_date': 'diagnosis_date',
        'worst_diagnosis_visit': 'diagnosis_visit'
    })
    
    # Sort by donor_id for consistent output
    subject_diagnosis = subject_diagnosis.sort_values('donor_id')
    
    # Set output file
    if output_dir:
        output_file = Path(output_dir) / 'adni_diagnosis_export.csv'
    else:
        output_file = 'adni_diagnosis_export.csv'
    
    # Save to CSV
    subject_diagnosis.to_csv(output_file, index=False)
    
    print(f"✅ ADNI diagnosis data exported to: {output_file}")
    print(f"📊 Total subjects: {len(subject_diagnosis):,}")
    
    # Show summary statistics
    print("\n📋 Diagnosis distribution:")
    diagnosis_counts = subject_diagnosis['most_severe_diagnosis'].value_counts()
    for diagnosis, count in diagnosis_counts.items():
        print(f"  {diagnosis}: {count:,}")
    
    # Show sample data
    print("\n📋 Sample data:")
    print(subject_diagnosis.head(10).to_string(index=False))
    
    # Check for missing data
    missing_diagnosis = subject_diagnosis['most_severe_diagnosis'].isna().sum()
    missing_dates = subject_diagnosis['diagnosis_date'].isna().sum()
    
    print(f"\n🔍 Data quality:")
    print(f"  Missing diagnosis: {missing_diagnosis}")
    print(f"  Missing diagnosis dates: {missing_dates}")
    
    return output_file

def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(description='Export ADNI diagnosis data for partner sharing')
    parser.add_argument('--input-file', help='Input ADNI H5AD file')
    parser.add_argument('--output-dir', help='Output directory for CSV export')
    
    args = parser.parse_args()
    
    export_adni_diagnosis(input_file=args.input_file, output_dir=args.output_dir)

if __name__ == '__main__':
    main()