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

def export_adni_diagnosis(input_file=None, output_dir=None, diagnosis_source='h5ad'):
    """Export ADNI diagnosis data to CSV.
    
    Args:
        input_file: Path to input file (H5AD or diagnosis CSV)
        output_dir: Output directory for CSV export
        diagnosis_source: 'h5ad' (default) or 'csv' (direct from source)
    """
    
    if diagnosis_source == 'csv':
        # Load diagnosis data directly from source CSV to include ALL subjects
        diagnosis_csv = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/AutoSync/subject_demographics/DXSUM_30Apr2025.csv'
        print(f"Loading ALL ADNI diagnosis data from source: {diagnosis_csv}")
        
        # Load and process all diagnosis data
        dx_df = pd.read_csv(diagnosis_csv, low_memory=False)
        dx_df['EXAMDATE'] = pd.to_datetime(dx_df['EXAMDATE'], errors='coerce')
        dx_df['RID'] = pd.to_numeric(dx_df['RID'], errors='coerce').astype('Int64')
        dx_df['DIAGNOSIS'] = pd.to_numeric(dx_df['DIAGNOSIS'], errors='coerce')
        
        # Diagnosis mapping
        diagnosis_mapping = {
            1: "Cognitively Normal", 
            2: "Mild Cognitive Impairment", 
            3: "Alzheimer's Disease",
            -4: "Missing/Unknown"
        }
        
        # Calculate worst diagnosis for each subject
        diagnosis_severity_order = {1: 1, 2: 2, 3: 3, -4: 0}
        
        subject_diagnosis_list = []
        # Don't require EXAMDATE - use USERDATE as backup
        dx_df_clean = dx_df.dropna(subset=['RID', 'DIAGNOSIS'])
        
        for rid, group in dx_df_clean.groupby('RID'):
            # Use the actual PTID from the data instead of reconstructing
            subject_id = group['PTID'].iloc[0]  # Get PTID from first row of this RID group
            
            # Get all valid diagnoses for this subject
            valid_diagnoses = []
            # Sort by EXAMDATE, with NaT values at the end, then by USERDATE
            group_sorted = group.sort_values(by=['EXAMDATE', 'USERDATE'], na_position='last')
            
            for _, row in group_sorted.iterrows():
                diag_code = row['DIAGNOSIS']
                diag_label = diagnosis_mapping.get(diag_code, f"Unknown code: {diag_code}")
                
                # Use EXAMDATE if available, otherwise USERDATE
                exam_date = row['EXAMDATE']
                if pd.isna(exam_date):
                    exam_date = pd.to_datetime(row['USERDATE'], errors='coerce')
                
                # Use VISCODE if VISCODE2 is empty
                visit_code = str(row.get('VISCODE2', ''))
                if not visit_code or visit_code == 'nan' or visit_code == '':
                    visit_code = str(row.get('VISCODE', 'N/A'))
                
                valid_diagnoses.append({
                    "visit_code": visit_code,
                    "exam_date": str(exam_date.date()) if pd.notna(exam_date) else 'N/A',
                    "diagnosis_code": int(diag_code) if pd.notna(diag_code) else None,
                    "diagnosis_label": diag_label
                })
            
            if valid_diagnoses:
                # Find worst diagnosis
                worst_entry = max(valid_diagnoses, 
                                key=lambda x: diagnosis_severity_order.get(x['diagnosis_code'], 0))
                
                subject_diagnosis_list.append({
                    'donor_id': subject_id,
                    'most_severe_diagnosis_code': worst_entry['diagnosis_code'],
                    'most_severe_diagnosis': worst_entry['diagnosis_label'],
                    'diagnosis_date': worst_entry['exam_date'],
                    'diagnosis_visit': worst_entry['visit_code']
                })
        
        subject_diagnosis = pd.DataFrame(subject_diagnosis_list)
        print(f"✅ Processed ALL subjects from diagnosis source: {len(subject_diagnosis):,}")
        
    else:
        # Original H5AD-based method (expression data subjects only)
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
    parser.add_argument('--source', choices=['h5ad', 'csv'], default='h5ad',
                       help='Data source: h5ad (expression subjects only) or csv (all subjects)')
    
    args = parser.parse_args()
    
    export_adni_diagnosis(input_file=args.input_file, output_dir=args.output_dir, 
                         diagnosis_source=args.source)

if __name__ == '__main__':
    main()