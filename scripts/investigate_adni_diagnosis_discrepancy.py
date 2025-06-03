#!/usr/bin/env python3
"""
Investigate ADNI diagnosis discrepancy for subject 941_S_4420.

The partner reports this subject should have diagnosis code 3, but our exported CSV shows code 2.
This script will examine all records for this subject and understand how worst_diagnosis_code is calculated.
"""

import scanpy as sc
import pandas as pd
import numpy as np

def investigate_diagnosis_discrepancy():
    """Investigate diagnosis discrepancy for subject 941_S_4420."""
    
    # Load ADNI dataset
    input_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
    print(f"Loading ADNI data from: {input_file}")
    adata = sc.read_h5ad(input_file)
    
    print(f"Total observations: {adata.n_obs:,}")
    print(f"Total unique subjects: {adata.obs['subject_id'].nunique():,}")
    
    # 1. Examine all records for subject 941_S_4420
    target_subject = '941_S_4420'
    subject_records = adata.obs[adata.obs['subject_id'] == target_subject].copy()
    
    print(f"\n" + "="*60)
    print(f"INVESTIGATION FOR SUBJECT {target_subject}")
    print(f"="*60)
    
    if len(subject_records) == 0:
        print(f"❌ No records found for subject {target_subject}")
        return
    
    print(f"✅ Found {len(subject_records)} records for subject {target_subject}")
    
    # 2. Look at diagnosis-related columns
    diagnosis_cols = [col for col in adata.obs.columns if 'diagnosis' in col.lower()]
    print(f"\n📋 Available diagnosis columns:")
    for col in diagnosis_cols:
        print(f"  - {col}")
    
    # Show all diagnosis data for this subject
    print(f"\n📊 All diagnosis data for subject {target_subject}:")
    relevant_cols = ['subject_id', 'sample_id', 'visit_id'] + diagnosis_cols
    available_cols = [col for col in relevant_cols if col in subject_records.columns]
    
    print(subject_records[available_cols].to_string(index=False))
    
    # 3. Check for multiple visits/records with different diagnosis codes
    if 'diagnosis_code' in subject_records.columns:
        unique_codes = subject_records['diagnosis_code'].unique()
        print(f"\n🔍 Unique diagnosis codes for this subject: {unique_codes}")
        
        # Show breakdown by visit if available
        if 'visit_id' in subject_records.columns:
            visit_diagnosis = subject_records[['visit_id', 'diagnosis_code', 'diagnosis_label']].drop_duplicates()
            print(f"\n📅 Diagnosis by visit:")
            print(visit_diagnosis.to_string(index=False))
    
    # 4. Examine worst_diagnosis_code calculation
    if 'worst_diagnosis_code' in subject_records.columns:
        worst_code = subject_records['worst_diagnosis_code'].iloc[0]
        print(f"\n🎯 Current worst_diagnosis_code: {worst_code}")
        
        # Check if there are any diagnosis codes higher than the worst_diagnosis_code
        if 'diagnosis_code' in subject_records.columns:
            all_codes = subject_records['diagnosis_code'].dropna()
            if len(all_codes) > 0:
                max_code = all_codes.max()
                print(f"🎯 Maximum diagnosis_code in records: {max_code}")
                
                if max_code != worst_code:
                    print(f"⚠️  DISCREPANCY FOUND: worst_diagnosis_code ({worst_code}) != max diagnosis_code ({max_code})")
                else:
                    print(f"✅ worst_diagnosis_code matches maximum diagnosis_code")
    
    # 5. Check a few other subjects to see if this is widespread
    print(f"\n" + "="*60)
    print(f"CHECKING OTHER SUBJECTS FOR SIMILAR ISSUES")
    print(f"="*60)
    
    # Get a sample of other subjects
    unique_subjects = adata.obs['subject_id'].unique()
    sample_subjects = np.random.choice(unique_subjects, size=min(10, len(unique_subjects)), replace=False)
    
    discrepancy_count = 0
    checked_subjects = 0
    
    for subject in sample_subjects:
        if subject == target_subject:
            continue
            
        subject_data = adata.obs[adata.obs['subject_id'] == subject]
        
        if 'worst_diagnosis_code' in subject_data.columns and 'diagnosis_code' in subject_data.columns:
            worst_code = subject_data['worst_diagnosis_code'].iloc[0]
            all_codes = subject_data['diagnosis_code'].dropna()
            
            if len(all_codes) > 0:
                max_code = all_codes.max()
                checked_subjects += 1
                
                if max_code != worst_code:
                    discrepancy_count += 1
                    print(f"⚠️  Subject {subject}: worst_diagnosis_code ({worst_code}) != max diagnosis_code ({max_code})")
    
    print(f"\n📊 Summary of random sample check:")
    print(f"  - Checked {checked_subjects} subjects")
    print(f"  - Found {discrepancy_count} subjects with discrepancies")
    print(f"  - Discrepancy rate: {discrepancy_count/checked_subjects*100:.1f}%" if checked_subjects > 0 else "  - No subjects to check")
    
    # 6. Look at how worst_diagnosis is supposed to be calculated
    print(f"\n" + "="*60)
    print(f"UNDERSTANDING DIAGNOSIS CODE LOGIC")
    print(f"="*60)
    
    # Check diagnosis code distribution
    if 'diagnosis_code' in adata.obs.columns:
        print(f"\n📊 Overall diagnosis code distribution:")
        code_counts = adata.obs['diagnosis_code'].value_counts().sort_index()
        for code, count in code_counts.items():
            print(f"  Code {code}: {count:,} records")
    
    if 'diagnosis_label' in adata.obs.columns:
        print(f"\n📊 Diagnosis labels:")
        label_counts = adata.obs['diagnosis_label'].value_counts()
        for label, count in label_counts.items():
            print(f"  {label}: {count:,} records")
    
    # Look for any documentation or comments about diagnosis severity
    print(f"\n🔍 Checking for diagnosis severity mapping...")
    
    # Create a mapping of codes to labels to understand severity order
    if 'diagnosis_code' in adata.obs.columns and 'diagnosis_label' in adata.obs.columns:
        code_label_map = adata.obs[['diagnosis_code', 'diagnosis_label']].drop_duplicates().sort_values('diagnosis_code')
        print(f"\n📋 Code to Label mapping (sorted by code):")
        print(code_label_map.to_string(index=False))
    
    return subject_records

if __name__ == '__main__':
    investigate_diagnosis_discrepancy()