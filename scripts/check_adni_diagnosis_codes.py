#!/usr/bin/env python3
"""
Check ADNI diagnosis codes to understand the distribution and look for subjects who progressed to AD (Code 3).
"""

import scanpy as sc
import pandas as pd
import numpy as np
import json

def check_diagnosis_distribution():
    """Check ADNI diagnosis code distribution and find examples of different codes."""
    
    # Load ADNI dataset
    input_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
    print(f"Loading ADNI data from: {input_file}")
    adata = sc.read_h5ad(input_file)
    
    print(f"\n" + "="*70)
    print(f"DIAGNOSIS CODE DISTRIBUTION ANALYSIS")
    print(f"="*70)
    
    # 1. Check distribution of worst_diagnosis_code
    worst_diagnosis_counts = adata.obs['worst_diagnosis_code'].value_counts().sort_index()
    print(f"Distribution of worst_diagnosis_code:")
    for code, count in worst_diagnosis_counts.items():
        label = adata.obs[adata.obs['worst_diagnosis_code'] == code]['worst_diagnosis_label'].iloc[0]
        print(f"  Code {code} ({label}): {count:,} subjects")
    
    # 2. Look for subjects with Code 3 (Alzheimer's Disease)
    ad_subjects = adata.obs[adata.obs['worst_diagnosis_code'] == 3]
    print(f"\n🔍 Subjects with Code 3 (Alzheimer's Disease): {len(ad_subjects)}")
    
    if len(ad_subjects) > 0:
        print(f"Sample AD subjects:")
        for i, (idx, row) in enumerate(ad_subjects.head(5).iterrows()):
            print(f"  {i+1}. {row['subject_id']} - {row['worst_diagnosis_label']} (visit: {row['worst_diagnosis_visit']}, date: {row['worst_diagnosis_date']})")
        
        # Check one AD subject's longitudinal data
        if 'longitudinal_diagnoses_adni' in adata.uns:
            sample_ad_subject = ad_subjects.iloc[0]['subject_id']
            sample_rid = sample_ad_subject.split('_S_')[-1]
            
            longitudinal_data = adata.uns['longitudinal_diagnoses_adni']
            if sample_rid in longitudinal_data:
                print(f"\n📅 Longitudinal history for AD subject {sample_ad_subject}:")
                subject_history = json.loads(longitudinal_data[sample_rid])
                for entry in subject_history:
                    print(f"  - {entry['visit_code']}: Code {entry['diagnosis_code']} ({entry['diagnosis_label']}) on {entry['exam_date']}")
    
    # 3. Check if subject 941_S_4420 could have progressed to AD in later visits not captured
    print(f"\n🎯 INVESTIGATING SUBJECT 941_S_4420 PROGRESSION:")
    target_subject = '941_S_4420'
    target_rid = '4420'
    
    if 'longitudinal_diagnoses_adni' in adata.uns:
        longitudinal_data = adata.uns['longitudinal_diagnoses_adni']
        if target_rid in longitudinal_data:
            subject_history = json.loads(longitudinal_data[target_rid])
            print(f"Complete diagnosis timeline for {target_subject}:")
            
            # Sort by date to see progression
            sorted_history = sorted(subject_history, key=lambda x: x['exam_date'])
            for entry in sorted_history:
                print(f"  {entry['exam_date']}: {entry['visit_code']} - Code {entry['diagnosis_code']} ({entry['diagnosis_label']})")
            
            # Check if there are any later visits that might not be in our data
            latest_date = max(entry['exam_date'] for entry in subject_history)
            print(f"\n  Latest recorded visit: {latest_date}")
            print(f"  Highest diagnosis code seen: {max(entry['diagnosis_code'] for entry in subject_history)}")
    
    # 4. Check if there's a pattern of missing progression data
    print(f"\n" + "="*70)
    print(f"POTENTIAL MISSING PROGRESSION DATA ANALYSIS")
    print(f"="*70)
    
    # Look at visit codes to understand data completeness
    if 'longitudinal_diagnoses_adni' in adata.uns:
        longitudinal_data = adata.uns['longitudinal_diagnoses_adni']
        
        all_visit_codes = set()
        progression_subjects = []
        
        for rid, history_str in longitudinal_data.items():
            history = json.loads(history_str)
            visit_codes = [entry['visit_code'] for entry in history]
            all_visit_codes.update(visit_codes)
            
            # Check for progression (increasing diagnosis codes)
            codes = [entry['diagnosis_code'] for entry in sorted(history, key=lambda x: x['exam_date'])]
            if len(set(codes)) > 1:  # Has multiple different codes
                progression_subjects.append((rid, codes))
        
        print(f"All visit codes found: {sorted(all_visit_codes)}")
        print(f"Subjects with diagnosis progression: {len(progression_subjects)}")
        
        if progression_subjects:
            print(f"\nSample progression patterns:")
            for i, (rid, codes) in enumerate(progression_subjects[:5]):
                print(f"  Subject {rid}: {codes}")

def main():
    check_diagnosis_distribution()

if __name__ == '__main__':
    main()