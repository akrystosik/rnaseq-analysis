#!/usr/bin/env python3
"""
Examine ADNI longitudinal diagnosis data to understand the discrepancy for subject 941_S_4420.

This script loads the H5AD file and examines the longitudinal diagnosis data stored in adata.uns
to understand how worst_diagnosis_code is calculated and why it differs from the individual record.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import json

def examine_longitudinal_diagnosis():
    """Examine longitudinal diagnosis data for subject 941_S_4420."""
    
    # Load ADNI dataset
    input_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
    print(f"Loading ADNI data from: {input_file}")
    adata = sc.read_h5ad(input_file)
    
    target_subject = '941_S_4420'
    target_rid = '4420'  # Extract RID from subject_id
    
    print(f"\n" + "="*70)
    print(f"EXAMINING LONGITUDINAL DIAGNOSIS DATA FOR SUBJECT {target_subject}")
    print(f"="*70)
    
    # 1. Check if longitudinal diagnosis data exists in adata.uns
    if 'longitudinal_diagnoses_adni' in adata.uns:
        longitudinal_data = adata.uns['longitudinal_diagnoses_adni']
        print(f"✅ Found longitudinal diagnosis data for {len(longitudinal_data)} subjects")
        
        # 2. Look for our target subject
        if target_rid in longitudinal_data:
            subject_history_raw = longitudinal_data[target_rid]
            print(f"\n🎯 Found diagnosis history for subject RID {target_rid}:")
            print(f"   Data type: {type(subject_history_raw)}")
            
            # Parse JSON string if needed
            if isinstance(subject_history_raw, str):
                try:
                    subject_history = json.loads(subject_history_raw)
                    print(f"   Successfully parsed JSON string")
                except json.JSONDecodeError as e:
                    print(f"   Error parsing JSON: {e}")
                    subject_history = []
            else:
                subject_history = subject_history_raw
            
            print(f"   Number of visits: {len(subject_history)}")
            
            # Display diagnosis entries
            print(f"\n📅 Complete diagnosis history:")
            if isinstance(subject_history, list):
                for i, entry in enumerate(subject_history[:10], 1):  # Show first 10 visits
                    print(f"   Visit {i}:")
                    if isinstance(entry, dict):
                        print(f"     - Visit Code: {entry.get('visit_code', 'N/A')}")
                        print(f"     - Exam Date: {entry.get('exam_date', 'N/A')}")
                        print(f"     - Diagnosis Code: {entry.get('diagnosis_code', 'N/A')}")
                        print(f"     - Diagnosis Label: {entry.get('diagnosis_label', 'N/A')}")
                    else:
                        print(f"     - Entry: {entry}")
                    print()
                if len(subject_history) > 10:
                    print(f"   ... and {len(subject_history) - 10} more visits")
            else:
                print(f"   Unexpected data structure: {type(subject_history)}")
                print(f"   Content: {subject_history}")
            
            # 3. Apply the worst diagnosis logic manually to verify
            print(f"🔍 MANUALLY CALCULATING WORST DIAGNOSIS:")
            diagnosis_severity_order = {1: 1, 2: 2, 3: 3, -4: 0}
            
            if isinstance(subject_history, list):
                valid_diagnoses = [entry for entry in subject_history if isinstance(entry, dict) and entry.get('diagnosis_code') is not None]
                print(f"   Valid diagnosis entries: {len(valid_diagnoses)}")
                
                if valid_diagnoses:
                    for entry in valid_diagnoses:
                        code = entry['diagnosis_code']
                        severity = diagnosis_severity_order.get(code, 0)
                        print(f"   - Code {code} ('{entry['diagnosis_label']}') has severity {severity}")
                    
                    # Find the worst diagnosis
                    worst_entry = max(valid_diagnoses, 
                                    key=lambda x: diagnosis_severity_order.get(x['diagnosis_code'], 0))
                    
                    print(f"\n🎯 CALCULATED WORST DIAGNOSIS:")
                    print(f"   - Code: {worst_entry['diagnosis_code']}")
                    print(f"   - Label: {worst_entry['diagnosis_label']}")
                    print(f"   - Visit: {worst_entry['visit_code']}")
                    print(f"   - Date: {worst_entry['exam_date']}")
                    
                    # 4. Compare with what's in the obs data
                    subject_obs = adata.obs[adata.obs['subject_id'] == target_subject].iloc[0]
                    stored_worst_code = subject_obs['worst_diagnosis_code']
                    stored_worst_label = subject_obs['worst_diagnosis_label']
                    
                    print(f"\n📊 COMPARISON WITH STORED VALUES:")
                    print(f"   - Calculated worst code: {worst_entry['diagnosis_code']}")
                    print(f"   - Stored worst code: {stored_worst_code}")
                    print(f"   - Match: {'✅ YES' if worst_entry['diagnosis_code'] == stored_worst_code else '❌ NO'}")
                    print(f"   - Calculated label: {worst_entry['diagnosis_label']}")
                    print(f"   - Stored label: {stored_worst_label}")
                    
                    if worst_entry['diagnosis_code'] != stored_worst_code:
                        print(f"\n⚠️  DISCREPANCY DETECTED!")
                        print(f"   This suggests an issue in the pipeline processing.")
                    else:
                        print(f"\n✅ Values match - the pipeline correctly calculated worst diagnosis as {stored_worst_code}")
                else:
                    print(f"   No valid diagnoses found")
            else:
                print(f"   Cannot process - subject_history is not a list")
        else:
            print(f"❌ Subject RID {target_rid} not found in longitudinal diagnosis data")
            print(f"   Available RIDs: {list(longitudinal_data.keys())[:10]}...")  # Show first 10
    else:
        print(f"❌ No longitudinal diagnosis data found in adata.uns")
    
    # 5. Show what the partner expects vs what we have
    print(f"\n" + "="*70)
    print(f"PARTNER EXPECTATION vs OUR DATA")
    print(f"="*70)
    
    subject_obs = adata.obs[adata.obs['subject_id'] == target_subject].iloc[0]
    
    print(f"Partner expects: Code 3 (Alzheimer's Disease)")
    print(f"Our export shows: Code {subject_obs['worst_diagnosis_code']} ({subject_obs['worst_diagnosis_label']})")
    print(f"Individual record shows: Code {subject_obs['diagnosis_code']} ({subject_obs['diagnosis']})")
    
    # 6. Check a few other subjects to see the pattern
    print(f"\n" + "="*70)
    print(f"CHECKING PATTERN IN OTHER SUBJECTS")
    print(f"="*70)
    
    if 'longitudinal_diagnoses_adni' in adata.uns:
        longitudinal_data = adata.uns['longitudinal_diagnoses_adni']
        
        # Sample a few subjects to check
        sample_rids = list(longitudinal_data.keys())[:5]
        
        for rid in sample_rids:
            subject_id = f"941_S_{rid}"  # Reconstruct subject_id
            subject_obs_data = adata.obs[adata.obs['subject_id'] == subject_id]
            
            if len(subject_obs_data) > 0:
                obs_row = subject_obs_data.iloc[0]
                history = longitudinal_data[rid]
                
                # Get all diagnosis codes from history
                codes_in_history = [entry['diagnosis_code'] for entry in history if entry.get('diagnosis_code') is not None]
                max_code_in_history = max(codes_in_history) if codes_in_history else None
                
                print(f"Subject {subject_id}:")
                print(f"  - Codes in history: {codes_in_history}")
                print(f"  - Max code in history: {max_code_in_history}")
                print(f"  - Stored worst_diagnosis_code: {obs_row['worst_diagnosis_code']}")
                print(f"  - Individual diagnosis_code: {obs_row['diagnosis_code']}")
                print()

def main():
    examine_longitudinal_diagnosis()

if __name__ == '__main__':
    main()