#!/usr/bin/env python3
"""Test the ADNI diagnosis fix by processing diagnosis data directly."""

import pandas as pd
import json

def test_adni_diagnosis_processing():
    """Test the updated ADNI diagnosis processing logic."""
    
    print("🧪 TESTING ADNI DIAGNOSIS PROCESSING FIX")
    print("=" * 50)
    
    # Load diagnosis CSV with same logic as pipeline
    diagnosis_csv = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/AutoSync/subject_demographics/DXSUM_30Apr2025.csv'
    print(f"📂 Loading diagnosis CSV: {diagnosis_csv}")
    
    dx_df_raw = pd.read_csv(diagnosis_csv, low_memory=False)
    dx_df_raw['EXAMDATE'] = pd.to_datetime(dx_df_raw['EXAMDATE'], errors='coerce')
    dx_df_raw['RID'] = pd.to_numeric(dx_df_raw['RID'], errors='coerce').astype('Int64')
    dx_df_raw['DIAGNOSIS'] = pd.to_numeric(dx_df_raw['DIAGNOSIS'], errors='coerce')
    
    print(f"✅ Loaded {len(dx_df_raw)} total diagnosis records")
    
    # Test subjects including 941_S_4420
    test_subjects = ['941_S_4420', '002_S_0413', '002_S_0729']
    test_rids = [4420, 413, 729]
    
    # Diagnosis mapping
    diag_map_codes = {
        1: "Cognitively Normal", 
        2: "Mild Cognitive Impairment", 
        3: "Alzheimer's Disease",
        -4: "Missing/Unknown"
    }
    
    # NEW LOGIC: Don't require EXAMDATE - use USERDATE as backup
    print(f"\n🔧 Testing NEW diagnosis processing logic...")
    dx_df_raw_cleaned = dx_df_raw.dropna(subset=['RID', 'DIAGNOSIS'])
    print(f"After dropping invalid records: {len(dx_df_raw_cleaned)} records")
    
    # Process each test subject
    all_subject_diagnoses = {}
    
    for rid in test_rids:
        subject_id = f"941_S_{rid}" if rid == 4420 else f"002_S_{rid:04d}"
        
        # Get records for this RID
        subject_records = dx_df_raw_cleaned[dx_df_raw_cleaned['RID'] == rid]
        
        if len(subject_records) == 0:
            print(f"\n❌ No records found for RID {rid} ({subject_id})")
            continue
            
        print(f"\n🎯 Processing RID {rid} ({subject_id}):")
        print(f"  Found {len(subject_records)} records")
        
        # Sort by EXAMDATE, with NaT values at the end, then by USERDATE
        group_sorted = subject_records.sort_values(by=['EXAMDATE', 'USERDATE'], na_position='last')
        
        subject_diagnoses = []
        for _, row in group_sorted.iterrows():
            diag_code = row['DIAGNOSIS']
            diag_label = diag_map_codes.get(diag_code, f"Unknown code: {diag_code}")
            
            # Use EXAMDATE if available, otherwise USERDATE
            exam_date_val = row['EXAMDATE']
            if pd.isna(exam_date_val):
                exam_date_val = pd.to_datetime(row['USERDATE'], errors='coerce')
            
            # Use VISCODE if VISCODE2 is empty
            visit_code = str(row.get('VISCODE2', ''))
            if not visit_code or visit_code == 'nan' or visit_code == '':
                visit_code = str(row.get('VISCODE', 'N/A'))
            
            diagnosis_entry = {
                "visit_code": visit_code,
                "exam_date": str(exam_date_val.date()) if pd.notna(exam_date_val) and hasattr(exam_date_val, 'date') else 'N/A',
                "diagnosis_code": int(diag_code) if pd.notna(diag_code) else None,
                "diagnosis_label": diag_label
            }
            
            subject_diagnoses.append(diagnosis_entry)
            print(f"    {diagnosis_entry['exam_date']}: Code {diagnosis_entry['diagnosis_code']} ({diagnosis_entry['diagnosis_label']}) at {diagnosis_entry['visit_code']}")
        
        if subject_diagnoses:
            all_subject_diagnoses[str(rid)] = subject_diagnoses
    
    # Test worst diagnosis calculation
    print(f"\n📊 TESTING WORST DIAGNOSIS CALCULATION:")
    diagnosis_severity_order = {1: 1, 2: 2, 3: 3, -4: 0}
    
    worst_diagnosis_map = {}
    for subject_rid_str, diagnosis_history in all_subject_diagnoses.items():
        subject_rid = int(subject_rid_str)
        valid_diagnoses = [entry for entry in diagnosis_history if entry.get('diagnosis_code') is not None]
        
        if valid_diagnoses:
            # Find diagnosis with highest severity
            worst_entry = max(valid_diagnoses, 
                            key=lambda x: diagnosis_severity_order.get(x['diagnosis_code'], 0))
            worst_diagnosis_map[subject_rid] = {
                'worst_diagnosis_code': worst_entry['diagnosis_code'],
                'worst_diagnosis_label': worst_entry['diagnosis_label'],
                'worst_diagnosis_visit': worst_entry['visit_code'],
                'worst_diagnosis_date': worst_entry['exam_date']
            }
            
            subject_id = f"941_S_{subject_rid}" if subject_rid == 4420 else f"002_S_{subject_rid:04d}"
            print(f"  {subject_id}: Worst = Code {worst_entry['diagnosis_code']} ({worst_entry['diagnosis_label']}) on {worst_entry['exam_date']}")
        else:
            print(f"  RID {subject_rid}: No valid diagnoses found")
    
    # Check specific case
    print(f"\n🎯 SPECIFIC TEST - 941_S_4420:")
    if 4420 in worst_diagnosis_map:
        result_4420 = worst_diagnosis_map[4420]
        print(f"  Worst diagnosis: Code {result_4420['worst_diagnosis_code']} ({result_4420['worst_diagnosis_label']})")
        print(f"  Date: {result_4420['worst_diagnosis_date']}")
        print(f"  Visit: {result_4420['worst_diagnosis_visit']}")
        
        # Check if this matches partner expectation
        if result_4420['worst_diagnosis_code'] == 3:
            print(f"  ✅ SUCCESS: Matches partner expectation (Code 3, AD)")
        else:
            print(f"  ❌ ISSUE: Does not match partner expectation (Code 3, AD)")
    else:
        print(f"  ❌ ERROR: Subject 941_S_4420 not found in results")
    
    print(f"\n📈 SUMMARY:")
    print(f"  Processed subjects: {len(all_subject_diagnoses)}")
    print(f"  Worst diagnoses calculated: {len(worst_diagnosis_map)}")
    
    return all_subject_diagnoses, worst_diagnosis_map

if __name__ == "__main__":
    all_diagnoses, worst_diagnoses = test_adni_diagnosis_processing()