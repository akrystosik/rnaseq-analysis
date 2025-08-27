#!/usr/bin/env python3
"""Update H5AD file with latest diagnosis data for subjects with expression data."""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import json

def update_h5ad_with_latest_diagnosis(
    h5ad_file=None, 
    diagnosis_csv=None,
    output_file=None,
    backup=True
):
    """Update H5AD file with latest diagnosis data."""
    
    # Default paths
    if h5ad_file is None:
        h5ad_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
    if diagnosis_csv is None:
        diagnosis_csv = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/AutoSync/subject_demographics/DXSUM_30Apr2025.csv'
    if output_file is None:
        output_file = h5ad_file  # Overwrite original
    
    print("🔄 UPDATING H5AD WITH LATEST DIAGNOSIS DATA")
    print("=" * 50)
    
    # Backup original file if requested
    if backup and output_file == h5ad_file:
        backup_file = f"{h5ad_file}.backup"
        print(f"📂 Creating backup: {backup_file}")
        import shutil
        shutil.copy2(h5ad_file, backup_file)
    
    # Load H5AD file
    print(f"📂 Loading H5AD: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    print(f"✅ Loaded: {adata.n_obs} samples × {adata.n_vars} genes")
    
    # Get list of expression subjects
    expression_subjects = set(adata.obs['subject_id'].unique())
    print(f"✅ Found {len(expression_subjects)} unique subjects with expression data")
    
    # Load and process diagnosis CSV
    print(f"📂 Loading diagnosis CSV: {diagnosis_csv}")
    dx_df = pd.read_csv(diagnosis_csv, low_memory=False)
    dx_df['EXAMDATE'] = pd.to_datetime(dx_df['EXAMDATE'], errors='coerce')
    dx_df['RID'] = pd.to_numeric(dx_df['RID'], errors='coerce').astype('Int64')
    dx_df['DIAGNOSIS'] = pd.to_numeric(dx_df['DIAGNOSIS'], errors='coerce')
    
    print(f"✅ Loaded {len(dx_df)} diagnosis records")
    
    # Diagnosis mapping
    diagnosis_mapping = {
        1: "Cognitively Normal", 
        2: "Mild Cognitive Impairment", 
        3: "Alzheimer's Disease",
        -4: "Missing/Unknown"
    }
    
    # Calculate worst diagnosis for each subject (same logic as export script)
    diagnosis_severity_order = {1: 1, 2: 2, 3: 3, -4: 0}
    
    dx_df_clean = dx_df.dropna(subset=['RID', 'DIAGNOSIS'])
    
    # Create mapping from subject_id to latest diagnosis
    subject_to_latest_diagnosis = {}
    
    for rid, group in dx_df_clean.groupby('RID'):
        # Use the actual PTID from the data
        subject_id = group['PTID'].iloc[0]
        
        # Only process if this subject has expression data
        if subject_id not in expression_subjects:
            continue
        
        # Get all valid diagnoses for this subject
        valid_diagnoses = []
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
            
            subject_to_latest_diagnosis[subject_id] = {
                'worst_diagnosis_code': worst_entry['diagnosis_code'],
                'worst_diagnosis_label': worst_entry['diagnosis_label'],
                'worst_diagnosis_visit': worst_entry['visit_code'],
                'worst_diagnosis_date': worst_entry['exam_date'],
                'all_diagnoses': valid_diagnoses  # Store full history
            }
    
    print(f"✅ Processed latest diagnosis for {len(subject_to_latest_diagnosis)} expression subjects")
    
    # Update H5AD observation data
    print("🔄 Updating H5AD observation data...")
    
    # Track changes
    changes = {'updated': 0, 'unchanged': 0, 'missing': 0}
    
    for idx, row in adata.obs.iterrows():
        subject_id = row['subject_id']
        
        if subject_id in subject_to_latest_diagnosis:
            latest_data = subject_to_latest_diagnosis[subject_id]
            
            # Check if data changed
            old_code = row.get('worst_diagnosis_code')
            new_code = latest_data['worst_diagnosis_code']
            
            if old_code != new_code:
                changes['updated'] += 1
                print(f"  📝 {subject_id}: {old_code} → {new_code} ({latest_data['worst_diagnosis_label']})")
            else:
                changes['unchanged'] += 1
            
            # Update the data
            adata.obs.loc[idx, 'worst_diagnosis_code'] = latest_data['worst_diagnosis_code']
            adata.obs.loc[idx, 'worst_diagnosis_label'] = latest_data['worst_diagnosis_label']
            adata.obs.loc[idx, 'worst_diagnosis_visit'] = latest_data['worst_diagnosis_visit']
            adata.obs.loc[idx, 'worst_diagnosis_date'] = latest_data['worst_diagnosis_date']
        else:
            changes['missing'] += 1
    
    # Update the longitudinal diagnosis data in adata.uns
    print("🔄 Updating longitudinal diagnosis data...")
    
    longitudinal_diagnoses = {}
    for subject_id, data in subject_to_latest_diagnosis.items():
        # Extract RID from subject_id
        if "_S_" in subject_id:
            rid = subject_id.split('_S_')[-1]
            longitudinal_diagnoses[rid] = json.dumps(data['all_diagnoses'])
    
    adata.uns['longitudinal_diagnoses_adni'] = longitudinal_diagnoses
    
    # Save updated H5AD file
    print(f"💾 Saving updated H5AD: {output_file}")
    adata.write_h5ad(output_file)
    
    # Summary
    print(f"\n📊 UPDATE SUMMARY:")
    print(f"✅ Updated subjects: {changes['updated']}")
    print(f"📌 Unchanged subjects: {changes['unchanged']}")
    print(f"❌ Missing diagnosis subjects: {changes['missing']}")
    print(f"🎯 Total subjects processed: {sum(changes.values())}")
    
    # Show examples of updated subjects
    if changes['updated'] > 0:
        print(f"\n🔍 SAMPLE UPDATES:")
        examples = 0
        for subject_id, latest_data in subject_to_latest_diagnosis.items():
            subject_obs = adata.obs[adata.obs['subject_id'] == subject_id]
            if len(subject_obs) > 0:
                old_code = subject_obs.iloc[0].get('worst_diagnosis_code')
                new_code = latest_data['worst_diagnosis_code']
                if old_code != new_code and examples < 5:
                    print(f"  {subject_id}: Code {old_code} → {new_code} ({latest_data['worst_diagnosis_date']})")
                    examples += 1
    
    return adata, changes

def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Update H5AD file with latest ADNI diagnosis data')
    parser.add_argument('--h5ad', help='Input H5AD file')
    parser.add_argument('--diagnosis-csv', help='Diagnosis CSV file')
    parser.add_argument('--output', help='Output H5AD file')
    parser.add_argument('--no-backup', action='store_true', help='Skip backup of original file')
    
    args = parser.parse_args()
    
    adata, changes = update_h5ad_with_latest_diagnosis(
        h5ad_file=args.h5ad,
        diagnosis_csv=args.diagnosis_csv,
        output_file=args.output,
        backup=not args.no_backup
    )
    
    print(f"\n🎉 H5AD update complete!")
    return changes

if __name__ == "__main__":
    main()