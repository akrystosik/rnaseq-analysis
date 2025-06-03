#!/usr/bin/env python3
"""Debug subject 941_S_4420 diagnosis parsing."""

import pandas as pd

# Load diagnosis data
diagnosis_csv = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/AutoSync/subject_demographics/DXSUM_30Apr2025.csv'
dx_df = pd.read_csv(diagnosis_csv, low_memory=False)

# Filter for subject 4420
subject_data = dx_df[dx_df['RID'] == 4420].copy()
print(f"Found {len(subject_data)} records for RID 4420:")
print("\nColumns in data:")
print(list(dx_df.columns))

print(f"\nRaw data for subject 4420:")
for idx, row in subject_data.iterrows():
    print(f"  Visit: VISCODE='{row['VISCODE']}', VISCODE2='{row['VISCODE2']}'")
    print(f"  Date: {row['EXAMDATE']}")
    print(f"  Diagnosis: {row['DIAGNOSIS']}")
    print(f"  Raw row: {row.to_dict()}")
    print()

# Process like the export script
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

dx_df_clean = dx_df.dropna(subset=['RID', 'EXAMDATE', 'DIAGNOSIS'])
subject_4420 = dx_df_clean[dx_df_clean['RID'] == 4420]

print(f"\nAfter cleaning, found {len(subject_4420)} valid records:")
for idx, row in subject_4420.iterrows():
    print(f"  Date: {row['EXAMDATE']}, Diagnosis: {row['DIAGNOSIS']}")

print(f"\nProcessed valid diagnoses:")
valid_diagnoses = []
for _, row in subject_4420.sort_values(by='EXAMDATE').iterrows():
    diag_code = row['DIAGNOSIS']
    diag_label = diagnosis_mapping.get(diag_code, f"Unknown code: {diag_code}")
    exam_date = row['EXAMDATE']
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
    print(f"  {exam_date.date()}: Code {diag_code} ({diag_label}) at visit {visit_code}")

# Find worst diagnosis
diagnosis_severity_order = {1: 1, 2: 2, 3: 3, -4: 0}
worst_entry = max(valid_diagnoses, key=lambda x: diagnosis_severity_order.get(x['diagnosis_code'], 0))
print(f"\nWorst diagnosis: Code {worst_entry['diagnosis_code']} ({worst_entry['diagnosis_label']}) on {worst_entry['exam_date']}")