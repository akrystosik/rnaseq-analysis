#!/usr/bin/env python3
"""Comprehensive test of ADNI diagnosis data consistency between H5AD and CSV export."""

import scanpy as sc
import pandas as pd

def test_diagnosis_consistency():
    """Test diagnosis data consistency across H5AD and CSV export."""
    
    print("🧪 COMPREHENSIVE ADNI DIAGNOSIS CONSISTENCY TEST")
    print("=" * 60)
    
    # Load H5AD file
    h5ad_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
    csv_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/adni_diagnosis_export.csv'
    
    print(f"📂 Loading H5AD: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    print(f"✅ H5AD loaded: {adata.n_obs} samples × {adata.n_vars} genes")
    
    print(f"📂 Loading CSV export: {csv_file}")
    csv_data = pd.read_csv(csv_file)
    print(f"✅ CSV loaded: {len(csv_data)} subjects")
    
    # Test cases - mix of different diagnosis codes and dates
    test_subjects = [
        '002_S_0413',  # Original problem case
        '002_S_0729',  # Original problem case  
        '002_S_1155',  # Original problem case
        '002_S_1261',  # Original problem case
        '002_S_1268',  # Original problem case
        '941_S_4420',  # Partner discrepancy case
        '002_S_0619',  # AD case
        '002_S_0782',  # MCI case
        '002_S_0559',  # CN case
        '941_S_1010',  # Different prefix test
    ]
    
    print(f"\n🎯 Testing {len(test_subjects)} subjects for consistency...")
    
    # Get H5AD diagnosis data
    h5ad_diagnosis = {}
    for subject in test_subjects:
        subject_data = adata.obs[adata.obs['subject_id'] == subject]
        if len(subject_data) > 0:
            sample = subject_data.iloc[0]
            h5ad_diagnosis[subject] = {
                'worst_diagnosis_code': sample.get('worst_diagnosis_code', 'MISSING'),
                'worst_diagnosis_label': sample.get('worst_diagnosis_label', 'MISSING'),
                'worst_diagnosis_date': sample.get('worst_diagnosis_date', 'MISSING'),
                'worst_diagnosis_visit': sample.get('worst_diagnosis_visit', 'MISSING'),
                'individual_diagnosis': sample.get('diagnosis', 'MISSING'),
                'individual_diagnosis_code': sample.get('diagnosis_code', 'MISSING')
            }
        else:
            h5ad_diagnosis[subject] = None
    
    # Get CSV diagnosis data  
    csv_diagnosis = {}
    for subject in test_subjects:
        subject_data = csv_data[csv_data['donor_id'] == subject]
        if len(subject_data) > 0:
            sample = subject_data.iloc[0]
            csv_diagnosis[subject] = {
                'most_severe_diagnosis_code': sample.get('most_severe_diagnosis_code', 'MISSING'),
                'most_severe_diagnosis': sample.get('most_severe_diagnosis', 'MISSING'),
                'diagnosis_date': sample.get('diagnosis_date', 'MISSING'),
                'diagnosis_visit': sample.get('diagnosis_visit', 'MISSING')
            }
        else:
            csv_diagnosis[subject] = None
    
    # Compare results
    print(f"\n📊 COMPARISON RESULTS:")
    print(f"{'Subject':<12} {'H5AD Code':<10} {'CSV Code':<9} {'H5AD Label':<25} {'CSV Label':<25} {'Match':<6}")
    print("-" * 95)
    
    matches = 0
    mismatches = 0
    missing_from_h5ad = 0
    missing_from_csv = 0
    
    for subject in test_subjects:
        h5ad_data = h5ad_diagnosis.get(subject)
        csv_data_subj = csv_diagnosis.get(subject)
        
        if h5ad_data is None and csv_data_subj is None:
            print(f"{subject:<12} {'MISSING':<10} {'MISSING':<9} {'MISSING':<25} {'MISSING':<25} {'N/A':<6}")
            continue
        elif h5ad_data is None:
            print(f"{subject:<12} {'MISSING':<10} {csv_data_subj['most_severe_diagnosis_code']:<9} {'MISSING':<25} {csv_data_subj['most_severe_diagnosis']:<25} {'❌':<6}")
            missing_from_h5ad += 1
            continue
        elif csv_data_subj is None:
            print(f"{subject:<12} {h5ad_data['worst_diagnosis_code']:<10} {'MISSING':<9} {h5ad_data['worst_diagnosis_label']:<25} {'MISSING':<25} {'❌':<6}")
            missing_from_csv += 1
            continue
        
        # Compare the data
        h5ad_code = h5ad_data['worst_diagnosis_code']
        csv_code = csv_data_subj['most_severe_diagnosis_code']
        h5ad_label = h5ad_data['worst_diagnosis_label']
        csv_label = csv_data_subj['most_severe_diagnosis']
        
        match = (h5ad_code == csv_code) and (h5ad_label == csv_label)
        match_symbol = '✅' if match else '❌'
        
        if match:
            matches += 1
        else:
            mismatches += 1
        
        print(f"{subject:<12} {h5ad_code:<10} {csv_code:<9} {h5ad_label:<25} {csv_label:<25} {match_symbol:<6}")
    
    print("-" * 95)
    print(f"\n📈 SUMMARY:")
    print(f"✅ Matches: {matches}")
    print(f"❌ Mismatches: {mismatches}")
    print(f"🔍 Missing from H5AD: {missing_from_h5ad}")
    print(f"🔍 Missing from CSV: {missing_from_csv}")
    
    # Test the original issue - individual vs worst diagnosis
    print(f"\n🔬 INDIVIDUAL vs WORST DIAGNOSIS TEST:")
    print(f"{'Subject':<12} {'Individual':<15} {'Worst':<25} {'Explanation'}")
    print("-" * 80)
    
    for subject in test_subjects[:5]:  # Test first 5
        h5ad_data = h5ad_diagnosis.get(subject)
        if h5ad_data:
            individual = h5ad_data['individual_diagnosis']
            worst = h5ad_data['worst_diagnosis_label']
            explanation = "Expected - worst from longitudinal" if individual != worst else "Same diagnosis"
            print(f"{subject:<12} {individual:<15} {worst:<25} {explanation}")
    
    # Test specific partner case
    print(f"\n🎯 PARTNER DISCREPANCY CASE (941_S_4420):")
    if '941_S_4420' in h5ad_diagnosis and h5ad_diagnosis['941_S_4420']:
        h5ad_4420 = h5ad_diagnosis['941_S_4420']
        csv_4420 = csv_diagnosis['941_S_4420']
        
        print(f"H5AD data:")
        print(f"  - Worst diagnosis code: {h5ad_4420['worst_diagnosis_code']}")
        print(f"  - Worst diagnosis label: {h5ad_4420['worst_diagnosis_label']}")
        print(f"  - Worst diagnosis date: {h5ad_4420['worst_diagnosis_date']}")
        print(f"  - Individual diagnosis: {h5ad_4420['individual_diagnosis']} (code: {h5ad_4420['individual_diagnosis_code']})")
        
        if csv_4420:
            print(f"CSV data:")
            print(f"  - Most severe code: {csv_4420['most_severe_diagnosis_code']}")
            print(f"  - Most severe label: {csv_4420['most_severe_diagnosis']}")
            print(f"  - Diagnosis date: {csv_4420['diagnosis_date']}")
            
            # Check if they match
            codes_match = h5ad_4420['worst_diagnosis_code'] == csv_4420['most_severe_diagnosis_code']
            labels_match = h5ad_4420['worst_diagnosis_label'] == csv_4420['most_severe_diagnosis']
            
            print(f"✅ Match: {codes_match and labels_match}")
            if not codes_match:
                print(f"❌ Code mismatch: H5AD={h5ad_4420['worst_diagnosis_code']} vs CSV={csv_4420['most_severe_diagnosis_code']}")
            if not labels_match:
                print(f"❌ Label mismatch: H5AD='{h5ad_4420['worst_diagnosis_label']}' vs CSV='{csv_4420['most_severe_diagnosis']}'")
    else:
        print("❌ Subject 941_S_4420 not found in H5AD data")
    
    # Overall assessment
    total_tests = matches + mismatches + missing_from_h5ad + missing_from_csv
    success_rate = (matches / total_tests * 100) if total_tests > 0 else 0
    
    print(f"\n🏆 OVERALL ASSESSMENT:")
    print(f"Success rate: {success_rate:.1f}% ({matches}/{total_tests})")
    
    if success_rate >= 90:
        print("🎉 EXCELLENT: Data consistency is very high")
    elif success_rate >= 70:
        print("⚠️  GOOD: Some issues but mostly consistent")
    else:
        print("❌ POOR: Significant consistency issues")
    
    return {
        'matches': matches,
        'mismatches': mismatches,
        'missing_from_h5ad': missing_from_h5ad,
        'missing_from_csv': missing_from_csv,
        'success_rate': success_rate
    }

if __name__ == "__main__":
    results = test_diagnosis_consistency()