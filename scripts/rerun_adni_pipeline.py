#!/usr/bin/env python3
"""
Re-run ADNI portion of pipeline with diagnosis fixes to regenerate H5AD.

This script runs the ADNI data processing with our fixes to create an updated
H5AD file containing the corrected diagnosis data.
"""

import sys
import os
import pandas as pd
from pathlib import Path

# Add the pipeline directory to path
sys.path.append('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2')

from standardize_datasets import process_adni_data

def rerun_adni_pipeline():
    """Re-run ADNI pipeline processing with diagnosis fixes."""
    
    print("🔄 RE-RUNNING ADNI PIPELINE WITH DIAGNOSIS FIXES")
    print("=" * 60)
    
    # Set up paths
    input_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni_microarray'
    output_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/adni_standardized_updated.h5ad'
    demographics_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/AutoSync/subject_demographics/PTDEMOG_25Apr2025.csv'
    diagnosis_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/AutoSync/subject_demographics/DXSUM_30Apr2025.csv'
    metadata_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json'
    
    print(f"📂 Input directory: {input_dir}")
    print(f"📂 Demographics file: {demographics_file}")
    print(f"📂 Diagnosis file: {diagnosis_file}")
    print(f"📂 Output file: {output_file}")
    
    # Check if files exist
    files_to_check = [
        (diagnosis_file, "Diagnosis file"),
        (demographics_file, "Demographics file"),
        (input_dir, "Input directory")
    ]
    
    for file_path, description in files_to_check:
        if not os.path.exists(file_path):
            print(f"❌ {description} not found: {file_path}")
            return False
            
    print("✅ All input files found")
    
    # Load the gencode mapping (required for processing)
    from rnaseq_utils import load_gencode_mapping
    
    print(f"\n🧬 Loading gene mapping...")
    try:
        gencode_map = load_gencode_mapping()
        if gencode_map is None:
            print("❌ Failed to load gencode mapping")
            return False
        print(f"✅ Loaded gencode mapping with {len(gencode_map)} entries")
    except Exception as e:
        print(f"❌ Error loading gencode mapping: {e}")
        return False
    
    try:
        # Run the ADNI processing with fixes
        print(f"\n🔄 Running ADNI data processing with diagnosis fixes...")
        result = process_adni_data(
            input_dir=input_dir,
            output_file=output_file,
            adni_demographics_file=demographics_file,
            adni_diagnosis_file=diagnosis_file,
            metadata_dir=metadata_dir,
            gencode_map=gencode_map
        )
        
        if result is not None:
            print(f"✅ ADNI processing completed successfully!")
            print(f"📊 Result: {result.n_obs} samples × {result.n_vars} genes")
            
            # Check specific subject 941_S_4420
            target_subject = '941_S_4420'
            subject_data = result.obs[result.obs['subject_id'] == target_subject]
            
            if len(subject_data) > 0:
                sample = subject_data.iloc[0]
                print(f"\n🎯 VERIFICATION - Subject {target_subject}:")
                print(f"  worst_diagnosis_code: {sample.get('worst_diagnosis_code', 'MISSING')}")
                print(f"  worst_diagnosis_label: {sample.get('worst_diagnosis_label', 'MISSING')}")
                print(f"  worst_diagnosis_date: {sample.get('worst_diagnosis_date', 'MISSING')}")
                print(f"  worst_diagnosis_visit: {sample.get('worst_diagnosis_visit', 'MISSING')}")
                
                # Check if it's Code 3 (AD) as expected
                if sample.get('worst_diagnosis_code') == 3:
                    print(f"  ✅ SUCCESS: Shows Code 3 (AD) as expected!")
                    
                    # Verify date is from 2023
                    diagnosis_date = str(sample.get('worst_diagnosis_date', ''))
                    if '2023' in diagnosis_date:
                        print(f"  ✅ SUCCESS: Date is from 2023 as expected!")
                    else:
                        print(f"  ⚠️  NOTE: Date is {diagnosis_date}, expected 2023")
                        
                    return output_file
                else:
                    print(f"  ❌ ISSUE: Expected Code 3, got {sample.get('worst_diagnosis_code')}")
                    return False
            else:
                print(f"❌ Subject {target_subject} not found in results")
                return False
        else:
            print(f"❌ ADNI processing failed")
            return False
            
    except Exception as e:
        print(f"❌ Error during ADNI processing: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_updated_h5ad(h5ad_file):
    """Test the updated H5AD file for diagnosis data."""
    
    print(f"\n🧪 TESTING UPDATED H5AD FILE")
    print("=" * 40)
    
    try:
        import scanpy as sc
        
        print(f"📂 Loading updated H5AD: {h5ad_file}")
        adata = sc.read_h5ad(h5ad_file)
        
        print(f"✅ Loaded: {adata.n_obs} samples × {adata.n_vars} genes")
        
        # Test key subjects
        test_subjects = ['941_S_4420', '002_S_0413', '002_S_0729']
        
        print(f"\n🎯 Testing key subjects:")
        for subject_id in test_subjects:
            subject_data = adata.obs[adata.obs['subject_id'] == subject_id]
            
            if len(subject_data) > 0:
                sample = subject_data.iloc[0]
                code = sample.get('worst_diagnosis_code', 'N/A')
                label = sample.get('worst_diagnosis_label', 'N/A')
                date = sample.get('worst_diagnosis_date', 'N/A')
                
                print(f"  {subject_id}: Code {code} ({label}) from {date}")
                
                if subject_id == '941_S_4420' and code == 3:
                    print(f"    ✅ 941_S_4420 correctly shows Code 3 (AD)!")
            else:
                print(f"  {subject_id}: ❌ Not found")
        
        return True
        
    except Exception as e:
        print(f"❌ Error testing H5AD: {e}")
        return False

if __name__ == "__main__":
    print("🚀 Starting ADNI pipeline re-run...")
    
    # Re-run ADNI processing
    result_file = rerun_adni_pipeline()
    
    if result_file:
        print(f"\n✅ ADNI pipeline re-run completed successfully!")
        print(f"📁 Updated H5AD file: {result_file}")
        
        # Test the updated file
        if test_updated_h5ad(result_file):
            print(f"\n🎉 Updated H5AD file verified successfully!")
            print(f"\n🔄 Ready for next steps:")
            print(f"  1. Export diagnosis CSV from updated H5AD")
            print(f"  2. Push changes to GitHub")
            print(f"  3. Share with partner")
        else:
            print(f"\n❌ H5AD verification failed")
            sys.exit(1)
    else:
        print(f"\n💥 ADNI pipeline re-run failed!")
        sys.exit(1)