#!/usr/bin/env python3
"""Test running just the ADNI processing part with our fixes."""

import sys
import os

# Add the pipeline directory to path
sys.path.append('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2')

from standardize_datasets import process_adni_data

def test_adni_pipeline_processing():
    """Test running the ADNI processing with our fixes."""
    
    print("🧪 TESTING ADNI PIPELINE PROCESSING WITH FIXES")
    print("=" * 50)
    
    # Set up paths
    input_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni_microarray'
    output_file = '/tmp/adni_test_output.h5ad'
    demographics_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/AutoSync/subject_demographics/PTDEMOG_25Apr2025.csv'
    diagnosis_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/AutoSync/subject_demographics/DXSUM_30Apr2025.csv'
    metadata_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json'
    
    print(f"📂 Input directory: {input_dir}")
    print(f"📂 Demographics file: {demographics_file}")
    print(f"📂 Diagnosis file: {diagnosis_file}")
    print(f"📂 Output file: {output_file}")
    
    # Check if files exist
    if not os.path.exists(diagnosis_file):
        print(f"❌ Diagnosis file not found: {diagnosis_file}")
        return False
    
    if not os.path.exists(demographics_file):
        print(f"❌ Demographics file not found: {demographics_file}")
        return False
        
    if not os.path.exists(input_dir):
        print(f"❌ Input directory not found: {input_dir}")
        return False
    
    print("✅ All input files found")
    
    try:
        # Run the ADNI processing
        print(f"\n🔄 Running ADNI data processing...")
        result = process_adni_data(
            input_dir=input_dir,
            output_file=output_file,
            adni_demographics_file=demographics_file,
            adni_diagnosis_file=diagnosis_file,
            metadata_dir=metadata_dir,
            gencode_map=None  # Will be handled by the function
        )
        
        if result is not None:
            print(f"✅ ADNI processing completed successfully!")
            print(f"📊 Result: {result.n_obs} samples × {result.n_vars} genes")
            
            # Check specific subject
            target_subject = '941_S_4420'
            subject_data = result.obs[result.obs['subject_id'] == target_subject]
            
            if len(subject_data) > 0:
                sample = subject_data.iloc[0]
                print(f"\n🎯 Subject {target_subject}:")
                print(f"  worst_diagnosis_code: {sample.get('worst_diagnosis_code', 'MISSING')}")
                print(f"  worst_diagnosis_label: {sample.get('worst_diagnosis_label', 'MISSING')}")
                print(f"  worst_diagnosis_date: {sample.get('worst_diagnosis_date', 'MISSING')}")
                print(f"  worst_diagnosis_visit: {sample.get('worst_diagnosis_visit', 'MISSING')}")
                
                # Check if it's Code 3 (AD)
                if sample.get('worst_diagnosis_code') == 3:
                    print(f"  ✅ SUCCESS: Shows Code 3 (AD) as expected!")
                    return True
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

if __name__ == "__main__":
    success = test_adni_pipeline_processing()
    if success:
        print(f"\n🎉 Test completed successfully!")
    else:
        print(f"\n💥 Test failed!")
        sys.exit(1)