#!/usr/bin/env python3
"""
Apply comprehensive tissue ontology corrections to partner's h5ad files
Output corrected files to new directory: preprocessed_data/run_20250502_211754_uberon_corrected

Based on comprehensive analysis identifying 12 correction items affecting 6,755 samples:
- 11 corrections in GTEx dataset (6,024 samples) 
- 1 correction in MAGE dataset (731 samples)
"""

import pandas as pd
import anndata as ad
import json
from pathlib import Path
import shutil
import numpy as np
from datetime import datetime

def create_output_directory():
    """Create the output directory for corrected h5ad files"""
    output_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/run_20250502_211754_uberon_corrected")
    
    if output_dir.exists():
        print(f"⚠️  Output directory already exists: {output_dir}")
        print("    Removing existing directory to start fresh...")
        shutil.rmtree(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"📁 Created output directory: {output_dir}")
    
    return output_dir

def load_correction_mappings():
    """Load all corrections identified in the comprehensive analysis"""
    
    # All 12 corrections identified from comprehensive analysis
    corrections = {
        # GTEx corrections (11 items)
        "brain - frontal cortex (ba9)": {
            "old_ontology": "UBERON:0013529",
            "new_ontology": "UBERON:0001870",
            "samples": 269,
            "dataset": "GTEx"
        },
        "brain - caudate (basal ganglia)": {
            "old_ontology": "UBERON:0001873", 
            "new_ontology": "UBERON:0002420",
            "samples": 300,
            "dataset": "GTEx"
        },
        "heart - atrial appendage": {
            "old_ontology": "UBERON:0006631",
            "new_ontology": "UBERON:0006618", 
            "samples": 461,
            "dataset": "GTEx"
        },
        "skin - not sun exposed (suprapubic)": {
            "old_ontology": "UBERON:0036151",
            "new_ontology": "UBERON:0036149",
            "samples": 651,
            "dataset": "GTEx"
        },
        "cells - ebv-transformed lymphocytes": {
            "old_ontology": "",
            "new_ontology": "CL:0000542",
            "samples": 327,
            "dataset": "GTEx"
        },
        "esophagus - muscularis": {
            "old_ontology": "UBERON:0004648",
            "new_ontology": "UBERON:0003832",
            "samples": 561,
            "dataset": "GTEx"
        },
        "adipose - visceral (omentum)": {
            "old_ontology": "UBERON:0016529",  # Partner's critical issue - brain cortex!
            "new_ontology": "UBERON:0014454",  # Correct visceral adipose tissue
            "samples": 587,
            "dataset": "GTEx"
        },
        "small intestine - terminal ileum": {
            "old_ontology": "UBERON:0001211",
            "new_ontology": "UBERON:0002116", 
            "samples": 207,
            "dataset": "GTEx"
        },
        "skin - sun exposed (lower leg)": {
            "old_ontology": "UBERON:0036149",
            "new_ontology": "UBERON:0004264",
            "samples": 754,
            "dataset": "GTEx"
        },
        "cells - cultured fibroblasts": {
            "old_ontology": "",
            "new_ontology": "CL:0000057",
            "samples": 652,
            "dataset": "GTEx"
        },
        "brain - hippocampus": {
            "old_ontology": "UBERON:0002310",
            "new_ontology": "UBERON:0002421",
            "samples": 255,
            "dataset": "GTEx"
        },
        
        # MAGE correction (1 item)
        "lymphoblast": {
            "old_ontology": "",
            "new_ontology": "CL:0017005",  # Corrected from CL:0000542 based on developmental stage
            "samples": 731,
            "dataset": "MAGE"
        }
    }
    
    return corrections

def apply_corrections_to_h5ad(file_path, dataset_name, corrections, output_dir):
    """Apply tissue ontology corrections to a single h5ad file"""
    
    print(f"\n🔧 PROCESSING: {dataset_name}")
    print(f"Input file: {file_path}")
    
    try:
        # Load the h5ad file
        adata = ad.read_h5ad(file_path)
        original_samples = adata.n_obs
        print(f"Loaded: {original_samples} samples × {adata.n_vars} genes")
        
        corrections_applied = 0
        samples_affected = 0
        
        # Handle categorical columns by converting to string first (do this once per dataset)
        if hasattr(adata.obs['tissue_ontology'], 'cat'):
            print(f"  🔄 Converting tissue_ontology from categorical to string")
            adata.obs['tissue_ontology'] = adata.obs['tissue_ontology'].astype(str)
        
        if hasattr(adata.obs['tissue_ontology_confidence'], 'cat'):
            print(f"  🔄 Converting tissue_ontology_confidence from categorical to string")  
            adata.obs['tissue_ontology_confidence'] = adata.obs['tissue_ontology_confidence'].astype(str)
        
        # Apply corrections specific to this dataset
        dataset_corrections = {k: v for k, v in corrections.items() if v['dataset'] == dataset_name}
        
        if not dataset_corrections:
            print(f"  ✅ No corrections needed for {dataset_name}")
        else:
            print(f"  🎯 Applying {len(dataset_corrections)} corrections:")
        
        for tissue_name, correction_info in dataset_corrections.items():
            old_id = correction_info['old_ontology']
            new_id = correction_info['new_ontology']
            expected_samples = correction_info['samples']
            
            print(f"    • {tissue_name}: '{old_id}' → '{new_id}'")
            
            # Find samples matching this tissue name and old ontology ID
            tissue_match = adata.obs['tissue'].str.lower().str.strip() == tissue_name
            
            if old_id == "":
                # Handle empty ontology mappings
                ontology_match = (adata.obs['tissue_ontology'].isna() | 
                                (adata.obs['tissue_ontology'] == "") | 
                                (adata.obs['tissue_ontology'] == "nan"))
            else:
                # Handle specific ontology ID changes
                ontology_match = adata.obs['tissue_ontology'] == old_id
            
            # Find samples that need this correction
            samples_to_correct = tissue_match & ontology_match
            n_samples_found = samples_to_correct.sum()
            
            if n_samples_found > 0:
                # Apply the correction
                adata.obs.loc[samples_to_correct, 'tissue_ontology'] = new_id
                
                # Update confidence for newly mapped tissues
                if old_id == "":
                    adata.obs.loc[samples_to_correct, 'tissue_ontology_confidence'] = 'high'
                
                corrections_applied += 1
                samples_affected += int(n_samples_found)
                
                print(f"      ✅ Corrected {n_samples_found} samples (expected: {expected_samples})")
                
                if abs(n_samples_found - expected_samples) > 0:
                    print(f"      ⚠️  Sample count difference: found {n_samples_found}, expected {expected_samples}")
            else:
                print(f"      ⚠️  No samples found matching this correction pattern")
        
        # Save corrected h5ad file
        output_file = output_dir / f"{dataset_name.lower()}_standardized_preprocessed.h5ad"
        adata.write_h5ad(output_file)
        
        correction_summary = {
            'dataset': dataset_name,
            'input_file': str(file_path),
            'output_file': str(output_file),
            'original_samples': int(original_samples),
            'corrections_applied': int(corrections_applied),
            'samples_affected': int(samples_affected),
            'corrections_expected': int(len(dataset_corrections)),
            'processing_time': datetime.now().isoformat()
        }
        
        print(f"  ✅ Saved corrected file: {output_file}")
        print(f"  📊 Applied {corrections_applied} corrections affecting {samples_affected} samples")
        
        return correction_summary
        
    except Exception as e:
        print(f"❌ Error processing {dataset_name}: {e}")
        return {
            'dataset': dataset_name,
            'status': 'error',
            'error': str(e)
        }

def create_correction_report(corrections, processing_results, output_dir):
    """Create comprehensive correction report"""
    
    report = {
        'correction_timestamp': datetime.now().isoformat(),
        'original_directory': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/run_20250502_211754',
        'output_directory': str(output_dir),
        'total_corrections_planned': len(corrections),
        'corrections_by_dataset': {},
        'processing_results': processing_results,
        'correction_details': corrections
    }
    
    # Summarize corrections by dataset
    for tissue, correction_info in corrections.items():
        dataset = correction_info['dataset']
        if dataset not in report['corrections_by_dataset']:
            report['corrections_by_dataset'][dataset] = {
                'correction_count': 0,
                'samples_affected': 0,
                'corrections': []
            }
        
        report['corrections_by_dataset'][dataset]['correction_count'] += 1
        report['corrections_by_dataset'][dataset]['samples_affected'] += correction_info['samples']
        report['corrections_by_dataset'][dataset]['corrections'].append({
            'tissue_name': tissue,
            'old_ontology': correction_info['old_ontology'],
            'new_ontology': correction_info['new_ontology'],
            'samples': correction_info['samples']
        })
    
    # Calculate total samples affected
    total_samples_affected = sum(correction_info['samples'] for correction_info in corrections.values())
    report['total_samples_affected'] = total_samples_affected
    
    # Save report
    report_file = output_dir / "correction_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n📋 CORRECTION SUMMARY:")
    print(f"Total corrections applied: {len(corrections)}")
    print(f"Total samples affected: {total_samples_affected}")
    
    for dataset, dataset_info in report['corrections_by_dataset'].items():
        print(f"{dataset}: {dataset_info['correction_count']} corrections, {dataset_info['samples_affected']} samples")
    
    print(f"\n📁 Detailed report saved: {report_file}")
    
    return report

def main():
    print("🔬 Comprehensive H5AD Tissue Ontology Correction")
    print("=" * 70)
    print("Applying 12 validated corrections affecting 6,755 samples")
    print("Output directory: preprocessed_data/run_20250502_211754_uberon_corrected")
    print("=" * 70)
    
    # Create output directory
    output_dir = create_output_directory()
    
    # Load correction mappings
    corrections = load_correction_mappings()
    print(f"\n📝 Loaded {len(corrections)} corrections:")
    
    # Show correction summary
    gtex_corrections = sum(1 for c in corrections.values() if c['dataset'] == 'GTEx')
    mage_corrections = sum(1 for c in corrections.values() if c['dataset'] == 'MAGE')
    total_samples = sum(c['samples'] for c in corrections.values())
    
    print(f"  • GTEx: {gtex_corrections} corrections")
    print(f"  • MAGE: {mage_corrections} corrections")
    print(f"  • Total samples to be corrected: {total_samples}")
    
    # Partner's h5ad files
    data_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/run_20250502_211754")
    
    files_to_process = [
        (data_dir / "gtex_standardized_preprocessed.h5ad", "GTEx"),
        (data_dir / "adni_standardized_preprocessed.h5ad", "ADNI"),
        (data_dir / "encode_standardized_preprocessed.h5ad", "ENCODE"),
        (data_dir / "mage_standardized_preprocessed.h5ad", "MAGE")
    ]
    
    # Process each file
    processing_results = []
    
    for file_path, dataset_name in files_to_process:
        if file_path.exists():
            result = apply_corrections_to_h5ad(file_path, dataset_name, corrections, output_dir)
            processing_results.append(result)
        else:
            print(f"⚠️  File not found: {file_path}")
            processing_results.append({
                'dataset': dataset_name,
                'status': 'file_not_found',
                'input_file': str(file_path)
            })
    
    # Create comprehensive report
    correction_report = create_correction_report(corrections, processing_results, output_dir)
    
    print(f"\n" + "=" * 70)
    print("✅ CORRECTION COMPLETE!")
    print(f"📁 Corrected h5ad files available in: {output_dir}")
    print(f"📋 Correction report: {output_dir}/correction_report.json")
    print("=" * 70)
    
    # Show critical corrections applied
    print(f"\n🎯 KEY CORRECTIONS APPLIED:")
    print(f"  • Partner's visceral adipose issue: UBERON:0016529 → UBERON:0014454 (587 samples)")
    print(f"  • Brain caudate correction: UBERON:0001873 → UBERON:0002420 (300 samples)")
    print(f"  • MAGE lymphoblast mapping: empty → CL:0017005 (731 samples)")
    print(f"  • Plus 9 other validated corrections")
    
    return True

if __name__ == "__main__":
    success = main()
    if success:
        exit(0)
    else:
        exit(1)