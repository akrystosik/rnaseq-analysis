#!/usr/bin/env python3
"""
Comprehensive review of all tissue and cell ontology mappings in partner's h5ad files
Compare against current corrected mapping and original backup to identify ALL needed updates
"""

import pandas as pd
import anndata as ad
import json
from pathlib import Path
import numpy as np

def load_mapping_files():
    """Load both current and backup mapping files"""
    current_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json"
    backup_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json_backup_20250814_184615/tissue_to_uberon.json"
    
    print("📁 Loading mapping files...")
    
    # Load current corrected mapping (lowercase keys)
    with open(current_file, 'r') as f:
        current_mapping = json.load(f)
    print(f"Current mapping: {len(current_mapping)} entries (corrected, lowercase)")
    
    # Load backup original mapping
    with open(backup_file, 'r') as f:
        backup_mapping = json.load(f)
    print(f"Backup mapping: {len(backup_mapping)} entries (original)")
    
    # Create normalized versions for comparison
    # Normalize backup mapping keys to lowercase for comparison
    backup_normalized = {k.lower(): v for k, v in backup_mapping.items()}
    
    return current_mapping, backup_normalized, backup_mapping

def analyze_h5ad_mappings(file_path, dataset_name, current_mapping, backup_mapping):
    """Analyze tissue/cell ontology mappings in a single h5ad file"""
    print(f"\n📋 ANALYZING: {dataset_name}")
    print(f"File: {file_path}")
    
    try:
        # Load the h5ad file
        adata = ad.read_h5ad(file_path)
        print(f"Loaded: {adata.n_obs} samples × {adata.n_vars} genes")
        
        # Get tissue-related columns
        tissue_cols = [col for col in adata.obs.columns if 'tissue' in col.lower()]
        print(f"Tissue columns: {tissue_cols}")
        
        analysis_results = {
            'dataset': dataset_name,
            'n_samples': adata.n_obs,
            'tissue_columns': tissue_cols,
            'mappings_analysis': {},
            'corrections_needed': [],
            'unmapped_tissues': []
        }
        
        # Analyze tissue and tissue_ontology columns
        if 'tissue' in adata.obs.columns and 'tissue_ontology' in adata.obs.columns:
            print("\n🔍 Analyzing tissue → tissue_ontology mappings:")
            
            # Get unique tissue to ontology pairs
            tissue_ontology_pairs = adata.obs[['tissue', 'tissue_ontology']].drop_duplicates()
            
            print(f"Unique tissue → ontology pairs: {len(tissue_ontology_pairs)}")
            
            for _, row in tissue_ontology_pairs.iterrows():
                tissue_name = str(row['tissue']).lower().strip()
                current_ontology = str(row['tissue_ontology']).strip()
                
                # Count samples with this mapping
                n_samples = ((adata.obs['tissue'].str.lower().str.strip() == tissue_name) & 
                           (adata.obs['tissue_ontology'].str.strip() == current_ontology)).sum()
                
                print(f"  '{tissue_name}' → '{current_ontology}' ({n_samples} samples)")
                
                # Check what the correct mapping should be
                correct_ontology = current_mapping.get(tissue_name, None)
                backup_ontology = backup_mapping.get(tissue_name, None)
                
                mapping_info = {
                    'tissue_name': tissue_name,
                    'current_h5ad_ontology': current_ontology,
                    'correct_ontology': correct_ontology,
                    'backup_ontology': backup_ontology,
                    'n_samples': int(n_samples)
                }
                
                # Determine if correction is needed
                if correct_ontology is None:
                    print(f"    ⚠️  Tissue '{tissue_name}' not found in mapping file")
                    analysis_results['unmapped_tissues'].append(mapping_info)
                elif current_ontology != correct_ontology:
                    if current_ontology == '' or pd.isna(current_ontology):
                        print(f"    🎯 NEEDS MAPPING: empty → {correct_ontology}")
                    else:
                        print(f"    🎯 NEEDS CORRECTION: {current_ontology} → {correct_ontology}")
                    
                    mapping_info['correction_type'] = 'mapping_needed' if current_ontology == '' else 'correction_needed'
                    analysis_results['corrections_needed'].append(mapping_info)
                else:
                    print(f"    ✅ Correct mapping")
                
                analysis_results['mappings_analysis'][tissue_name] = mapping_info
        
        # Check for tissue_ontology values that don't match any known correct mappings
        if 'tissue_ontology' in adata.obs.columns:
            print(f"\n🔍 Checking for invalid ontology IDs in tissue_ontology column:")
            
            unique_ontologies = adata.obs['tissue_ontology'].unique()
            all_correct_ontologies = set(current_mapping.values())
            
            for ontology_id in unique_ontologies:
                if ontology_id and ontology_id != '' and ontology_id not in all_correct_ontologies:
                    n_samples = (adata.obs['tissue_ontology'] == ontology_id).sum()
                    print(f"  ⚠️  Unknown ontology ID: '{ontology_id}' ({n_samples} samples)")
                    
                    # Try to find which tissue this ontology is associated with
                    associated_tissues = adata.obs[adata.obs['tissue_ontology'] == ontology_id]['tissue'].unique()
                    print(f"      Associated tissues: {list(associated_tissues)}")
        
        return analysis_results
        
    except Exception as e:
        print(f"❌ Error analyzing {dataset_name}: {e}")
        return {
            'dataset': dataset_name,
            'status': 'error',
            'error': str(e)
        }

def main():
    print("🔬 Comprehensive H5AD Mapping Review")
    print("=" * 60)
    
    # Load mapping files
    current_mapping, backup_normalized, backup_original = load_mapping_files()
    
    # Compare mapping files to show what changed
    print(f"\n📊 MAPPING FILE CHANGES:")
    changes_found = 0
    for tissue, current_id in current_mapping.items():
        backup_id = backup_normalized.get(tissue)
        if backup_id != current_id:
            changes_found += 1
            print(f"  {tissue}: {backup_id} → {current_id}")
    
    print(f"Total mapping changes made: {changes_found}")
    
    # Partner's h5ad files
    data_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/run_20250502_211754")
    
    files_to_analyze = [
        (data_dir / "gtex_standardized_preprocessed.h5ad", "GTEx"),
        (data_dir / "adni_standardized_preprocessed.h5ad", "ADNI"),
        (data_dir / "encode_standardized_preprocessed.h5ad", "ENCODE"),
        (data_dir / "mage_standardized_preprocessed.h5ad", "MAGE")
    ]
    
    # Analyze each file
    all_analyses = []
    total_corrections_needed = 0
    
    for file_path, dataset_name in files_to_analyze:
        if file_path.exists():
            analysis = analyze_h5ad_mappings(file_path, dataset_name, current_mapping, backup_normalized)
            all_analyses.append(analysis)
            
            if 'corrections_needed' in analysis:
                dataset_corrections = len(analysis['corrections_needed'])
                total_corrections_needed += dataset_corrections
                print(f"  Corrections needed for {dataset_name}: {dataset_corrections}")
        else:
            print(f"⚠️  File not found: {file_path}")
    
    # Comprehensive summary
    print(f"\n" + "=" * 60)
    print("COMPREHENSIVE SUMMARY")
    print(f"Total datasets analyzed: {len(all_analyses)}")
    print(f"Total correction items needed: {total_corrections_needed}")
    
    # Detailed summary by dataset
    for analysis in all_analyses:
        if 'corrections_needed' in analysis and analysis['corrections_needed']:
            print(f"\n📋 {analysis['dataset']} - {len(analysis['corrections_needed'])} corrections needed:")
            for correction in analysis['corrections_needed']:
                tissue = correction['tissue_name']
                current_ont = correction['current_h5ad_ontology']
                correct_ont = correction['correct_ontology']
                samples = correction['n_samples']
                
                if current_ont == '':
                    print(f"  🎯 '{tissue}': ADD mapping → {correct_ont} ({samples} samples)")
                else:
                    print(f"  🎯 '{tissue}': {current_ont} → {correct_ont} ({samples} samples)")
    
    # Summary of correction types
    correction_types = {}
    for analysis in all_analyses:
        for correction in analysis.get('corrections_needed', []):
            old_id = correction['current_h5ad_ontology']
            new_id = correction['correct_ontology']
            key = f"{old_id} → {new_id}"
            
            if key not in correction_types:
                correction_types[key] = {'count': 0, 'samples': 0, 'datasets': set()}
            
            correction_types[key]['count'] += 1
            correction_types[key]['samples'] += correction['n_samples']
            correction_types[key]['datasets'].add(analysis['dataset'])
    
    print(f"\n📊 CORRECTION SUMMARY BY TYPE:")
    for correction_type, info in correction_types.items():
        datasets = ', '.join(info['datasets'])
        print(f"  {correction_type}: {info['samples']} samples across {len(info['datasets'])} datasets ({datasets})")
    
    # Save comprehensive analysis
    comprehensive_results = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'mapping_comparison': {
            'current_mapping_entries': len(current_mapping),
            'backup_mapping_entries': len(backup_normalized),
            'changes_made': changes_found
        },
        'h5ad_analyses': all_analyses,
        'correction_summary': {
            'total_corrections_needed': total_corrections_needed,
            'correction_types': correction_types
        }
    }
    
    output_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/debug_uberon/comprehensive_h5ad_mapping_review.json"
    with open(output_file, 'w') as f:
        json.dump(comprehensive_results, f, indent=2, default=str)
    
    print(f"\n✅ Comprehensive analysis complete!")
    print(f"📁 Results saved to: {output_file}")

if __name__ == "__main__":
    main()