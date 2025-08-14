#!/usr/bin/env python3
"""
Investigate current tissue ontology annotations in partner's h5ad files
Identify which tissues need UBERON ID updates based on our validated corrections
"""

import pandas as pd
import anndata as ad
import json
from pathlib import Path

def investigate_h5ad_tissue_annotations(file_path, dataset_name):
    """Investigate tissue annotations in a single h5ad file"""
    print(f"\n📋 INVESTIGATING: {dataset_name}")
    print(f"File: {file_path}")
    
    try:
        # Load the h5ad file
        adata = ad.read_h5ad(file_path)
        print(f"Loaded AnnData: {adata.n_obs} samples × {adata.n_vars} genes")
        
        # Check available columns
        obs_columns = list(adata.obs.columns)
        tissue_related_cols = [col for col in obs_columns if 'tissue' in col.lower()]
        
        print(f"Available columns: {len(obs_columns)}")
        print(f"Tissue-related columns: {tissue_related_cols}")
        
        results = {
            'dataset': dataset_name,
            'file_path': str(file_path),
            'n_samples': adata.n_obs,
            'tissue_columns': tissue_related_cols,
            'tissue_annotations': {},
            'needs_correction': []
        }
        
        # Examine each tissue-related column
        for col in tissue_related_cols:
            print(f"\n🔍 Column: {col}")
            unique_values = adata.obs[col].unique()
            print(f"Unique values: {len(unique_values)}")
            
            # Store unique values and counts
            value_counts = adata.obs[col].value_counts()
            results['tissue_annotations'][col] = {
                'unique_count': len(unique_values),
                'values': list(unique_values)[:20],  # Limit for readability
                'value_counts': value_counts.head(10).to_dict()
            }
            
            # Show sample of values
            for val in unique_values[:10]:
                count = (adata.obs[col] == val).sum()
                print(f"  {val}: {count} samples")
            
            if len(unique_values) > 10:
                print(f"  ... and {len(unique_values) - 10} more values")
        
        return results
        
    except Exception as e:
        print(f"❌ Error loading {dataset_name}: {e}")
        return {
            'dataset': dataset_name,
            'error': str(e)
        }

def main():
    print("🔬 Partner H5AD Tissue Annotation Investigation")
    print("=" * 60)
    
    # Load our validated corrections
    corrections = {
        'brain - caudate (basal ganglia)': {
            'old': 'UBERON:0005382',
            'new': 'UBERON:0002420'
        },
        'brain - hippocampus': {
            'old': 'UBERON:0002310', 
            'new': 'UBERON:0002421'
        },
        'esophagus - muscularis': {
            'old': 'UBERON:0004648',
            'new': 'UBERON:0003832'
        },
        'heart - atrial appendage': {
            'old': 'UBERON:0006631',
            'new': 'UBERON:0006618'
        },
        'lymphoblast': {
            'old': 'CL:0000542',
            'new': 'CL:0017005'
        },
        'pbmc': {
            'old': 'UBERON:0000178',
            'new': 'CL:2000001'
        }
    }
    
    print(f"Validated corrections to apply: {len(corrections)}")
    
    # Partner's h5ad files
    data_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/run_20250502_211754")
    
    files_to_investigate = [
        (data_dir / "gtex_standardized_preprocessed.h5ad", "GTEx"),
        (data_dir / "adni_standardized_preprocessed.h5ad", "ADNI"),
        (data_dir / "encode_standardized_preprocessed.h5ad", "ENCODE"),
        (data_dir / "mage_standardized_preprocessed.h5ad", "MAGE")
    ]
    
    all_results = []
    
    # Investigate each file
    for file_path, dataset_name in files_to_investigate:
        if file_path.exists():
            result = investigate_h5ad_tissue_annotations(file_path, dataset_name)
            all_results.append(result)
        else:
            print(f"⚠️  File not found: {file_path}")
            all_results.append({
                'dataset': dataset_name,
                'error': 'File not found',
                'file_path': str(file_path)
            })
    
    print(f"\n" + "=" * 60)
    print("SUMMARY OF FINDINGS")
    
    # Identify potential corrections needed
    corrections_needed = []
    
    for result in all_results:
        if 'error' in result:
            continue
            
        dataset = result['dataset']
        print(f"\n📊 {dataset}:")
        print(f"  Samples: {result['n_samples']}")
        print(f"  Tissue columns: {result['tissue_columns']}")
        
        # Check if any tissue ontology values match our old incorrect IDs
        for col, col_data in result['tissue_annotations'].items():
            for tissue_name, correction in corrections.items():
                old_id = correction['old']
                new_id = correction['new']
                
                # Check if old incorrect ID is present
                if old_id in col_data['values']:
                    corrections_needed.append({
                        'dataset': dataset,
                        'column': col,
                        'tissue_name': tissue_name,
                        'old_id': old_id,
                        'new_id': new_id,
                        'sample_count': col_data['value_counts'].get(old_id, 0)
                    })
                    print(f"    🎯 CORRECTION NEEDED: {tissue_name}")
                    print(f"       Column: {col}")
                    print(f"       Old ID: {old_id} → New ID: {new_id}")
                    print(f"       Affected samples: {col_data['value_counts'].get(old_id, 0)}")
    
    # Save results
    investigation_results = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'validated_corrections': corrections,
        'h5ad_investigation': all_results,
        'corrections_needed': corrections_needed
    }
    
    output_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/debug_uberon/partner_h5ad_investigation.json"
    with open(output_file, 'w') as f:
        json.dump(investigation_results, f, indent=2)
    
    print(f"\n✅ Investigation complete!")
    print(f"📁 Results saved to: {output_file}")
    print(f"🔧 Corrections needed: {len(corrections_needed)}")
    
    if corrections_needed:
        print(f"\n🎯 SUMMARY OF CORRECTIONS NEEDED:")
        for correction in corrections_needed:
            print(f"  {correction['dataset']}: {correction['tissue_name']} ({correction['sample_count']} samples)")

if __name__ == "__main__":
    main()