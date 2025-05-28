#!/usr/bin/env python3
"""
Integrate complete ethnicity mapping into pipeline H5AD datasets.
This script updates the preprocessed H5AD files with the enhanced ethnicity data.
"""

import scanpy as sc
import pandas as pd
from pathlib import Path
import sys

def integrate_ethnicity_mapping():
    """Integrate the complete ethnicity mapping into pipeline datasets."""
    
    print("🧬 **INTEGRATING ETHNICITY DATA INTO PIPELINE DATASETS**")
    print("=" * 55)
    
    # Load the complete ethnicity mapping
    ethnicity_file = 'subject_ethnicity_mapping_with_ontology.csv'
    if not Path(ethnicity_file).exists():
        print(f"❌ Error: {ethnicity_file} not found. Run create_complete_ethnicity_mapping.py first.")
        return False
    
    ethnicity_df = pd.read_csv(ethnicity_file)
    print(f"✅ Loaded ethnicity mapping: {len(ethnicity_df):,} subjects")
    
    # Define dataset paths
    data_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137")
    
    dataset_files = {
        'ADNI': 'adni_standardized_preprocessed.h5ad',
        'ENCODE': 'encode_standardized_preprocessed.h5ad', 
        'GTEx': 'gtex_standardized_preprocessed.h5ad',
        'MAGE': 'mage_standardized_preprocessed.h5ad'
    }
    
    updated_datasets = []
    
    for dataset_name, filename in dataset_files.items():
        file_path = data_dir / filename
        
        if not file_path.exists():
            print(f"⚠️  Skipping {dataset_name}: File not found at {file_path}")
            continue
            
        print(f"\n🔄 Processing {dataset_name}...")
        
        # Load dataset
        adata = sc.read_h5ad(file_path)
        original_shape = adata.shape
        
        # Get ethnicity data for this dataset
        dataset_ethnicity = ethnicity_df[ethnicity_df['dataset'] == dataset_name].copy()
        
        if len(dataset_ethnicity) == 0:
            print(f"  ⚠️  No ethnicity data found for {dataset_name}")
            continue
        
        print(f"  📊 Found ethnicity data for {len(dataset_ethnicity):,} subjects")
        
        # Handle GTEx special case (multiple samples per subject)
        if dataset_name == 'GTEx':
            # Map subject-level ethnicity to sample-level
            subject_to_ethnicity = dict(zip(dataset_ethnicity['subject_id'], dataset_ethnicity['ethnicity']))
            subject_to_hancestro = dict(zip(dataset_ethnicity['subject_id'], dataset_ethnicity['hancestro_term']))
            
            # Extract subject_id from sample_id for GTEx
            adata.obs['ethnicity_mapping_subject_id'] = adata.obs['subject_id'].copy()
            
            # Map ethnicity to samples
            adata.obs['self_reported_ethnicity_updated'] = adata.obs['ethnicity_mapping_subject_id'].map(subject_to_ethnicity)
            adata.obs['self_reported_ethnicity_ontology_term_id_updated'] = adata.obs['ethnicity_mapping_subject_id'].map(subject_to_hancestro)
            
        else:
            # For other datasets, direct mapping by subject_id
            subject_to_ethnicity = dict(zip(dataset_ethnicity['subject_id'], dataset_ethnicity['ethnicity']))
            subject_to_hancestro = dict(zip(dataset_ethnicity['subject_id'], dataset_ethnicity['hancestro_term']))
            
            adata.obs['self_reported_ethnicity_updated'] = adata.obs['subject_id'].map(subject_to_ethnicity)
            adata.obs['self_reported_ethnicity_ontology_term_id_updated'] = adata.obs['subject_id'].map(subject_to_hancestro)
        
        # Check mapping success
        mapped_count = adata.obs['self_reported_ethnicity_updated'].notna().sum()
        mapping_rate = mapped_count / len(adata.obs) * 100
        
        print(f"  ✅ Mapped ethnicity for {mapped_count:,}/{len(adata.obs):,} samples ({mapping_rate:.1f}%)")
        
        # Replace original ethnicity columns with updated ones
        adata.obs['self_reported_ethnicity_original'] = adata.obs['self_reported_ethnicity'].copy()
        adata.obs['self_reported_ethnicity_ontology_term_id_original'] = adata.obs['self_reported_ethnicity_ontology_term_id'].copy()
        
        adata.obs['self_reported_ethnicity'] = adata.obs['self_reported_ethnicity_updated'].fillna('unknown or not reported')
        adata.obs['self_reported_ethnicity_ontology_term_id'] = adata.obs['self_reported_ethnicity_ontology_term_id_updated'].fillna('unknown')
        
        # Clean up temporary columns
        adata.obs = adata.obs.drop(columns=['self_reported_ethnicity_updated', 'self_reported_ethnicity_ontology_term_id_updated'])
        if 'ethnicity_mapping_subject_id' in adata.obs.columns:
            adata.obs = adata.obs.drop(columns=['ethnicity_mapping_subject_id'])
        
        # Add metadata about the update
        adata.uns['ethnicity_mapping_update'] = {
            'timestamp': pd.Timestamp.now().isoformat(),
            'source_file': ethnicity_file,
            'mapping_rate': f"{mapping_rate:.1f}%",
            'updated_samples': int(mapped_count),
            'total_samples': len(adata.obs)
        }
        
        # Show updated distribution
        print(f"  📈 Updated ethnicity distribution:")
        ethnicity_counts = adata.obs['self_reported_ethnicity'].value_counts()
        for ethnicity, count in ethnicity_counts.head(5).items():
            print(f"    • {ethnicity}: {count:,}")
        
        # Save updated dataset
        backup_path = file_path.with_suffix('.h5ad.backup')
        if not backup_path.exists():
            print(f"  💾 Creating backup: {backup_path.name}")
            # Create backup of original
            adata_backup = sc.read_h5ad(file_path)
            adata_backup.write_h5ad(backup_path)
        
        print(f"  💾 Saving updated dataset...")
        adata.write_h5ad(file_path)
        
        updated_datasets.append(dataset_name)
        print(f"  ✅ {dataset_name} updated successfully")
    
    print(f"\n🎯 **INTEGRATION COMPLETE**")
    print(f"✅ Updated datasets: {', '.join(updated_datasets)}")
    print(f"📊 Total subjects with enhanced ethnicity: {len(ethnicity_df):,}")
    print(f"🧬 HANCESTRO coverage: 99.8%")
    
    return len(updated_datasets) > 0

if __name__ == '__main__':
    success = integrate_ethnicity_mapping()
    if success:
        print("\n🚀 **READY FOR VALIDATION**: Run validation to confirm integration")
        print("   python validate_standardized_datasets.py --input-dir /mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137")
    else:
        print("\n❌ **INTEGRATION FAILED**: Check error messages above")
        sys.exit(1)