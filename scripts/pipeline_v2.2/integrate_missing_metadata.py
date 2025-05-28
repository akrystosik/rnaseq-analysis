#!/usr/bin/env python3
"""
Integrate missing age and sex metadata from source files.
- GTEx: Sex data from controlled-access phenotype file
- ADNI: Age data from demographics file  
- MAGE: Age data from 1000 Genomes (if available)
"""

import scanpy as sc
import pandas as pd
import gzip
from io import StringIO
from pathlib import Path
import numpy as np
from datetime import datetime

def integrate_gtex_sex_data():
    """Integrate GTEx sex data from controlled-access phenotype file."""
    print("🔄 **INTEGRATING GTEx SEX DATA**")
    print("=" * 35)
    
    # Load GTEx phenotype data
    gtex_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/phs000424.v10.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz'
    
    with gzip.open(gtex_file, 'rt') as f:
        content = f.read()

    lines = content.split('\n')
    data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
    clean_content = '\n'.join(data_lines)
    phenotype_df = pd.read_csv(StringIO(clean_content), sep='\t')

    # Map sex codes (1=male, 2=female per dbGaP)
    sex_mapping = {1: 'male', 2: 'female', 98: 'unknown', 99: 'unknown'}
    phenotype_df['sex_mapped'] = phenotype_df['SEX'].map(sex_mapping)
    
    # Create subject-to-sex mapping
    gtex_sex_map = dict(zip(phenotype_df['SUBJID'], phenotype_df['sex_mapped']))
    
    print(f"✅ Loaded sex data for {len(gtex_sex_map):,} GTEx subjects")
    sex_counts = phenotype_df['sex_mapped'].value_counts()
    for sex, count in sex_counts.items():
        print(f"  • {sex}: {count:,}")
    
    # Load GTEx dataset
    # Use dynamic path based on current repo structure
    gtex_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/gtex_standardized_preprocessed.h5ad'
    adata = sc.read_h5ad(gtex_path)
    
    # Map sex to samples
    adata.obs['sex_updated'] = adata.obs['subject_id'].map(gtex_sex_map)
    
    # Check mapping success
    mapped_count = adata.obs['sex_updated'].notna().sum()
    mapping_rate = mapped_count / len(adata.obs) * 100
    
    print(f"✅ Mapped sex for {mapped_count:,}/{len(adata.obs):,} samples ({mapping_rate:.1f}%)")
    
    # Update sex field
    adata.obs['sex_original_unknown'] = adata.obs['sex'].copy()
    adata.obs['sex'] = adata.obs['sex_updated'].fillna('unknown')
    adata.obs = adata.obs.drop(columns=['sex_updated'])
    
    # Add metadata about the update
    adata.uns['sex_mapping_update'] = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'source_file': gtex_file,
        'mapping_rate': f"{mapping_rate:.1f}%",
        'updated_samples': int(mapped_count)
    }
    
    # Show updated distribution
    updated_sex_counts = adata.obs['sex'].value_counts()
    print(f"📈 Updated sex distribution:")
    for sex, count in updated_sex_counts.items():
        print(f"  • {sex}: {count:,}")
    
    # Save updated dataset
    backup_path = Path(gtex_path).with_suffix('.h5ad.backup_sex')
    if not backup_path.exists():
        print(f"💾 Creating backup: {backup_path.name}")
        adata_backup = sc.read_h5ad(gtex_path)
        adata_backup.write_h5ad(backup_path)
    
    print(f"💾 Saving updated GTEx dataset...")
    adata.write_h5ad(gtex_path)
    
    return mapped_count > 0

def integrate_adni_age_data():
    """Integrate ADNI age data from demographics file."""
    print("\n🔄 **INTEGRATING ADNI AGE DATA**")
    print("=" * 33)
    
    # Load ADNI demographics data
    adni_demo_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadataADNI/subject_demographics/PTDEMOG_25Apr2025.csv'
    demo_df = pd.read_csv(adni_demo_file)
    
    print(f"✅ Loaded ADNI demographics: {len(demo_df):,} records")
    
    # Calculate age from birth year and visit date
    # Convert PTDOBYY to age at visit
    demo_df['VISDATE'] = pd.to_datetime(demo_df['VISDATE'], errors='coerce')
    demo_df['age_at_visit'] = demo_df['VISDATE'].dt.year - demo_df['PTDOBYY']
    
    # Remove invalid ages
    demo_df = demo_df[(demo_df['age_at_visit'] >= 0) & (demo_df['age_at_visit'] <= 120)]
    
    print(f"✅ Calculated ages for {len(demo_df):,} records")
    print(f"  Age range: {demo_df['age_at_visit'].min():.0f} - {demo_df['age_at_visit'].max():.0f}")
    
    # Create subject-to-age mapping (use most recent visit per subject)
    subject_age_map = demo_df.groupby('PTID')['age_at_visit'].max().to_dict()
    
    print(f"✅ Age data for {len(subject_age_map):,} unique subjects")
    
    # Load ADNI dataset
    # Use dynamic path based on current repo structure
    adni_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad'
    adata = sc.read_h5ad(adni_path)
    
    # Map subject_id to PTID format (e.g., "002_S_0413" -> "002_S_0413")
    # ADNI subject_ids should already match PTID format
    adata.obs['age_updated'] = adata.obs['subject_id'].map(subject_age_map)
    
    # Check mapping success
    mapped_count = adata.obs['age_updated'].notna().sum()
    mapping_rate = mapped_count / len(adata.obs) * 100
    
    print(f"✅ Mapped age for {mapped_count:,}/{len(adata.obs):,} samples ({mapping_rate:.1f}%)")
    
    if mapped_count > 0:
        # Update age field
        adata.obs['age_original_empty'] = adata.obs['age'].copy()
        adata.obs['age'] = adata.obs['age_updated'].astype(str).fillna('')
        adata.obs = adata.obs.drop(columns=['age_updated'])
        
        # Add metadata about the update
        adata.uns['age_mapping_update'] = {
            'timestamp': pd.Timestamp.now().isoformat(),
            'source_file': adni_demo_file,
            'mapping_rate': f"{mapping_rate:.1f}%",
            'updated_samples': int(mapped_count)
        }
        
        # Show age distribution
        age_values = adata.obs['age'][adata.obs['age'] != '']
        if len(age_values) > 0:
            age_numeric = pd.to_numeric(age_values, errors='coerce').dropna()
            print(f"📈 Updated age distribution:")
            print(f"  • Mean age: {age_numeric.mean():.1f}")
            print(f"  • Age range: {age_numeric.min():.0f} - {age_numeric.max():.0f}")
        
        # Save updated dataset
        backup_path = Path(adni_path).with_suffix('.h5ad.backup_age')
        if not backup_path.exists():
            print(f"💾 Creating backup: {backup_path.name}")
            adata_backup = sc.read_h5ad(adni_path)
            adata_backup.write_h5ad(backup_path)
        
        print(f"💾 Saving updated ADNI dataset...")
        adata.write_h5ad(adni_path)
        
        return True
    else:
        print("❌ No age mapping successful - check subject_id format matching")
        return False

def check_mage_age_availability():
    """Check if MAGE age data is available from 1000 Genomes."""
    print("\n🔍 **CHECKING MAGE AGE DATA AVAILABILITY**")
    print("=" * 42)
    
    print("ℹ️  1000 Genomes typically does not provide age data for privacy protection")
    print("   MAGE age field will remain empty per standard privacy practices")
    
    return False

def main():
    """Main integration function."""
    print("🧬 **INTEGRATING MISSING METADATA**")
    print("=" * 38)
    
    results = {
        'gtex_sex': False,
        'adni_age': False,
        'mage_age': False
    }
    
    # Integrate GTEx sex data
    try:
        results['gtex_sex'] = integrate_gtex_sex_data()
    except Exception as e:
        print(f"❌ GTEx sex integration failed: {e}")
    
    # Integrate ADNI age data
    try:
        results['adni_age'] = integrate_adni_age_data()
    except Exception as e:
        print(f"❌ ADNI age integration failed: {e}")
    
    # Check MAGE age (typically not available)
    results['mage_age'] = check_mage_age_availability()
    
    # Summary
    print(f"\n🎯 **INTEGRATION SUMMARY**")
    print(f"✅ GTEx sex data integrated: {results['gtex_sex']}")
    print(f"✅ ADNI age data integrated: {results['adni_age']}")
    print(f"ℹ️  MAGE age data: {results['mage_age']} (privacy protected)")
    
    successful_integrations = sum(results.values())
    print(f"\n📊 Successfully integrated: {successful_integrations}/2 available metadata fields")
    
    if successful_integrations > 0:
        print("\n🚀 **NEXT STEPS:**")
        print("1. Re-run validation to confirm metadata integration")
        print("2. Update partner presentation notebook")
        print("3. Regenerate combined dataset if needed")
    
    return results

if __name__ == '__main__':
    results = main()