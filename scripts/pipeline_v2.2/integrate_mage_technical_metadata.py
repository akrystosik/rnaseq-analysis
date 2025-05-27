#!/usr/bin/env python3
"""
Integrate MAGE technical metadata from the comprehensive sample metadata file.
This adds valuable technical quality metrics like RIN scores, batch information, etc.
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

def integrate_mage_technical_metadata():
    """Integrate MAGE technical metadata including RIN scores and batch information."""
    print("🔄 **INTEGRATING MAGE TECHNICAL METADATA**")
    print("=" * 42)
    
    # Load MAGE technical metadata
    mage_metadata_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/sample.metadata.MAGE.v1.0.txt'
    metadata_df = pd.read_csv(mage_metadata_file, sep='\t')
    
    print(f"✅ Loaded MAGE technical metadata: {len(metadata_df):,} samples")
    print(f"   Columns: {list(metadata_df.columns)}")
    
    # Load MAGE dataset
    mage_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137/mage_standardized_preprocessed.h5ad'
    adata = sc.read_h5ad(mage_path)
    
    print(f"✅ Loaded MAGE dataset: {adata.shape}")
    
    # Create mapping using sample_kgpID (has 100% overlap)
    metadata_mapping = metadata_df.set_index('sample_kgpID')
    
    # Check mapping success
    adata_subjects = set(adata.obs['subject_id'])
    metadata_subjects = set(metadata_df['sample_kgpID'])
    overlap = len(adata_subjects & metadata_subjects)
    
    print(f"✅ Subject overlap: {overlap:,}/{len(adata_subjects):,} ({overlap/len(adata_subjects)*100:.1f}%)")
    
    if overlap > 0:
        # Map technical metadata
        print("\\n🔧 **Integrating Technical Metadata:**")
        
        # RIN Score (RNA Integrity Number)
        adata.obs['rna_integrity_number'] = adata.obs['subject_id'].map(metadata_mapping['RIN'])
        rin_mapped = adata.obs['rna_integrity_number'].notna().sum()
        print(f"   • RIN scores: {rin_mapped:,} samples mapped")
        if rin_mapped > 0:
            rin_stats = adata.obs['rna_integrity_number'].describe()
            print(f"     Range: {rin_stats['min']:.1f} - {rin_stats['max']:.1f}, Mean: {rin_stats['mean']:.1f}")
        
        # Batch information
        adata.obs['batch'] = adata.obs['subject_id'].map(metadata_mapping['batch'])
        batch_mapped = adata.obs['batch'].notna().sum()
        print(f"   • Batch info: {batch_mapped:,} samples mapped")
        if batch_mapped > 0:
            batch_counts = adata.obs['batch'].value_counts()
            print(f"     Batches: {len(batch_counts)} total, samples per batch: {batch_counts.min()}-{batch_counts.max()}")
        
        # Continental group (more detailed than our current ethnicity grouping)
        adata.obs['continental_group'] = adata.obs['subject_id'].map(metadata_mapping['continentalGroup'])
        continental_mapped = adata.obs['continental_group'].notna().sum()
        print(f"   • Continental groups: {continental_mapped:,} samples mapped")
        if continental_mapped > 0:
            continental_counts = adata.obs['continental_group'].value_counts()
            print(f"     Groups: {dict(continental_counts)}")
        
        # RNA concentration and quality metrics
        adata.obs['rna_concentration_ng_ul'] = adata.obs['subject_id'].map(metadata_mapping['RNAQubitConc_ng-ul'])
        adata.obs['rna_total_amount_ng'] = adata.obs['subject_id'].map(metadata_mapping['RNAQubitTotalAmount_ng'])
        conc_mapped = adata.obs['rna_concentration_ng_ul'].notna().sum()
        print(f"   • RNA quality metrics: {conc_mapped:,} samples mapped")
        
        # Number of reads (sequencing depth)
        adata.obs['num_reads'] = adata.obs['subject_id'].map(metadata_mapping['numReads'])
        reads_mapped = adata.obs['num_reads'].notna().sum()
        print(f"   • Sequencing depth: {reads_mapped:,} samples mapped")
        if reads_mapped > 0:
            reads_stats = adata.obs['num_reads'].describe()
            print(f"     Read count range: {reads_stats['min']:,.0f} - {reads_stats['max']:,.0f}")
            print(f"     Mean reads: {reads_stats['mean']:,.0f}")
        
        # Add metadata about the update
        adata.uns['mage_technical_metadata_update'] = {
            'timestamp': pd.Timestamp.now().isoformat(),
            'source_file': mage_metadata_file,
            'mapped_samples': int(overlap),
            'technical_fields_added': [
                'rna_integrity_number', 'batch', 'continental_group', 
                'rna_concentration_ng_ul', 'rna_total_amount_ng', 'num_reads'
            ]
        }
        
        # Save updated dataset
        backup_path = Path(mage_path).with_suffix('.h5ad.backup_technical')
        if not backup_path.exists():
            print(f"\\n💾 Creating backup: {backup_path.name}")
            adata_backup = sc.read_h5ad(mage_path)
            adata_backup.write_h5ad(backup_path)
        
        print(f"💾 Saving updated MAGE dataset...")
        adata.write_h5ad(mage_path)
        
        # Summary
        print(f"\\n🎯 **MAGE TECHNICAL METADATA INTEGRATION COMPLETE**")
        print(f"   ✅ RIN scores: Added for quality assessment")
        print(f"   ✅ Batch information: Added for technical variation analysis")
        print(f"   ✅ Continental groups: Enhanced population stratification")
        print(f"   ✅ RNA quality metrics: Concentration and total amounts")
        print(f"   ✅ Sequencing depth: Read count information")
        print(f"   🌟 Result: MAGE now has comprehensive technical metadata")
        
        return True
    else:
        print("❌ No overlap found - check subject ID mapping")
        return False

def main():
    """Main function."""
    success = integrate_mage_technical_metadata()
    
    if success:
        print("\\n🚀 **NEXT STEPS:**")
        print("1. Re-run validation to confirm technical metadata integration")
        print("2. Update partner presentation to highlight technical metadata")
        print("3. Document technical quality metrics for partner")
    
    return success

if __name__ == '__main__':
    main()