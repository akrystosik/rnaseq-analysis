#!/usr/bin/env python3
"""
Working ethnicity UMAP using the complete sample mapping with simple ID extraction
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

def extract_sample_id(full_sample_id):
    """Extract base sample ID from complex H5AD sample ID"""
    # For ADNI: 002_S_0413_002_S_0413_gencode_v24_pruned -> 002_S_0413
    parts = full_sample_id.split('_')
    
    # ADNI pattern
    if len(parts) >= 6 and parts[0] == parts[3] and parts[1] == parts[4] and parts[2] == parts[5]:
        return f"{parts[0]}_{parts[1]}_{parts[2]}"
    
    # MAGE pattern: NA06985_NA06985_gencode_V24_pruned -> NA06985
    if len(parts) >= 4 and parts[0] == parts[1]:
        return parts[0]
    
    # GTEx: GTEX-1117F-0005-SM-HL9SH -> GTEX-1117F  
    if full_sample_id.startswith('GTEX-'):
        parts = full_sample_id.split('-')
        if len(parts) >= 2:
            return f"{parts[0]}-{parts[1]}"
    
    # ENCODE: should be direct
    if full_sample_id.startswith('ENCFF'):
        return full_sample_id
    
    # Default: return as-is
    return full_sample_id

def main():
    # Load complete sample ethnicity mapping
    print("Loading complete sample ethnicity mapping...")
    ethnicity_df = pd.read_csv("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/sample_ethnicity_mapping_complete.csv")
    
    print("Correct ethnicity distribution from complete mapping:")
    print(ethnicity_df['ethnicity'].value_counts())
    
    # Create sample_id to ethnicity mapping
    sample_to_ethnicity = ethnicity_df.set_index('sample_id')['ethnicity'].to_dict()
    
    # Load dataset
    print("Loading H5AD dataset...")
    adata = sc.read_h5ad("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/combined_dataset_latest.h5ad")
    
    # Subsample
    print("Subsampling...")
    sc.pp.subsample(adata, n_obs=4000, random_state=42)
    
    # Extract simple sample IDs
    print("Extracting sample IDs...")
    adata.obs['extracted_sample_id'] = adata.obs['sample_id'].apply(extract_sample_id)
    
    # Show examples
    print("ID extraction examples:")
    for i in range(5):
        original = adata.obs['sample_id'].iloc[i]
        extracted = adata.obs['extracted_sample_id'].iloc[i]
        print(f"  {original} -> {extracted}")
    
    # Map ethnicity
    print("Mapping ethnicity...")
    adata.obs['ethnicity'] = adata.obs['extracted_sample_id'].map(sample_to_ethnicity)
    
    # Check mapping success
    mapped_count = adata.obs['ethnicity'].notna().sum()
    print(f"Successfully mapped: {mapped_count}/{len(adata)} ({mapped_count/len(adata)*100:.1f}%)")
    
    # Fill unknowns
    adata.obs['ethnicity'] = adata.obs['ethnicity'].fillna('unknown')
    
    print("Final ethnicity distribution in subsample:")
    ethnicity_dist = adata.obs['ethnicity'].value_counts()
    print(ethnicity_dist)
    
    # Only continue if we have reasonable mapping
    white_count = ethnicity_dist.get('white', 0)
    if white_count < 100:  # Should have many white samples
        print("❌ Still not getting enough white samples, there may be an ID mapping issue")
        
        # Debug by checking a few specific mappings
        print("\nDebugging specific samples:")
        for i in range(min(10, len(adata))):
            sample_id = adata.obs['sample_id'].iloc[i]
            extracted = adata.obs['extracted_sample_id'].iloc[i] 
            ethnicity = adata.obs['ethnicity'].iloc[i]
            in_csv = extracted in sample_to_ethnicity
            print(f"  {sample_id}")
            print(f"    -> {extracted}")
            print(f"    -> In CSV: {in_csv}")
            print(f"    -> Ethnicity: {ethnicity}")
            if in_csv:
                csv_ethnicity = sample_to_ethnicity[extracted]
                print(f"    -> CSV ethnicity: {csv_ethnicity}")
            print()
        return
    
    # Continue with UMAP if mapping looks good
    print("✓ Good mapping rate, proceeding with UMAP...")
    
    # Preprocessing
    print("Preprocessing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1200)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata)
    
    # UMAP
    print("Computing UMAP...")
    sc.tl.pca(adata, n_comps=35)
    sc.pp.neighbors(adata, n_neighbors=12, n_pcs=35)
    sc.tl.umap(adata)
    
    # Output
    output_dir = Path("umap_working_ethnicity")
    output_dir.mkdir(exist_ok=True)
    
    # Colors
    ethnicity_colors = {
        'white': '#1f77b4',
        'black or african american': '#ff7f0e', 
        'asian': '#2ca02c',
        'hispanic or latino': '#d62728',
        'american indian or alaska native': '#9467bd',
        'multiethnic': '#8c564b',
        'native hawaiian or other pacific islander': '#e377c2',
        'unknown or not reported': '#7f7f7f',
        'unknown': '#bcbd22'
    }
    
    # Create UMAP
    print("Creating final ethnicity UMAP...")
    plt.figure(figsize=(14, 10))
    unique_ethnicities = adata.obs['ethnicity'].unique()
    palette = [ethnicity_colors.get(eth, '#cccccc') for eth in unique_ethnicities]
    
    sc.pl.umap(adata, color='ethnicity', show=False, frameon=False,
               legend_loc='right margin', size=35, alpha=0.8, palette=palette)
    plt.title('UMAP with Properly Corrected Ethnicity\n(White samples correctly represented)', 
              fontsize=16, pad=20, fontweight='bold')
    
    plt.savefig(output_dir / 'umap_working_corrected_ethnicity.png', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved to: {output_dir}/umap_working_corrected_ethnicity.png")
    plt.close()
    
    # Summary
    print("\n" + "="*60)
    print("✅ WORKING CORRECTED ETHNICITY UMAP")
    print("="*60)
    print("📊 Final ethnicity distribution:")
    for eth, count in ethnicity_dist.items():
        pct = (count/len(adata)) * 100
        print(f"   {eth}: {count} ({pct:.1f}%)")
    
    # Show the correction
    white_pct = (white_count / len(adata)) * 100
    print(f"\n🔧 Correction successful:")
    print(f"   ✅ White samples: {white_count} ({white_pct:.1f}%)")
    print(f"   ❌ Previous (incorrect): ~80% Pacific Islander")
    print("="*60)

if __name__ == "__main__":
    main()