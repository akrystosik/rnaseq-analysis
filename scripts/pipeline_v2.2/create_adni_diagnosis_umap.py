#!/usr/bin/env python3
"""
Create UMAP visualization for ADNI dataset colored by worst Alzheimer's diagnosis
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_diagnosis_data():
    """Load ADNI diagnosis data"""
    diagnosis_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/adni_diagnosis_export.csv"
    if Path(diagnosis_path).exists():
        print(f"Loading diagnosis data from: {diagnosis_path}")
        diagnosis_df = pd.read_csv(diagnosis_path)
        print(f"Diagnosis data shape: {diagnosis_df.shape}")
        print(f"Diagnosis distribution:")
        print(diagnosis_df['most_severe_diagnosis'].value_counts())
        return diagnosis_df
    else:
        print("Diagnosis file not found")
        return None

def load_adni_dataset(file_path, diagnosis_df=None):
    """Load ADNI dataset and add diagnosis information"""
    print(f"Loading ADNI dataset from: {file_path}")
    adata = sc.read_h5ad(file_path)
    print(f"ADNI dataset shape: {adata.shape}")
    
    print("Sample metadata columns:")
    for col in adata.obs.columns:
        print(f"  - {col}: {adata.obs[col].dtype}")
        if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category':
            unique_vals = adata.obs[col].nunique()
            print(f"    Unique values: {unique_vals}")
            if unique_vals <= 20:
                print(f"    Values: {list(adata.obs[col].unique())}")
    
    # Add diagnosis information if available
    if diagnosis_df is not None:
        print("Adding diagnosis information...")
        
        # Check if we have subject_id or donor_id
        if 'subject_id' in adata.obs.columns:
            id_col = 'subject_id'
        elif 'donor_id' in adata.obs.columns:
            id_col = 'donor_id'
        else:
            print("No suitable ID column found for diagnosis matching")
            return adata
        
        # Create diagnosis mapping
        diagnosis_map = diagnosis_df.set_index('donor_id')['most_severe_diagnosis'].to_dict()
        diagnosis_code_map = diagnosis_df.set_index('donor_id')['most_severe_diagnosis_code'].to_dict()
        
        # Map diagnosis to samples
        adata.obs['diagnosis'] = adata.obs[id_col].map(diagnosis_map)
        adata.obs['diagnosis_code'] = adata.obs[id_col].map(diagnosis_code_map)
        
        # Handle missing values
        missing_diagnosis = adata.obs['diagnosis'].isna().sum()
        print(f"Samples with missing diagnosis: {missing_diagnosis}")
        adata.obs['diagnosis'] = adata.obs['diagnosis'].fillna('Unknown')
        adata.obs['diagnosis_code'] = adata.obs['diagnosis_code'].fillna(0)
        
        print("Final diagnosis distribution in dataset:")
        print(adata.obs['diagnosis'].value_counts())
        
    return adata

def preprocess_for_umap(adata, n_top_genes=2000):
    """Preprocess ADNI data for UMAP generation"""
    print("Preprocessing ADNI data for UMAP...")
    
    # Make a copy to avoid modifying original
    adata_proc = adata.copy()
    
    # Basic preprocessing
    sc.pp.normalize_total(adata_proc, target_sum=1e4)
    sc.pp.log1p(adata_proc)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata_proc, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_top_genes)
    adata_proc.raw = adata_proc
    adata_proc = adata_proc[:, adata_proc.var.highly_variable]
    
    # Scale data
    sc.pp.scale(adata_proc, max_value=10)
    
    return adata_proc

def create_adni_diagnosis_umap(adata, output_dir):
    """Create UMAP for ADNI data colored by diagnosis"""
    print("Computing PCA...")
    sc.tl.pca(adata, svd_solver='arpack')
    
    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    print("Computing UMAP embedding...")
    sc.tl.umap(adata)
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Define diagnosis colors
    diagnosis_colors = {
        'Cognitively Normal': '#2ca02c',        # Green
        'Mild Cognitive Impairment': '#ff7f0e', # Orange  
        'Alzheimer\'s Disease': '#d62728',      # Red
        'Unknown': '#8c564b'                    # Brown
    }
    
    # Create diagnosis UMAP
    if 'diagnosis' in adata.obs.columns:
        print("Creating UMAP colored by diagnosis...")
        plt.figure(figsize=(12, 10))
        
        # Get unique diagnoses and create color palette
        unique_diagnoses = adata.obs['diagnosis'].unique()
        palette = [diagnosis_colors.get(diag, '#cccccc') for diag in unique_diagnoses]
        
        sc.pl.umap(adata, color='diagnosis', show=False, frameon=False,
                   legend_loc='right margin', size=50, alpha=0.8, palette=palette)
        plt.title('ADNI UMAP colored by Worst Alzheimer\'s Diagnosis', 
                 fontsize=18, pad=20, fontweight='bold')
        
        diagnosis_path = output_dir / 'adni_umap_by_diagnosis.png'
        plt.savefig(diagnosis_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"✓ Saved diagnosis UMAP to: {diagnosis_path}")
        plt.close()
    
    # Also create a version colored by sex for comparison
    if 'sex' in adata.obs.columns:
        print("Creating UMAP colored by sex...")
        plt.figure(figsize=(12, 10))
        
        sex_colors = {'male': '#1f77b4', 'female': '#ff7f0e', 'unknown': '#8c564b'}
        unique_sex = adata.obs['sex'].unique()
        palette = [sex_colors.get(s, '#cccccc') for s in unique_sex]
        
        sc.pl.umap(adata, color='sex', show=False, frameon=False,
                   legend_loc='right margin', size=50, alpha=0.8, palette=palette)
        plt.title('ADNI UMAP colored by Sex', fontsize=18, pad=20, fontweight='bold')
        
        sex_path = output_dir / 'adni_umap_by_sex.png'
        plt.savefig(sex_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"✓ Saved sex UMAP to: {sex_path}")
        plt.close()
    
    # Create combined plot if we have both diagnosis and sex
    if 'diagnosis' in adata.obs.columns and 'sex' in adata.obs.columns:
        print("Creating combined ADNI UMAP...")
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
        
        # Diagnosis plot
        unique_diagnoses = adata.obs['diagnosis'].unique()
        palette = [diagnosis_colors.get(diag, '#cccccc') for diag in unique_diagnoses]
        sc.pl.umap(adata, color='diagnosis', ax=axes[0], show=False, frameon=False,
                   size=40, alpha=0.8, palette=palette)
        axes[0].set_title('Alzheimer\'s Diagnosis', fontsize=16, fontweight='bold')
        
        # Sex plot
        unique_sex = adata.obs['sex'].unique()
        palette = [sex_colors.get(s, '#cccccc') for s in unique_sex]
        sc.pl.umap(adata, color='sex', ax=axes[1], show=False, frameon=False,
                   size=40, alpha=0.8, palette=palette)
        axes[1].set_title('Sex', fontsize=16, fontweight='bold')
        
        plt.suptitle('ADNI Dataset UMAP Analysis', fontsize=20, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        combined_path = output_dir / 'adni_combined_umap.png'
        plt.savefig(combined_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"✓ Saved combined ADNI UMAP to: {combined_path}")
        plt.close()
    
    # Create diagnosis summary statistics
    if 'diagnosis' in adata.obs.columns:
        print("\n" + "="*50)
        print("ADNI DIAGNOSIS SUMMARY")
        print("="*50)
        print(f"Total ADNI samples: {adata.n_obs}")
        print(f"Diagnosis distribution:")
        diag_counts = adata.obs['diagnosis'].value_counts()
        for diag, count in diag_counts.items():
            percentage = (count / adata.n_obs) * 100
            print(f"  {diag}: {count} ({percentage:.1f}%)")
        
        if 'sex' in adata.obs.columns:
            print(f"\nSex distribution:")
            sex_counts = adata.obs['sex'].value_counts()
            for sex, count in sex_counts.items():
                percentage = (count / adata.n_obs) * 100
                print(f"  {sex}: {count} ({percentage:.1f}%)")
        
        # Cross-tabulation
        if 'sex' in adata.obs.columns:
            print(f"\nDiagnosis by sex crosstab:")
            crosstab = pd.crosstab(adata.obs['diagnosis'], adata.obs['sex'], margins=True)
            print(crosstab)
        
        print("="*50)
    
    return adata

def main():
    # Load diagnosis data
    diagnosis_df = load_diagnosis_data()
    
    # Use the correct ADNI dataset file (latest_v2.2 has proper sex data)
    adni_paths = [
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/latest_v2.2/adni_standardized_v2.h5ad",
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/adni_standardized_preprocessed.h5ad",
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/adni_standardized_v2.h5ad"
    ]
    
    for adni_path in adni_paths:
        if Path(adni_path).exists():
            print(f"\n{'='*60}")
            print(f"Processing ADNI dataset: {adni_path}")
            print(f"{'='*60}")
            
            try:
                # Load ADNI dataset
                adata = load_adni_dataset(adni_path, diagnosis_df)
                
                # Preprocess for UMAP
                adata_proc = preprocess_for_umap(adata)
                
                # Create output directory
                output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/adni_diagnosis_umap_results"
                
                # Create ADNI diagnosis UMAP
                adata_final = create_adni_diagnosis_umap(adata_proc, output_dir)
                
                print(f"\n✓ Successfully created ADNI diagnosis UMAPs")
                print(f"  Output directory: {output_dir}")
                
                # Only process the first successful dataset
                break
                
            except Exception as e:
                print(f"Error processing {adni_path}: {str(e)}")
                import traceback
                traceback.print_exc()
                continue
    
    print("\n" + "="*60)
    print("ADNI DIAGNOSIS UMAP GENERATION COMPLETE!")
    print("="*60)

if __name__ == "__main__":
    main()