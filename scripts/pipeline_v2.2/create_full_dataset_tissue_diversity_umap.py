#!/usr/bin/env python3
"""
Create comprehensive tissue diversity UMAP using the full dataset (21,000+ samples)
This expands beyond the 10,000 sample limit to capture maximum tissue representation
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import time

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set scanpy settings for full dataset processing
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_full_dataset():
    """Load the complete combined dataset"""
    print("Loading full combined dataset...")
    
    # Try different dataset paths
    dataset_paths = [
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/latest_v2.2/combined_dataset_all_genes_sparse.h5ad",
        "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/combined_dataset_latest.h5ad"
    ]
    
    for path in dataset_paths:
        if Path(path).exists():
            print(f"Loading dataset from: {path}")
            adata = sc.read_h5ad(path)
            print(f"✓ Dataset loaded: {adata.shape}")
            
            # Show dataset composition
            if 'dataset' in adata.obs.columns:
                print("Dataset composition:")
                dataset_counts = adata.obs['dataset'].value_counts()
                for dataset, count in dataset_counts.items():
                    percentage = (count / adata.n_obs) * 100
                    print(f"  {dataset}: {count} ({percentage:.1f}%)")
            
            # Show tissue diversity
            if 'tissue' in adata.obs.columns:
                tissue_counts = adata.obs['tissue'].value_counts()
                print(f"Total tissue types: {len(tissue_counts)}")
                print("Top 10 tissues by sample count:")
                for tissue, count in tissue_counts.head(10).items():
                    print(f"  {tissue}: {count}")
                    
            return adata
    
    raise FileNotFoundError("Could not find combined dataset file")

def preprocess_full_dataset(adata, n_top_genes=2000):
    """Preprocess full dataset for UMAP with memory-efficient approach"""
    print(f"Preprocessing full dataset ({adata.n_obs} samples, {adata.n_vars} genes)...")
    
    # Make a copy to avoid modifying original
    adata_proc = adata.copy()
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    print("  Step 1/4: TPM normalization...")
    sc.pp.normalize_total(adata_proc, target_sum=1e4)
    
    print("  Step 2/4: Log transformation...")
    sc.pp.log1p(adata_proc)
    
    print("  Step 3/4: Highly variable gene selection...")
    sc.pp.highly_variable_genes(adata_proc, min_mean=0.0125, max_mean=3, 
                               min_disp=0.5, n_top_genes=n_top_genes)
    
    # Store raw data and filter to highly variable genes
    adata_proc.raw = adata_proc
    adata_proc = adata_proc[:, adata_proc.var.highly_variable]
    print(f"  Reduced to {adata_proc.n_vars} highly variable genes")
    
    print("  Step 4/4: Z-score scaling...")
    sc.pp.scale(adata_proc, max_value=10)
    
    print("✓ Preprocessing complete")
    return adata_proc

def compute_full_dataset_umap(adata, n_pcs=40, n_neighbors=15):
    """Compute UMAP for full dataset with optimized parameters"""
    print("Computing UMAP for full dataset...")
    
    print(f"  Step 1/3: PCA ({n_pcs} components)...")
    start_time = time.time()
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)
    print(f"    PCA completed in {time.time() - start_time:.1f}s")
    
    print(f"  Step 2/3: Neighborhood graph (k={n_neighbors})...")
    start_time = time.time()
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    print(f"    Neighborhood graph completed in {time.time() - start_time:.1f}s")
    
    print("  Step 3/3: UMAP embedding...")
    start_time = time.time()
    sc.tl.umap(adata, min_dist=0.3, spread=1.0, random_state=42)
    print(f"    UMAP completed in {time.time() - start_time:.1f}s")
    
    print("✓ UMAP computation complete")
    return adata

def create_comprehensive_tissue_colors():
    """Create comprehensive color palette for all tissue types"""
    # Use a large qualitative palette
    colors = plt.cm.tab20.colors + plt.cm.tab20b.colors + plt.cm.tab20c.colors + plt.cm.Set3.colors
    
    # Convert to hex
    hex_colors = ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) for r, g, b in colors]
    
    return hex_colors

def create_full_dataset_tissue_umap(adata, output_dir):
    """Create comprehensive tissue diversity UMAP visualization"""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print("Creating full dataset tissue diversity UMAP...")
    
    # Get tissue information
    if 'tissue' not in adata.obs.columns:
        print("❌ No tissue column found in dataset")
        return
    
    # Remove any samples with missing tissue information
    valid_tissue_mask = adata.obs['tissue'].notna() & (adata.obs['tissue'] != '') & (adata.obs['tissue'] != 'nan')
    adata_clean = adata[valid_tissue_mask].copy()
    
    print(f"Samples with valid tissue information: {adata_clean.n_obs}/{adata.n_obs}")
    
    # Get tissue statistics
    tissue_counts = adata_clean.obs['tissue'].value_counts()
    print(f"Total tissue types: {len(tissue_counts)}")
    
    # Create comprehensive color palette
    unique_tissues = tissue_counts.index.tolist()
    tissue_colors = create_comprehensive_tissue_colors()
    
    # Ensure we have enough colors
    if len(unique_tissues) > len(tissue_colors):
        # Generate additional colors using a color map
        additional_colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(unique_tissues) - len(tissue_colors)))
        additional_hex = ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) for r, g, b in additional_colors[:, :3]]
        tissue_colors.extend(additional_hex)
    
    # Create tissue color mapping
    tissue_color_map = {tissue: tissue_colors[i] for i, tissue in enumerate(unique_tissues)}
    
    # Create main tissue diversity plot
    fig, ax = plt.subplots(figsize=(20, 14))
    
    # Plot points by tissue
    for tissue in unique_tissues:
        tissue_mask = adata_clean.obs['tissue'] == tissue
        tissue_data = adata_clean[tissue_mask]
        
        if tissue_data.n_obs > 0:
            x = tissue_data.obsm['X_umap'][:, 0]
            y = tissue_data.obsm['X_umap'][:, 1]
            
            ax.scatter(x, y, 
                      c=tissue_color_map[tissue], 
                      s=8, 
                      alpha=0.7,
                      label=f"{tissue} ({tissue_data.n_obs})",
                      rasterized=True)  # Rasterize for better performance with many points
    
    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.set_title(f'Tissue Diversity: Full Multi-Dataset RNA-seq Collection\n{adata_clean.n_obs:,} samples across {len(unique_tissues)} tissue types', 
                fontsize=18, pad=20, fontweight='bold')
    
    # Remove axis ticks for cleaner look
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Add dataset composition info
    if 'dataset' in adata_clean.obs.columns:
        dataset_info = []
        dataset_counts = adata_clean.obs['dataset'].value_counts()
        for dataset, count in dataset_counts.items():
            percentage = (count / adata_clean.n_obs) * 100
            dataset_info.append(f"{dataset}: {count:,} ({percentage:.1f}%)")
        
        info_text = f"Total samples: {adata_clean.n_obs:,}\nTissue types: {len(unique_tissues)}\n\nDataset composition:\n" + "\n".join(dataset_info)
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Create legend in separate subplot for better organization
    fig2, ax2 = plt.subplots(figsize=(8, max(12, len(unique_tissues) * 0.3)))
    ax2.axis('off')
    
    # Sort tissues by sample count for legend
    sorted_tissues = tissue_counts.index.tolist()
    
    # Create legend entries
    legend_elements = []
    for tissue in sorted_tissues:
        count = tissue_counts[tissue]
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=tissue_color_map[tissue], 
                                        markersize=8, 
                                        label=f"{tissue} ({count})"))
    
    # Split legend into columns if too many tissues
    ncol = max(1, len(legend_elements) // 40)  # Max 40 entries per column
    ax2.legend(handles=legend_elements, loc='center', fontsize=9, ncol=ncol)
    ax2.set_title('Tissue Types with Sample Counts', fontsize=14, fontweight='bold')
    
    # Save main plot
    main_path = output_dir / 'full_dataset_tissue_diversity_umap.png'
    fig.savefig(main_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved main tissue UMAP to: {main_path}")
    plt.close(fig)
    
    # Save legend
    legend_path = output_dir / 'full_dataset_tissue_diversity_legend.png'
    fig2.savefig(legend_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved tissue legend to: {legend_path}")
    plt.close(fig2)
    
    # Create PDF versions for presentations
    pdf_main_path = output_dir / 'full_dataset_tissue_diversity_umap.pdf'
    pdf_legend_path = output_dir / 'full_dataset_tissue_diversity_legend.pdf'
    
    # Recreate for PDF
    fig, ax = plt.subplots(figsize=(20, 14))
    for tissue in unique_tissues:
        tissue_mask = adata_clean.obs['tissue'] == tissue
        tissue_data = adata_clean[tissue_mask]
        
        if tissue_data.n_obs > 0:
            x = tissue_data.obsm['X_umap'][:, 0]
            y = tissue_data.obsm['X_umap'][:, 1]
            
            ax.scatter(x, y, 
                      c=tissue_color_map[tissue], 
                      s=8, 
                      alpha=0.7,
                      label=f"{tissue} ({tissue_data.n_obs})")
    
    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.set_title(f'Tissue Diversity: Full Multi-Dataset RNA-seq Collection\n{adata_clean.n_obs:,} samples across {len(unique_tissues)} tissue types', 
                fontsize=18, pad=20, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])
    
    if 'dataset' in adata_clean.obs.columns:
        dataset_info = []
        dataset_counts = adata_clean.obs['dataset'].value_counts()
        for dataset, count in dataset_counts.items():
            percentage = (count / adata_clean.n_obs) * 100
            dataset_info.append(f"{dataset}: {count:,} ({percentage:.1f}%)")
        
        info_text = f"Total samples: {adata_clean.n_obs:,}\nTissue types: {len(unique_tissues)}\n\nDataset composition:\n" + "\n".join(dataset_info)
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    fig.savefig(pdf_main_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved PDF tissue UMAP to: {pdf_main_path}")
    plt.close(fig)
    
    # Summary statistics
    print("\n" + "="*70)
    print("FULL DATASET TISSUE DIVERSITY ANALYSIS SUMMARY")
    print("="*70)
    print(f"📊 Total samples analyzed: {adata_clean.n_obs:,}")
    print(f"🧬 Total genes used: {adata.n_vars:,}")
    print(f"🏥 Total tissue types: {len(unique_tissues)}")
    
    if 'dataset' in adata_clean.obs.columns:
        print(f"\n📋 Dataset composition:")
        for dataset, count in dataset_counts.items():
            percentage = (count / adata_clean.n_obs) * 100
            print(f"   {dataset}: {count:,} samples ({percentage:.1f}%)")
    
    print(f"\n🔬 Top 10 tissues by sample count:")
    for tissue, count in tissue_counts.head(10).items():
        percentage = (count / adata_clean.n_obs) * 100
        print(f"   {tissue}: {count} ({percentage:.1f}%)")
    
    print("="*70)
    
    return adata_clean

def main():
    # Load full dataset
    adata = load_full_dataset()
    
    print(f"\n{'='*60}")
    print(f"FULL DATASET TISSUE DIVERSITY UMAP ANALYSIS")
    print(f"{'='*60}")
    print(f"Processing {adata.n_obs:,} samples with {adata.n_vars:,} genes")
    
    # Preprocess
    adata_proc = preprocess_full_dataset(adata, n_top_genes=2000)
    
    # Compute UMAP
    adata_final = compute_full_dataset_umap(adata_proc, n_pcs=40, n_neighbors=15)
    
    # Create visualizations
    output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/full_dataset_tissue_analysis"
    adata_final = create_full_dataset_tissue_umap(adata_final, output_dir)
    
    print(f"\n✅ FULL DATASET ANALYSIS COMPLETE!")
    print(f"📁 Output directory: {output_dir}")
    print(f"📈 Processed {adata_final.n_obs:,} samples across {adata_final.obs['tissue'].nunique()} tissue types")
    print("="*60)

if __name__ == "__main__":
    main()