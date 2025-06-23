#!/usr/bin/env python3
"""
Create GTEx UMAP with all 54 tissues visible in legend for diversity presentation
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import matplotlib.patches as mpatches

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=150, facecolor='white')

def create_full_tissue_gtex_umap():
    """Create GTEx UMAP with all 54 tissues in legend"""
    
    # Load the processed data
    gtex_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest_v2.2/gtex_standardized_preprocessed.h5ad"
    print(f"Loading GTEx dataset from: {gtex_path}")
    adata = sc.read_h5ad(gtex_path)
    
    # Preprocess for UMAP
    print("Preprocessing data...")
    adata_proc = adata.copy()
    sc.pp.normalize_total(adata_proc, target_sum=1e4)
    sc.pp.log1p(adata_proc)
    sc.pp.highly_variable_genes(adata_proc, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
    adata_proc.raw = adata_proc
    adata_proc = adata_proc[:, adata_proc.var.highly_variable]
    sc.pp.scale(adata_proc, max_value=10)
    
    # Compute UMAP
    print("Computing UMAP...")
    sc.tl.pca(adata_proc, svd_solver='arpack')
    sc.pp.neighbors(adata_proc, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_proc)
    
    # Use the detailed tissue column (54 tissues)
    tissue_col = 'tissue'  # This has 54 unique values
    tissues = sorted(adata_proc.obs[tissue_col].unique())
    
    print(f"Found {len(tissues)} unique tissues")
    
    # Generate distinct colors for all tissues
    # Use a combination of color palettes to get 54 distinct colors
    colors1 = sns.color_palette("husl", 20)
    colors2 = sns.color_palette("Set3", 12)
    colors3 = sns.color_palette("tab20", 20)
    colors4 = sns.color_palette("Dark2", 8)
    
    all_colors = colors1 + colors2 + colors3 + colors4
    # Take only what we need
    tissue_colors = {tissue: all_colors[i % len(all_colors)] for i, tissue in enumerate(tissues)}
    
    # Create a wide figure to accommodate the legend
    print("Creating full tissue diversity figure...")
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(20, 12))  # Very wide for legend
    
    # Plot each tissue
    for tissue in tissues:
        mask = adata_proc.obs[tissue_col] == tissue
        color = tissue_colors[tissue]
        
        ax.scatter(adata_proc.obsm['X_umap'][mask, 0], 
                  adata_proc.obsm['X_umap'][mask, 1],
                  c=[color], s=8, alpha=0.8, label=tissue, rasterized=True)
    
    # Clean up the plot
    ax.set_xlabel('UMAP 1', fontsize=16, fontweight='bold')
    ax.set_ylabel('UMAP 2', fontsize=16, fontweight='bold')
    ax.set_title('GTEx Human Tissue Diversity: 54 Tissue Types', fontsize=20, fontweight='bold', pad=20)
    
    # Remove axes ticks and spines for cleaner look
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    # Create legend with all tissues - organize in columns
    # Clean up tissue names for better readability
    clean_tissue_names = []
    for tissue in tissues:
        # Capitalize each word and clean up formatting
        clean_name = tissue.replace(' - ', ' • ').title()
        clean_name = clean_name.replace('Ebv-Transformed', 'EBV-transformed')
        clean_name = clean_name.replace('Ba9', 'BA9').replace('Ba24', 'BA24')
        clean_name = clean_name.replace('(Ba9)', '(BA9)').replace('(Ba24)', '(BA24)')
        clean_name = clean_name.replace('C-1', 'C1')
        clean_tissue_names.append(clean_name)
    
    # Create legend elements
    legend_elements = [mpatches.Patch(color=tissue_colors[tissue], label=clean_name) 
                      for tissue, clean_name in zip(tissues, clean_tissue_names)]
    
    # Position legend in multiple columns to the right
    legend = ax.legend(handles=legend_elements, 
                      loc='center left', 
                      bbox_to_anchor=(1.02, 0.5),
                      fontsize=9, 
                      frameon=True, 
                      fancybox=True, 
                      shadow=True,
                      ncol=2,  # Two columns to fit better
                      columnspacing=1.5,
                      handletextpad=0.5)
    
    # Add sample count annotation
    ax.text(0.03, 0.97, f'{adata_proc.n_obs:,} samples\n{len(tissues)} tissue types\nShowing tissue diversity', 
            transform=ax.transAxes, fontsize=14, verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.7", facecolor='white', alpha=0.9, edgecolor='black'))
    
    plt.tight_layout()
    
    # Save high-quality figure
    output_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/gtex_tissue_umap_results")
    output_path = output_dir / 'gtex_all_54_tissues_presentation.png'
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', 
                edgecolor='none', format='png')
    print(f"✓ Saved full tissue diversity UMAP to: {output_path}")
    
    # Also save as PDF
    pdf_path = output_dir / 'gtex_all_54_tissues_presentation.pdf'
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white', 
                edgecolor='none', format='pdf')
    print(f"✓ Saved vector PDF to: {pdf_path}")
    
    plt.close()
    
    # Create an alternative layout with legend below the plot
    print("Creating alternative layout with legend below...")
    fig, ax = plt.subplots(figsize=(16, 14))  # Taller for bottom legend
    
    # Plot each tissue again
    for tissue in tissues:
        mask = adata_proc.obs[tissue_col] == tissue
        color = tissue_colors[tissue]
        
        ax.scatter(adata_proc.obsm['X_umap'][mask, 0], 
                  adata_proc.obsm['X_umap'][mask, 1],
                  c=[color], s=10, alpha=0.8, label=tissue, rasterized=True)
    
    # Clean up the plot
    ax.set_xlabel('UMAP 1', fontsize=16, fontweight='bold')
    ax.set_ylabel('UMAP 2', fontsize=16, fontweight='bold')
    ax.set_title('GTEx Human Tissue Diversity: 54 Tissue Types', fontsize=20, fontweight='bold', pad=20)
    
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    # Legend below the plot in multiple columns
    legend_elements = [mpatches.Patch(color=tissue_colors[tissue], label=clean_name) 
                      for tissue, clean_name in zip(tissues, clean_tissue_names)]
    
    # Position legend below plot with many columns
    legend = ax.legend(handles=legend_elements, 
                      loc='upper center', 
                      bbox_to_anchor=(0.5, -0.05),
                      fontsize=8, 
                      frameon=True, 
                      fancybox=True, 
                      shadow=True,
                      ncol=6,  # Six columns to fit horizontally
                      columnspacing=1.0,
                      handletextpad=0.3)
    
    # Add sample count annotation
    ax.text(0.03, 0.97, f'{adata_proc.n_obs:,} samples\n{len(tissues)} tissue types', 
            transform=ax.transAxes, fontsize=14, verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.7", facecolor='white', alpha=0.9, edgecolor='black'))
    
    plt.tight_layout()
    
    # Save alternative layout
    alt_path = output_dir / 'gtex_all_54_tissues_bottom_legend.png'
    plt.savefig(alt_path, dpi=300, bbox_inches='tight', facecolor='white', 
                edgecolor='none', format='png')
    print(f"✓ Saved bottom legend version to: {alt_path}")
    
    alt_pdf = output_dir / 'gtex_all_54_tissues_bottom_legend.pdf'
    plt.savefig(alt_pdf, bbox_inches='tight', facecolor='white', 
                edgecolor='none', format='pdf')
    print(f"✓ Saved bottom legend PDF to: {alt_pdf}")
    
    plt.close()
    
    # Print tissue list for reference
    print("\nAll 54 GTEx tissues included:")
    print("="*60)
    for i, (tissue, clean_name) in enumerate(zip(tissues, clean_tissue_names), 1):
        print(f"{i:2d}. {clean_name}")
    
    print("\n" + "="*60)
    print("FULL TISSUE DIVERSITY UMAP GENERATION COMPLETE!")
    print("="*60)
    print("Files created:")
    print(f"  - gtex_all_54_tissues_presentation.png (side legend)")
    print(f"  - gtex_all_54_tissues_presentation.pdf (side legend)")
    print(f"  - gtex_all_54_tissues_bottom_legend.png (bottom legend)")  
    print(f"  - gtex_all_54_tissues_bottom_legend.pdf (bottom legend)")
    print("="*60)

if __name__ == "__main__":
    create_full_tissue_gtex_umap()