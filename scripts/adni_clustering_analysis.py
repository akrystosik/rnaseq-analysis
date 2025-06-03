#!/usr/bin/env python3
"""
ADNI Gene Expression Clustering Analysis

This script performs clustering analysis on ADNI gene expression data to:
1. Visualize sample clustering patterns using PCA and UMAP
2. Apply various clustering algorithms (K-means, hierarchical)
3. Analyze cluster separation and relationship to disease labels
4. Generate publication-quality visualizations
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score
import umap
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_adni_data(data_path=None):
    """Load ADNI preprocessed data."""
    if data_path is None:
        data_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest_v2.2/adni_standardized_preprocessed.h5ad'
    
    print(f"Loading ADNI data from: {data_path}")
    adata = sc.read_h5ad(data_path)
    
    print(f"✅ Loaded ADNI data: {adata.n_obs} samples × {adata.n_vars} genes")
    
    # Show available metadata columns
    print("\n📋 Available metadata columns:")
    for col in adata.obs.columns:
        print(f"  {col}")
    
    # Show diagnosis distribution
    if 'worst_diagnosis_label' in adata.obs.columns:
        print("\n📊 Diagnosis distribution:")
        diagnosis_counts = adata.obs['worst_diagnosis_label'].value_counts()
        for diagnosis, count in diagnosis_counts.items():
            print(f"  {diagnosis}: {count}")
    
    return adata

def preprocess_for_clustering(adata, n_top_genes=2000):
    """Preprocess data for clustering analysis."""
    print(f"\n🔧 Preprocessing data for clustering...")
    
    # Copy for analysis
    adata_work = adata.copy()
    
    # Filter genes by variability - use simple cell_ranger flavor instead
    sc.pp.highly_variable_genes(adata_work, n_top_genes=n_top_genes, flavor='cell_ranger')
    
    print(f"Selected {n_top_genes} highly variable genes for clustering")
    
    # Keep only highly variable genes
    adata_work = adata_work[:, adata_work.var.highly_variable]
    
    # Log transform if not already done
    if not hasattr(adata_work, 'raw'):
        sc.pp.log1p(adata_work)
    
    # Scale data
    sc.pp.scale(adata_work, max_value=10)
    
    return adata_work

def perform_dimensionality_reduction(adata):
    """Perform PCA and UMAP dimensionality reduction."""
    print("\n📐 Performing dimensionality reduction...")
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    
    # UMAP
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)
    
    print("✅ Completed PCA and UMAP")
    
    return adata

def perform_clustering(adata, methods=['kmeans', 'hierarchical']):
    """Perform various clustering methods."""
    print("\n🎯 Performing clustering analysis...")
    
    # Get PCA representation for clustering
    X_pca = adata.obsm['X_pca'][:, :30]  # Use first 30 PCs
    
    results = {}
    
    # K-means clustering (try different k values)
    if 'kmeans' in methods:
        print("  Running K-means clustering...")
        for k in [2, 3, 4, 5]:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            clusters = kmeans.fit_predict(X_pca)
            silhouette = silhouette_score(X_pca, clusters)
            
            adata.obs[f'kmeans_k{k}'] = clusters.astype(str)
            results[f'kmeans_k{k}'] = {
                'clusters': clusters,
                'silhouette': silhouette,
                'n_clusters': k
            }
            print(f"    K={k}: Silhouette score = {silhouette:.3f}")
    
    # Hierarchical clustering
    if 'hierarchical' in methods:
        print("  Running hierarchical clustering...")
        for k in [2, 3, 4, 5]:
            hierarchical = AgglomerativeClustering(n_clusters=k, linkage='ward')
            clusters = hierarchical.fit_predict(X_pca)
            silhouette = silhouette_score(X_pca, clusters)
            
            adata.obs[f'hierarchical_k{k}'] = clusters.astype(str)
            results[f'hierarchical_k{k}'] = {
                'clusters': clusters,
                'silhouette': silhouette,
                'n_clusters': k
            }
            print(f"    K={k}: Silhouette score = {silhouette:.3f}")
    
    # Find best clustering
    best_method = max(results.keys(), key=lambda x: results[x]['silhouette'])
    print(f"\n🏆 Best clustering: {best_method} (Silhouette = {results[best_method]['silhouette']:.3f})")
    
    return adata, results, best_method

def analyze_cluster_disease_relationship(adata, cluster_col, diagnosis_col='worst_diagnosis_label'):
    """Analyze relationship between clusters and disease labels."""
    if diagnosis_col not in adata.obs.columns:
        print(f"Warning: {diagnosis_col} not found in data")
        return None
    
    print(f"\n🔍 Analyzing relationship between {cluster_col} and {diagnosis_col}...")
    
    # Create contingency table
    contingency = pd.crosstab(adata.obs[cluster_col], adata.obs[diagnosis_col])
    print("\n📊 Cluster vs Diagnosis contingency table:")
    print(contingency)
    
    # Calculate percentages within each cluster
    cluster_purity = contingency.div(contingency.sum(axis=1), axis=0) * 100
    print("\n📈 Diagnosis composition within each cluster (%):")
    print(cluster_purity.round(1))
    
    # Calculate adjusted rand index if both are available
    if len(adata.obs[cluster_col].unique()) > 1 and len(adata.obs[diagnosis_col].unique()) > 1:
        # Convert to numeric for ARI calculation
        from sklearn.preprocessing import LabelEncoder
        le_cluster = LabelEncoder()
        le_diagnosis = LabelEncoder()
        
        cluster_numeric = le_cluster.fit_transform(adata.obs[cluster_col])
        diagnosis_numeric = le_diagnosis.fit_transform(adata.obs[diagnosis_col])
        
        ari = adjusted_rand_score(cluster_numeric, diagnosis_numeric)
        print(f"\n🎯 Adjusted Rand Index (cluster agreement with diagnosis): {ari:.3f}")
        
        return {
            'contingency': contingency,
            'cluster_purity': cluster_purity,
            'ari': ari
        }
    
    return {
        'contingency': contingency,
        'cluster_purity': cluster_purity
    }

def create_visualizations(adata, best_method, output_dir='adni_clustering_plots'):
    """Create comprehensive visualizations."""
    print(f"\n🎨 Creating visualizations...")
    
    Path(output_dir).mkdir(exist_ok=True)
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Figure 1: PCA plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # PCA by clusters
    sc.pl.pca(adata, color=best_method, ax=axes[0,0], show=False, frameon=False)
    axes[0,0].set_title(f'PCA colored by {best_method}', fontsize=14, fontweight='bold')
    
    # PCA by diagnosis
    if 'worst_diagnosis_label' in adata.obs.columns:
        sc.pl.pca(adata, color='worst_diagnosis_label', ax=axes[0,1], show=False, frameon=False)
        axes[0,1].set_title('PCA colored by Diagnosis', fontsize=14, fontweight='bold')
    
    # UMAP by clusters
    sc.pl.umap(adata, color=best_method, ax=axes[1,0], show=False, frameon=False)
    axes[1,0].set_title(f'UMAP colored by {best_method}', fontsize=14, fontweight='bold')
    
    # UMAP by diagnosis
    if 'worst_diagnosis_label' in adata.obs.columns:
        sc.pl.umap(adata, color='worst_diagnosis_label', ax=axes[1,1], show=False, frameon=False)
        axes[1,1].set_title('UMAP colored by Diagnosis', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/adni_clustering_overview.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/adni_clustering_overview.pdf', bbox_inches='tight')
    plt.close()
    
    # Figure 2: High-quality UMAP plots
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    # UMAP with clusters - larger points and better colors
    umap_data = adata.obsm['X_umap']
    clusters = adata.obs[best_method]
    
    # Plot clusters
    for cluster in clusters.unique():
        mask = clusters == cluster
        axes[0].scatter(umap_data[mask, 0], umap_data[mask, 1], 
                       s=50, alpha=0.7, label=f'Cluster {cluster}')
    
    axes[0].set_xlabel('UMAP 1', fontsize=12)
    axes[0].set_ylabel('UMAP 2', fontsize=12)
    axes[0].set_title(f'ADNI Samples - {best_method}', fontsize=16, fontweight='bold')
    axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[0].grid(True, alpha=0.3)
    
    # Plot diagnosis if available
    if 'worst_diagnosis_label' in adata.obs.columns:
        diagnosis = adata.obs['worst_diagnosis_label']
        diagnosis_colors = {'Cognitively Normal': '#2E8B57', 
                          'Mild Cognitive Impairment': '#FF8C00', 
                          'Alzheimer\'s Disease': '#DC143C'}
        
        for diag in diagnosis.unique():
            if pd.notna(diag):
                mask = diagnosis == diag
                color = diagnosis_colors.get(diag, '#808080')
                axes[1].scatter(umap_data[mask, 0], umap_data[mask, 1], 
                              s=50, alpha=0.7, label=diag, color=color)
        
        axes[1].set_xlabel('UMAP 1', fontsize=12)
        axes[1].set_ylabel('UMAP 2', fontsize=12)
        axes[1].set_title('ADNI Samples - Diagnosis Labels', fontsize=16, fontweight='bold')
        axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/adni_umap_detailed.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/adni_umap_detailed.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"✅ Visualizations saved to {output_dir}/")
    return output_dir

def generate_report(adata, results, best_method, cluster_analysis, output_dir):
    """Generate a comprehensive analysis report."""
    
    report = f"""
# ADNI Gene Expression Clustering Analysis Report

## Dataset Overview
- **Samples**: {adata.n_obs:,}
- **Genes**: {adata.n_vars:,} (highly variable genes used for clustering)
- **Date**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## Clustering Results

### Best Clustering Method: {best_method}
- **Silhouette Score**: {results[best_method]['silhouette']:.3f}
- **Number of Clusters**: {results[best_method]['n_clusters']}

### All Methods Comparison:
"""
    
    for method, result in results.items():
        report += f"- **{method}**: Silhouette = {result['silhouette']:.3f}\n"
    
    if cluster_analysis and 'ari' in cluster_analysis:
        report += f"""
## Cluster-Diagnosis Relationship
- **Adjusted Rand Index**: {cluster_analysis['ari']:.3f}
  - 0.0 = random clustering relative to diagnosis
  - 1.0 = perfect agreement with diagnosis
  - Values > 0.3 suggest meaningful relationship

### Diagnosis Distribution by Cluster:
"""
        
        for cluster in cluster_analysis['cluster_purity'].index:
            report += f"\n**Cluster {cluster}:**\n"
            for diagnosis, percentage in cluster_analysis['cluster_purity'].loc[cluster].items():
                report += f"- {diagnosis}: {percentage:.1f}%\n"
    
    report += f"""
## Files Generated
- `adni_clustering_overview.png/pdf`: Multi-panel overview
- `adni_umap_detailed.png/pdf`: High-quality UMAP plots
- `adni_clustering_report.md`: This report

## Interpretation Notes
1. **Cluster Separation**: Look for clear visual separation in UMAP plots
2. **Disease Association**: Check if clusters align with diagnosis categories
3. **Silhouette Score**: Values > 0.5 indicate good cluster separation
4. **Adjusted Rand Index**: Values > 0.3 suggest clusters relate to disease labels

## Next Steps
- Investigate genes driving cluster separation (differential expression analysis)
- Analyze other metadata factors (age, sex, APOE status if available)
- Consider more sophisticated clustering methods (leiden, louvain)
"""
    
    # Save report
    report_path = f"{output_dir}/adni_clustering_report.md"
    with open(report_path, 'w') as f:
        f.write(report)
    
    print(f"📄 Analysis report saved to: {report_path}")
    return report_path

def main():
    """Main analysis workflow."""
    print("🧬 ADNI Gene Expression Clustering Analysis")
    print("=" * 50)
    
    # Load data
    adata = load_adni_data()
    
    # Preprocess
    adata_processed = preprocess_for_clustering(adata, n_top_genes=2000)
    
    # Dimensionality reduction
    adata_processed = perform_dimensionality_reduction(adata_processed)
    
    # Clustering
    adata_processed, results, best_method = perform_clustering(adata_processed)
    
    # Analyze cluster-disease relationship
    cluster_analysis = analyze_cluster_disease_relationship(adata_processed, best_method)
    
    # Create visualizations
    output_dir = create_visualizations(adata_processed, best_method)
    
    # Generate report
    generate_report(adata_processed, results, best_method, cluster_analysis, output_dir)
    
    print(f"\n🎉 Analysis complete! Check {output_dir}/ for results.")
    
    return adata_processed, results, best_method

if __name__ == "__main__":
    adata_processed, results, best_method = main()