#!/usr/bin/env python3
"""
Heatmap Visualization for RNA-seq Analysis
-----------------------------------------

This module generates heatmap visualizations for RNA-seq data:
- Correlation heatmaps between datasets
- Expression heatmaps for genes
- Cluster heatmaps
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
from datetime import datetime

class HeatmapVisualizer:
    def __init__(self, config_file=None):
        """Initialize the heatmap visualizer with configuration"""
        # Load configuration
        if config_file is None:
            script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
            config_file = script_dir.parent / "config" / "settings.json"
            
        with open(config_file, 'r') as f:
            self.config = json.load(f)
            
        # Set up paths
        self.base_dir = Path(self.config['base_dir'])
        self.analysis_dir = Path(self.config['analysis_dir'])
        
        # Create output directory for visualizations
        self.viz_dir = self.analysis_dir / "visualizations"
        self.viz_dir.mkdir(parents=True, exist_ok=True)
        
        # Load housekeeping genes
        self.housekeeping_genes = self.config['housekeeping_genes']
        
        # Set default colormap
        self.cmap = "viridis"
    
    def _setup_dated_directory(self):
        """Set up a dated directory for visualizations"""
        today = datetime.now().strftime("%Y%m%d")
        dated_dir = self.viz_dir / today
        dated_dir.mkdir(exist_ok=True)
        return dated_dir
    
    def generate_correlation_heatmap(self, report, output_file=None, correlation_type='pearson', use_log=False):
        """Generate correlation heatmap from analysis report"""
        # Setup output directory
        if output_file is None:
            dated_dir = self._setup_dated_directory()
            output_file = dated_dir / f"correlation_heatmap_{correlation_type}.png"
        
        # Extract data
        cell_lines = list(report['dataset_summaries'].keys())
        
        # Create correlation matrix
        matrix_size = len(cell_lines)
        corr_matrix = np.ones((matrix_size, matrix_size))
        
        # Fill matrix with correlations
        for i, cell1 in enumerate(cell_lines):
            for j, cell2 in enumerate(cell_lines):
                if i == j:
                    # Diagonal - use replicate correlation if available
                    if report['dataset_summaries'][cell1].get('replicate_correlation'):
                        rep_corr = report['dataset_summaries'][cell1]['replicate_correlation']
                        corr_value = rep_corr.get(correlation_type, rep_corr.get('correlation', 1.0))
                        corr_matrix[i, j] = corr_value if corr_value is not None else 1.0
                elif i < j:
                    # Upper triangle
                    if cell1 in report['cross_correlations'] and cell2 in report['cross_correlations'][cell1]:
                        corr_data = report['cross_correlations'][cell1][cell2]
                        corr_value = corr_data.get(correlation_type, corr_data.get('correlation'))
                        if corr_value is not None:
                            corr_matrix[i, j] = corr_value
                            corr_matrix[j, i] = corr_value
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Plot heatmap
        mask = np.zeros_like(corr_matrix, dtype=bool)
        mask[np.triu_indices_from(mask, k=0)] = False  # Keep only upper triangle
        
        # Use RdBu_r colormap for correlations
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        
        sns.heatmap(corr_matrix, mask=mask, cmap=cmap, vmin=-1, vmax=1,
                  square=True, linewidths=.5, cbar_kws={"shrink": .5},
                  annot=True, fmt=".2f")
        
        # Add RNA type annotation
        rna_types = [report['dataset_summaries'][cell]['rna_type'] for cell in cell_lines]
        is_tissue = [report['dataset_summaries'][cell].get('is_tissue', False) for cell in cell_lines]
        
        # Create custom tick labels
        tick_labels = []
        for i, cell in enumerate(cell_lines):
            prefix = "T: " if is_tissue[i] else "C: "  # T for tissue, C for cell line
            tick_labels.append(f"{prefix}{cell} ({rna_types[i]})")
        
        # Apply labels
        plt.xticks(np.arange(len(cell_lines)) + 0.5, tick_labels, rotation=45, ha='right')
        plt.yticks(np.arange(len(cell_lines)) + 0.5, tick_labels, rotation=0)
        
        # Add title with transformation info
        transform_text = "Log-transformed " if use_log else ""
        plt.title(f'Dataset {transform_text}{correlation_type.capitalize()} Correlation Matrix')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        print(f"Gene expression heatmap saved to {output_file}")
        return output_file
    
    def generate_cluster_heatmap(self, cluster_report, output_file=None):
        """Generate heatmap showing expression clusters"""
        # Setup output directory
        if output_file is None:
            dated_dir = self._setup_dated_directory()
            output_file = dated_dir / "expression_clusters.png"
        
        # Extract data
        cell_lines = cluster_report['correlation_matrix']['cell_lines']
        matrix = np.array(cluster_report['correlation_matrix']['matrix'])
        
        # Get cluster assignments
        cluster_assignments = {}
        for i, cluster in enumerate(cluster_report['clusters']):
            cluster_id = cluster['id']
            for cell in cluster['cells']:
                cluster_assignments[cell] = cluster_id
                
        # Add outliers
        for outlier in cluster_report['outliers']:
            cluster_assignments[outlier['cell']] = 'Outlier'
        
        # Sort cell lines by cluster
        ordered_cells = []
        ordered_clusters = []
        
        # Add clustered cells first (in order of clusters)
        for cluster in cluster_report['clusters']:
            for cell in cluster['cells']:
                ordered_cells.append(cell)
                ordered_clusters.append(cluster['id'])
        
        # Add outliers at the end
        for outlier in cluster_report['outliers']:
            ordered_cells.append(outlier['cell'])
            ordered_clusters.append('Outlier')
        
        # Create reordered matrix
        indices = [cell_lines.index(cell) for cell in ordered_cells]
        reordered_matrix = matrix[np.ix_(indices, indices)]
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Plot heatmap
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        
        ax = sns.heatmap(reordered_matrix, cmap=cmap, vmin=-1, vmax=1,
                       square=True, linewidths=0.5, cbar_kws={"shrink": .5},
                       annot=True, fmt=".2f")
        
        # Set tick labels
        ax.set_xticks(np.arange(len(ordered_cells)) + 0.5)
        ax.set_yticks(np.arange(len(ordered_cells)) + 0.5)
        ax.set_xticklabels(ordered_cells, rotation=45, ha='right')
        ax.set_yticklabels(ordered_cells)
        
        # Add cluster boundaries
        current_cluster = ordered_clusters[0]
        boundaries = [0]
        
        for i, cluster in enumerate(ordered_clusters[1:], 1):
            if cluster != current_cluster:
                boundaries.append(i - 0.5)
                current_cluster = cluster
        
        # Draw lines for cluster boundaries
        for b in boundaries[1:]:
            plt.axhline(y=b, color='black', linestyle='-', linewidth=2)
            plt.axvline(x=b, color='black', linestyle='-', linewidth=2)
        
        # Add cluster labels to right side
        for i, cluster in enumerate(ordered_clusters):
            plt.text(len(ordered_cells) + 0.5, i + 0.5, cluster, 
                   ha='left', va='center', fontweight='bold')
        
        plt.title('Dataset Expression Clusters')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        print(f"Cluster heatmap saved to {output_file}")
        return output_file
    
    def generate_all_visualizations(self, report, cluster_report, expr_matrix, output_dir=None):
        """Generate all visualizations for a complete analysis"""
        # Setup output directory
        if output_dir is None:
            dated_dir = self._setup_dated_directory()
        else:
            dated_dir = Path(output_dir)
            dated_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate correlation heatmap
        self.generate_correlation_heatmap(
            report, 
            output_file=dated_dir / "correlation_heatmap.png",
            correlation_type='pearson',
            use_log=report['analysis_parameters'].get('use_log', False)
        )
        
        # Generate gene expression heatmap
        self.generate_gene_expression_heatmap(
            expr_matrix,
            output_file=dated_dir / "gene_expression_heatmap.png",
            use_log=report['analysis_parameters'].get('use_log', False)
        )
        
        # Generate cluster heatmap
        self.generate_cluster_heatmap(
            cluster_report,
            output_file=dated_dir / "expression_clusters.png"
        )
        
        print(f"All visualizations saved to {dated_dir}")
        return dated_dir
        plt.close()
        
        print(f"Correlation heatmap saved to {output_file}")
        return output_file
    
    def generate_gene_expression_heatmap(self, expr_matrix, output_file=None, use_log=False):
        """Generate heatmap of gene expression across datasets"""
        # Setup output directory
        if output_file is None:
            dated_dir = self._setup_dated_directory()
            output_file = dated_dir / "gene_expression_heatmap.png"
        
        # Extract data
        datasets = expr_matrix['dataset_id'].tolist()
        
        # Get gene columns (exclude metadata)
        gene_cols = [col for col in expr_matrix.columns if col in self.housekeeping_genes]
        
        # Create expression matrix for plotting
        plot_data = expr_matrix[gene_cols].values
        
        # Apply log transformation if requested
        if use_log:
            plot_data = np.log2(plot_data + 0.1)
        
        # Create figure
        plt.figure(figsize=(max(10, len(gene_cols) * 0.5), max(8, len(datasets) * 0.4)))
        
        # Plot heatmap
        cmap = "YlGnBu" if not use_log else "viridis"
        
        # Determine vmin and vmax
        if use_log:
            vmin = np.log2(0.1)  # log2 of minimum value
            vmax = np.max(plot_data)
        else:
            vmin = 0
            vmax = np.percentile(plot_data[plot_data > 0], 95)  # 95th percentile of non-zero values
        
        ax = sns.heatmap(plot_data, cmap=cmap, vmin=vmin, vmax=vmax,
                       xticklabels=gene_cols, yticklabels=datasets,
                       linewidths=0.5, cbar_kws={"shrink": .5})
        
        # Add RNA type and tissue annotations
        rna_types = expr_matrix['rna_type'].tolist()
        is_tissue = expr_matrix['is_tissue'].tolist()
        
        # Adjust y-tick labels
        new_labels = []
        for i, dataset in enumerate(datasets):
            prefix = "T: " if is_tissue[i] else "C: "  # T for tissue, C for cell line
            new_labels.append(f"{prefix}{dataset} ({rna_types[i]})")
            
        ax.set_yticklabels(new_labels)
        
        # Add title with transformation info
        transform_text = "Log-transformed " if use_log else ""
        plt.title(f'{transform_text}Gene Expression Across Datasets')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)