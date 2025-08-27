#!/usr/bin/env python3
"""
Expression Analyzer for RNA-seq Gene Quantification Data
-------------------------------------------------------

This module analyzes gene expression patterns, including:
- Calculating expression correlations between samples
- Identifying expression clusters
- Generating gene expression summaries
"""

import json
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os
from pathlib import Path

class ExpressionAnalyzer:
    def __init__(self, config_file=None):
        """Initialize the expression analyzer with configuration"""
        # Load configuration
        if config_file is None:
            script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
            config_file = script_dir.parent / "config" / "settings.json"
            
        with open(config_file, 'r') as f:
            self.config = json.load(f)
            
        # Set up paths
        self.base_dir = Path(self.config['base_dir'])
        self.analysis_dir = Path(self.config['analysis_dir'])
        
        # Ensure analysis directory exists
        self.analysis_dir.mkdir(parents=True, exist_ok=True)
        
        # Load housekeeping genes
        self.housekeeping_genes = self.config['housekeeping_genes']
        self.rna_seq_types = self.config['rna_seq_types']
    
    def find_common_housekeeping_genes(self, results1, results2, rep1='rep1', rep2='rep1', tpm_threshold=1.0):
        """Find housekeeping genes present in both datasets above threshold"""
        if rep1 not in results1 or rep2 not in results2:
            return []
            
        genes1 = set(gene for gene in results1[rep1].keys() 
                    if gene in self.housekeeping_genes and results1[rep1][gene] >= tpm_threshold)
        genes2 = set(gene for gene in results2[rep2].keys() 
                    if gene in self.housekeeping_genes and results2[rep2][gene] >= tpm_threshold)
        
        common_genes = genes1 & genes2
        
        if not common_genes:
            print(f"Warning: No common housekeeping genes found between datasets")
        else:
            print(f"Found {len(common_genes)} common housekeeping genes")
            
        return list(common_genes)
    
    def calculate_replicate_correlation(self, cell_results):
        """Calculate correlation between replicates using only common genes"""
        replicates = cell_results.get('replicates', {})
        rep_keys = [rep for rep in replicates.keys() if rep.startswith('rep')]
        
        if len(rep_keys) < 2:
            # Need at least 2 replicates
            return None
                
        # Get values for first two replicates
        rep1 = rep_keys[0]
        rep2 = rep_keys[1]
        
        # Find common housekeeping genes
        common_genes = self.find_common_housekeeping_genes(replicates, replicates, rep1, rep2)
        
        if len(common_genes) < 3:
            print(f"Too few common genes for correlation calculation")
            return None
            
        values1 = [replicates[rep1].get(gene, 0) for gene in common_genes]
        values2 = [replicates[rep2].get(gene, 0) for gene in common_genes]
        
        # Calculate Pearson and Spearman correlations
        pearson_corr = np.corrcoef(values1, values2)[0, 1]
        spearman_corr = spearmanr(values1, values2).correlation
        
        # Return correlation along with common gene count
        return {
            'correlation': pearson_corr,  # For backward compatibility
            'pearson': pearson_corr,
            'spearman': spearman_corr,
            'common_genes': len(common_genes),
            'total_genes': len(self.housekeeping_genes),
            'percent_common': len(common_genes) / len(self.housekeeping_genes) * 100
        }
    
    def calculate_expression_correlation(self, expr1, expr2, tpm_threshold=1.0, use_log=False):
        """Calculate correlations between expression profiles with option for log transformation"""
        # Find genes present in both profiles above threshold
        common_genes = [gene for gene in self.housekeeping_genes 
                      if gene in expr1 and gene in expr2 
                      and expr1[gene] >= tpm_threshold and expr2[gene] >= tpm_threshold]
        
        if len(common_genes) < 3:
            print(f"Too few common genes for correlation: {len(common_genes)}")        
            return {
                'pearson': None,
                'spearman': None,
                'correlation': None,  # For backward compatibility
                'common_genes': len(common_genes),
                'percent_common': len(common_genes) / len(self.housekeeping_genes) * 100
            }
                
        values1 = [expr1[gene] for gene in common_genes]
        values2 = [expr2[gene] for gene in common_genes]

        # Apply log transformation if requested
        if use_log:
            # Add small value to avoid log(0)
            values1 = np.log2([v + 0.1 for v in values1])
            values2 = np.log2([v + 0.1 for v in values2])
                
        # Calculate both correlation types
        pearson_corr = np.corrcoef(values1, values2)[0, 1]
        spearman_corr = spearmanr(values1, values2).correlation
        
        return {
            'pearson': pearson_corr,
            'spearman': spearman_corr,
            'correlation': pearson_corr,  # For backward compatibility
            'log_transformed': use_log,
            'common_genes': len(common_genes),
            'percent_common': len(common_genes) / len(self.housekeeping_genes) * 100
        }
        
    def calculate_cell_correlations(self, all_results, use_log=False, correlation_type='pearson'):
        """Calculate correlations between all cell lines"""
        cell_lines = list(all_results.keys())
        correlations = {}
        
        for i, cell1 in enumerate(cell_lines):
            correlations[cell1] = {}
            
            cell1_data = all_results[cell1]
            cell1_rna_type = cell1_data['metadata']['rna_type']
            cell1_avg = self._calculate_avg_expression(cell1_data)
            
            for j, cell2 in enumerate(cell_lines):
                if i == j:
                    continue
                    
                cell2_data = all_results[cell2]
                cell2_rna_type = cell2_data['metadata']['rna_type']
                cell2_avg = self._calculate_avg_expression(cell2_data)
                
                # Calculate correlation with common genes info
                corr_data = self.calculate_expression_correlation(cell1_avg, cell2_avg, use_log=use_log)
                
                correlations[cell1][cell2] = {
                    'correlation': corr_data[correlation_type],  # Use the requested correlation type
                    'pearson': corr_data['pearson'],
                    'spearman': corr_data['spearman'],
                    'common_genes': corr_data.get('common_genes', 0),
                    'percent_common': corr_data.get('percent_common', 0),
                    'cell1_rna_type': cell1_rna_type,
                    'cell2_rna_type': cell2_rna_type,
                    'same_rna_type': cell1_rna_type == cell2_rna_type,
                    'log_transformed': use_log
                }
            
        return correlations
      
    def _calculate_avg_expression(self, cell_data):
        """Calculate average expression across replicates"""
        replicates = cell_data.get('replicates', {})
        
        if not replicates:
            return {}
            
        avg_expression = {}
        
        for gene in self.housekeeping_genes:
            values = []
            for rep_data in replicates.values():
                if gene in rep_data:
                    values.append(rep_data[gene])
            
            if values:
                avg_expression[gene] = np.mean(values)
            else:
                avg_expression[gene] = 0
                
        return avg_expression
    
    def generate_expression_matrix(self, all_results):
        """Generate a matrix of gene expression for all datasets"""
        # Get unique list of genes across all datasets
        all_genes = set(self.housekeeping_genes)
        
        # Prepare data for matrix
        matrix_data = []
        
        for cell_line, results in all_results.items():
            # Calculate average expression for this cell line
            avg_expr = self._calculate_avg_expression(results)
            
            # Add metadata
            row_data = {
                'dataset_id': cell_line,
                'rna_type': results['metadata']['rna_type'],
                'is_tissue': results['metadata'].get('is_tissue', False)
            }
            
            # Add expression values
            for gene in all_genes:
                row_data[gene] = avg_expr.get(gene, 0)
            
            matrix_data.append(row_data)
        
        # Convert to DataFrame
        expr_matrix = pd.DataFrame(matrix_data)
        
        return expr_matrix
    
    def identify_expression_clusters(self, all_results, min_correlation=0.9, use_log=False, correlation_type='pearson'):
        """Identify clusters of datasets with similar expression patterns"""
        # Create expression matrix
        expr_matrix = self.generate_expression_matrix(all_results)
        
        # Get gene columns (exclude metadata)
        gene_cols = [col for col in expr_matrix.columns if col in self.housekeeping_genes]
        
        # Calculate correlation matrix
        cell_lines = expr_matrix['dataset_id'].tolist()
        n_cells = len(cell_lines)
        corr_matrix = np.zeros((n_cells, n_cells))
        
        # Calculate correlation for each pair
        for i in range(n_cells):
            for j in range(i, n_cells):
                if i == j:
                    corr_matrix[i, j] = 1.0
                else:
                    # Get expression vectors
                    vec1 = expr_matrix.iloc[i][gene_cols].values
                    vec2 = expr_matrix.iloc[j][gene_cols].values
                    
                    # Apply log transformation if requested
                    if use_log:
                        # Add small value to avoid log(0)
                        vec1_log = np.log2(vec1 + 0.1)
                        vec2_log = np.log2(vec2 + 0.1)
                        
                        # Calculate correlation
                        if correlation_type == 'spearman':
                            corr = spearmanr(vec1_log, vec2_log).correlation
                        else:
                            corr = np.corrcoef(vec1_log, vec2_log)[0, 1]
                    else:
                        # Calculate correlation without log transform
                        if correlation_type == 'spearman':
                            corr = spearmanr(vec1, vec2).correlation
                        else:
                            corr = np.corrcoef(vec1, vec2)[0, 1]
                    
                    corr_matrix[i, j] = corr
                    corr_matrix[j, i] = corr
        
        # Clustering algorithm (basic approach)
        clusters = []
        unassigned = set(range(n_cells))
        
        while unassigned:
            # Start a new cluster with first unassigned cell
            current = min(unassigned)
            current_cluster = {current}
            unassigned.remove(current)
            
            # Find cells that correlate well with current cluster
            changed = True
            while changed and unassigned:
                changed = False
                to_add = set()
                
                for cell_idx in unassigned:
                    # Check correlation with all cells in current cluster
                    correlates_with_all = True
                    for cluster_idx in current_cluster:
                        if corr_matrix[cell_idx, cluster_idx] < min_correlation:
                            correlates_with_all = False
                            break
                    
                    if correlates_with_all:
                        to_add.add(cell_idx)
                
                # Add correlating cells to cluster
                if to_add:
                    current_cluster.update(to_add)
                    unassigned.difference_update(to_add)
                    changed = True
            
            # Convert indices to cell line names
            name_cluster = [cell_lines[idx] for idx in current_cluster]
            clusters.append(name_cluster)
        
        # Create cluster report
        cluster_report = {
            'clusters': [],
            'outliers': [],
            'correlation_matrix': {
                'cell_lines': cell_lines,
                'matrix': corr_matrix.tolist()
            }
        }
        
        # Label each cluster and identify key characteristics
        for i, cluster in enumerate(clusters):
            if len(cluster) >= 2:
                # Get RNA types in this cluster
                rna_types = {}
                for cell in cluster:
                    row_idx = cell_lines.index(cell)
                    cell_type = expr_matrix.iloc[row_idx]['rna_type']
                    rna_types[cell_type] = rna_types.get(cell_type, 0) + 1
                
                cluster_report['clusters'].append({
                    'id': f"Group {chr(65+i)}",  # A, B, C, etc.
                    'cells': cluster,
                    'size': len(cluster),
                    'rna_types': rna_types
                })
            else:
                # Single cell "outliers"
                row_idx = cell_lines.index(cluster[0])
                cluster_report['outliers'].append({
                    'cell': cluster[0],
                    'rna_type': expr_matrix.iloc[row_idx]['rna_type']
                })
        
        return cluster_report
    
    def generate_gene_summary(self, all_results, use_log=False):
        """Generate summary of gene expression across datasets"""
        # Generate expression matrix
        expr_matrix = self.generate_expression_matrix(all_results)
        
        # Summarize for each gene
        summary = {}
        
        for gene in self.housekeeping_genes:
            gene_summary = {
                'expression_by_cell': {},
                'expression_by_rna_type': {rna_type: [] for rna_type in self.rna_seq_types},
                'variability': None,
                'consistency': None
            }
            
            # Get expression values for this gene
            values = []
            
            for i, row in expr_matrix.iterrows():
                cell_line = row['dataset_id']
                rna_type = row['rna_type']
                expr_value = row[gene]
                
                # Add to dataset-specific summary
                gene_summary['expression_by_cell'][cell_line] = {
                    'mean': expr_value,
                    'rna_type': rna_type
                }
                
                # Add to RNA type group
                gene_summary['expression_by_rna_type'][rna_type].append(expr_value)
                values.append(expr_value)
            
            # When calculating variability metrics
            if len(values) >= 3:
                # Original CV calculation
                original_cv = np.std(values) / np.mean(values) if np.mean(values) > 0 else None
                
                # Log-transformed CV calculation
                log_values = np.log2([v + 0.1 for v in values])
                log_cv = np.std(log_values) / np.mean(log_values) if np.mean(log_values) > 0 else None
                
                gene_summary['variability'] = {
                    'mean': float(np.mean(values)),
                    'std': float(np.std(values)),
                    'cv': float(original_cv) if original_cv is not None else None,
                    'log_mean': float(np.mean(log_values)),
                    'log_std': float(np.std(log_values)),
                    'log_cv': float(log_cv) if log_cv is not None else None
                }
                
                # Assess consistency
                cv = log_cv if use_log else original_cv
                if cv is not None:
                    if cv < 0.2:
                        gene_summary['consistency'] = "high"
                    elif cv < 0.5:
                        gene_summary['consistency'] = "medium"
                    else:
                        gene_summary['consistency'] = "low"
            
            # Add this gene's summary to the overall summary dictionary
            summary[gene] = gene_summary
        
        return summary
    
    def generate_comprehensive_report(self, all_results, use_log=False, correlation_type='pearson'):
        """Generate a comprehensive analysis report"""
        # Calculate correlations between datasets
        cross_correlations = self.calculate_cell_correlations(
            all_results, use_log=use_log, correlation_type=correlation_type
        )
        
        # Generate gene expression summary
        gene_summary = self.generate_gene_summary(all_results, use_log=use_log)
        
        # Group by RNA type
        rna_type_groups = {}
        for cell_line, results in all_results.items():
            rna_type = results['metadata']['rna_type']
            if rna_type not in rna_type_groups:
                rna_type_groups[rna_type] = []
            rna_type_groups[rna_type].append(cell_line)
        
        # Dataset summaries
        dataset_summaries = {}
        for cell_line, results in all_results.items():
            # Calculate replicate correlation if available
            rep_corr = self.calculate_replicate_correlation(results)
            
            dataset_summaries[cell_line] = {
                'replicate_correlation': rep_corr,
                'has_replicates': len(results['replicates']) > 1,
                'rna_type': results['metadata']['rna_type'],
                'is_tissue': results['metadata'].get('is_tissue', False)
            }
        
        # Final report
        report = {
            'dataset_summaries': dataset_summaries,
            'cross_correlations': cross_correlations,
            'gene_expression_summary': gene_summary,
            'rna_type_groups': rna_type_groups,
            'analysis_parameters': {
                'use_log': use_log,
                'correlation_type': correlation_type
            }
        }
        
        # Save report
        report_file = self.analysis_dir / "expression_analysis_report.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        return report