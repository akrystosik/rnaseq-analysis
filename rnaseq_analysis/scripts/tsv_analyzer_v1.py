#!/bin/bash
# RNA-seq TSV Analysis Installation Script
# This script sets up the TSV analyzer for comparing RNA-seq datasets

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis"
SCRIPT_DIR="${BASE_DIR}/scripts"


#!/usr/bin/env python3
"""
TSV Analyzer for RNA-seq Gene Quantification Data
-------------------------------------------------

This script analyzes gene quantification TSV files from ENCODE RNA-seq experiments
to identify comparable datasets across cell lines. Unlike the bigWig-based approach,
this directly uses the processed gene-level quantifications (TPM values) which:

1. Is more efficient - eliminates processing raw signal tracks
2. Is more direct - works with the same gene quantifications used for genomic integration
3. Produces more consistent results - avoids technical variations in signal extraction
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import requests
import json
import argparse
from pathlib import Path
import re
from collections import defaultdict

class TSVAnalyzer:
    def __init__(self, base_dir="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis"):
        self.base_dir = Path(base_dir)
        self.raw_data_dir = self.base_dir / "raw_data/gene_quantification"
        self.analysis_dir = self.base_dir / "analysis"
        self.metadata_dir = self.base_dir / "metadata"
        
        # Create directories
        self.raw_data_dir.mkdir(parents=True, exist_ok=True)
        self.analysis_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_dir.mkdir(parents=True, exist_ok=True)
        
        # Define housekeeping genes to focus on
        self.housekeeping_genes = [
            'ACTB', 'GAPDH', 'TUBB', 'B2M', 'PPIA', 'TBP', 'HPRT1', 'RPL13A'
        ]
        
        # RNA-seq types
        self.rna_seq_types = ['total', 'polyA+', 'polyA-']
        
        # Define cell line datasets (experiment IDs)
        self.cell_lines = {
            'K562': {'experiment_id': 'ENCSR792OIJ', 'rna_type': 'total'},
            'Panc1': {'experiment_id': 'ENCSR128CYL', 'rna_type': 'total'},
            'NCI-H460': {'experiment_id': 'ENCSR164OCT', 'rna_type': 'total'},
            'A549': {'experiment_id': 'ENCSR414IGI', 'rna_type': 'total'},
            'SK-N-MC': {'experiment_id': '', 'rna_type': 'total'},
            'K562_AEL': {'experiment_id': 'ENCSR000AEL', 'rna_type': 'total'},
            'HepG2_total': {'experiment_id': 'ENCSR245ATJ', 'rna_type': 'total'},
            'HepG2_polyA_plus': {'experiment_id': 'ENCSR000CPD', 'rna_type': 'polyA+'},
            'HepG2_polyA_minus': {'experiment_id': 'ENCSR000CPE', 'rna_type': 'polyA-'},
            'HepG2_total_ZGR': {'experiment_id': 'ENCSR181ZGR', 'rna_type': 'total'},
            'Caki2_584JXD': {'experiment_id': 'ENCSR584JXD', 'rna_type': 'total'},
            'A549_polyA_minus': {'experiment_id': 'ENCSR000CQC', 'rna_type': 'polyA-'},
            'A549_polyA_plus': {'experiment_id': 'ENCSR000CON', 'rna_type': 'polyA+'},
            'GM23248_BPP': {'experiment_id': 'ENCSR797BPP', 'rna_type': 'total'}
        }
                
        # Gene ID pattern regex for matching gene symbols
        self.gene_id_pattern = re.compile(r'(ENSG\d+\.\d+)\|(\w+)')
        
        
    def get_tsv_files(self, experiment_id):
        """Find gene quantification TSV files for a given experiment ID"""
        # Fetch metadata if not already cached
        metadata = self.fetch_experiment_metadata(experiment_id)
        
        if not metadata:
            return []
        
        # Look for gene quantification files
        tsv_files = []
        for file_info in metadata.get('files', []):
            if (file_info.get('output_type') == 'gene quantifications' and 
                file_info.get('file_format') == 'tsv'):
                
                rep_info = file_info.get('biological_replicates', [])
                rep_num = rep_info[0] if rep_info else 0
                
                tsv_files.append({
                    'accession': file_info.get('accession'),
                    'url': f"https://www.encodeproject.org/files/{file_info.get('accession')}/@@download/{file_info.get('accession')}.tsv",
                    'rep': rep_num
                })
        
        return tsv_files
    
    def fetch_experiment_metadata(self, experiment_id):
        """Fetch experiment metadata from ENCODE"""
        metadata_file = self.metadata_dir / f"{experiment_id}.json"
        
        # Try to load from cache first
        if metadata_file.exists():
            try:
                with open(metadata_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Error loading cached metadata: {e}")
        
        # Fetch from ENCODE API
        try:
            print(f"Fetching metadata for experiment {experiment_id}...")
            url = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
            response = requests.get(url)
            if response.status_code == 200:
                metadata = response.json()
                # Save to cache
                with open(metadata_file, 'w') as f:
                    json.dump(metadata, f, indent=2)
                return metadata
            else:
                print(f"Error fetching metadata: HTTP {response.status_code}")
        except Exception as e:
            print(f"Error fetching metadata: {e}")
        
        return None
    
    def download_tsv(self, accession, url, cell_line):
        """Download a TSV file"""
        output_dir = self.raw_data_dir / cell_line
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{accession}.tsv"
        
        if output_file.exists():
            print(f"File {accession}.tsv already exists for {cell_line}")
            return output_file
        
        try:
            print(f"Downloading {accession}.tsv for {cell_line}...")
            response = requests.get(url)
            if response.status_code == 200:
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded {accession}.tsv")
                return output_file
            else:
                print(f"Error downloading file: HTTP {response.status_code}")
        except Exception as e:
            print(f"Error downloading {accession}: {str(e)}")
        
        return None
    
    def parse_tsv(self, file_path):
        """Parse a gene quantification TSV file"""
        try:
            df = pd.read_csv(file_path, sep='\t')
            
            # Check columns
            required_cols = ['gene_id', 'TPM']
            if not all(col in df.columns for col in required_cols):
                # Try to find alternative columns
                if 'gene' in df.columns and 'TPM' in df.columns:
                    df = df.rename(columns={'gene': 'gene_id'})
                elif 'gene_id' in df.columns and 'FPKM' in df.columns and 'TPM' not in df.columns:
                    # Some older ENCODE files have FPKM but not TPM
                    df = df.rename(columns={'FPKM': 'TPM'})
                    print(f"Using FPKM as TPM for {file_path}")
                else:
                    print(f"Required columns not found in {file_path}")
                    print(f"Available columns: {df.columns.tolist()}")
                    return None
            
            return df[['gene_id', 'TPM']]
            
        except Exception as e:
            print(f"Error parsing {file_path}: {str(e)}")
        
        return None
    
    def extract_gene_data(self, df, gene_symbol):
        """Extract data for a specific gene from the dataframe using Ensembl IDs"""
        if df is None:
            return 0
        
        # Map of gene symbols to their Ensembl IDs
        gene_ensembl_map = {
            'ACTB': 'ENSG00000075624',
            'GAPDH': 'ENSG00000111640',
            'TUBB': 'ENSG00000196230',
            'B2M': 'ENSG00000166710',
            'PPIA': 'ENSG00000196262',
            'TBP': 'ENSG00000112592',
            'HPRT1': 'ENSG00000165704',
            'RPL13A': 'ENSG00000142541'
        }
        
        # Try using the Ensembl ID if we have it
        if gene_symbol in gene_ensembl_map:
            ensembl_id = gene_ensembl_map[gene_symbol]
            # Match the Ensembl ID, including the version number (if any)
            mask = df['gene_id'].str.startswith(ensembl_id, na=False)
            matched = df[mask]
            
            if not matched.empty:
                return matched['TPM'].values[0]
        
        # Fallbacks if the ensembl ID approach doesn't work
        
        # Try direct string match with gene symbol
        mask = df['gene_id'].str.contains(gene_symbol, case=False, na=False)
        matched = df[mask]
        
        if not matched.empty:
            return matched['TPM'].values[0]
        
        # Try with transcript ID column
        if 'transcript_id(s)' in df.columns:
            mask = df['transcript_id(s)'].str.contains(gene_symbol, case=False, na=False)
            matched = df[mask]
            
            if not matched.empty:
                return matched['TPM'].values[0]
        
        print(f"Gene {gene_symbol} not found in dataframe")
        return 0    
    
    def analyze_cell_line(self, cell_line):
        """Analyze gene expression for a cell line"""
        if cell_line not in self.cell_lines:
            print(f"Cell line {cell_line} not found")
            return None
        
        experiment_id = self.cell_lines[cell_line]['experiment_id']
        tsv_files = self.get_tsv_files(experiment_id)
        
        if not tsv_files:
            print(f"No TSV files found for {cell_line} ({experiment_id})")
            return None
        
        results = {
            'metadata': {
                'experiment_id': experiment_id,
                'rna_type': self.cell_lines[cell_line]['rna_type']
            },
            'replicates': {}
        }
        
        # Process each TSV file
        for file_info in tsv_files:
            accession = file_info['accession']
            url = file_info['url']
            rep = f"rep{file_info['rep']}"
            
            # Download and parse TSV
            file_path = self.download_tsv(accession, url, cell_line)
            if not file_path:
                continue
                
            df = self.parse_tsv(file_path)
            if df is None:
                continue
            
            # Extract housekeeping gene data
            replicate_data = {}
            for gene in self.housekeeping_genes:
                tpm = self.extract_gene_data(df, gene)
                replicate_data[gene] = tpm
            
            results['replicates'][rep] = replicate_data
        
        return results
    
    def analyze_all_cell_lines(self):
        """Analyze all cell lines"""
        all_results = {}
        
        for cell_line in self.cell_lines:
            print(f"\nAnalyzing {cell_line}...")
            results = self.analyze_cell_line(cell_line)
            if results and results['replicates']:
                all_results[cell_line] = results
        
        return all_results

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
    # 2. Update the calculate_replicate_correlation method
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
        
        # Return correlation along with common gene count
        corr = np.corrcoef(values1, values2)[0, 1]
        return {
            'correlation': corr,
            'common_genes': len(common_genes),
            'total_genes': len(self.housekeeping_genes),
            'percent_common': len(common_genes) / len(self.housekeeping_genes) * 100
        }


    def _calculate_expression_correlation(self, expr1, expr2, tpm_threshold=1.0, use_log=False):
        """Calculate Pearson and Spearman correlations between expression profiles with option for log transformation"""
        # Find genes present in both profiles above threshold
        common_genes = [gene for gene in self.housekeeping_genes 
                    if gene in expr1 and gene in expr2 
                    and expr1[gene] >= tpm_threshold and expr2[gene] >= tpm_threshold]
        
        if len(common_genes) < 3:
            print(f"Too few common genes for correlation: {len(common_genes)}")        
            return {
                'pearson': None,
                'spearman': None,
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
                
        from scipy.stats import spearmanr
        
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
                corr_data = self._calculate_expression_correlation(cell1_avg, cell2_avg, use_log=use_log)
                
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
    

    def identify_expression_clusters(self, all_results, min_correlation=0.9, use_log=False, correlation_type='pearson'):
        """Identify clusters of cell lines with similar expression patterns"""
        # Calculate normalized data
        normalized_data = self._generate_normalized_data(all_results)
        
        # Create cell line expression vectors
        expression_vectors = {}
        
        # Create expression vectors using housekeeping genes
        for cell_line, data in normalized_data.items():
            vector = []
            for gene in self.housekeeping_genes:
                if gene in data['normalized_expression']:
                    vector.append(data['normalized_expression'][gene])
                else:
                    vector.append(0)
            expression_vectors[cell_line] = np.array(vector)
        
        # Calculate correlation matrix
        cell_lines = sorted(expression_vectors.keys())
        n_cells = len(cell_lines)
        corr_matrix = np.zeros((n_cells, n_cells))

        #  calculating correlations
        for i in range(n_cells):
            for j in range(i, n_cells):
                if i == j:
                    corr_matrix[i, j] = 1.0
                else:
                    vec1 = expression_vectors[cell_lines[i]]
                    vec2 = expression_vectors[cell_lines[j]]
                    
                    # Apply log transformation if requested
                    if use_log:
                        # Add small value to avoid log(0)
                        vec1_log = np.log2(vec1 + 0.1)
                        vec2_log = np.log2(vec2 + 0.1)
                        
                        if correlation_type == 'spearman':
                            from scipy.stats import spearmanr
                            corr = spearmanr(vec1_log, vec2_log).correlation
                        else:
                            corr = np.corrcoef(vec1_log, vec2_log)[0, 1]
                    else:
                        if correlation_type == 'spearman':
                            from scipy.stats import spearmanr
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
                    cell_type = all_results[cell]['metadata']['rna_type']
                    rna_types[cell_type] = rna_types.get(cell_type, 0) + 1
                
                cluster_report['clusters'].append({
                    'id': f"Group {chr(65+i)}",  # A, B, C, etc.
                    'cells': cluster,
                    'size': len(cluster),
                    'rna_types': rna_types
                })
            else:
                # Single cell "outliers"
                cluster_report['outliers'].append({
                    'cell': cluster[0],
                    'rna_type': all_results[cluster[0]]['metadata']['rna_type']
                })
        
        return cluster_report
    
    def _generate_normalized_data(self, all_results):
        """Generate normalized gene expression data"""
        normalized_data = {}
        
        for cell_line, results in all_results.items():
            normalized_data[cell_line] = {
                'metadata': results['metadata'],
                'normalized_expression': {}
            }
            
            # Calculate average expression for each gene
            for gene in self.housekeeping_genes:
                values = []
                for rep_data in results['replicates'].values():
                    if gene in rep_data:
                        values.append(rep_data[gene])
                
                if values:
                    normalized_data[cell_line]['normalized_expression'][gene] = np.mean(values)
                else:
                    normalized_data[cell_line]['normalized_expression'][gene] = 0
        
        return normalized_data
    
    def _group_by_rna_type(self, all_results):
        """Group cell lines by RNA-seq type"""
        groups = {rna_type: [] for rna_type in self.rna_seq_types}
        
        for cell_line, results in all_results.items():
            rna_type = results['metadata']['rna_type']
            groups[rna_type].append(cell_line)
            
        return groups
    
    def generate_gene_summary(self, all_results, use_log=False):
        """Generate summary of gene expression across cell lines"""
        summary = {}
        
        for gene in self.housekeeping_genes:
            gene_summary = {
                'expression_by_cell': {},
                'expression_by_rna_type': {rna_type: [] for rna_type in self.rna_seq_types},
                'variability': None,
                'consistency': None
            }
            
            # Calculate average expression for each cell line
            values = []
            
            for cell_line, results in all_results.items():
                rna_type = results['metadata']['rna_type']
                
                # Calculate average across replicates
                gene_values = []
                for rep_data in results['replicates'].values():
                    if gene in rep_data:
                        gene_values.append(rep_data[gene])
                
                if gene_values:
                    avg_value = np.mean(gene_values)
                    gene_summary['expression_by_cell'][cell_line] = {
                        'mean': avg_value,
                        'rna_type': rna_type
                    }
                    
                    # Add to RNA type group
                    gene_summary['expression_by_rna_type'][rna_type].append(avg_value)
                    values.append(avg_value)
            
            # When calculating variability metrics
            if len(values) >= 3:
                # Original CV calculation
                original_cv = np.std(values) / np.mean(values) if np.mean(values) > 0 else None
                
                # Log-transformed CV calculation
                log_values = np.log2([v + 0.1 for v in values])
                log_cv = np.std(log_values) / np.mean(log_values) if np.mean(log_values) > 0 else None
                
                gene_summary['variability'] = {
                    'mean': np.mean(values),
                    'std': np.std(values),
                    'cv': original_cv,
                    'log_mean': np.mean(log_values),
                    'log_std': np.std(log_values),
                    'log_cv': log_cv
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
    
    def generate_report(self, all_results, use_log=False, correlation_type='pearson'):
        """Generate comprehensive analysis report"""
        report = {
            'cell_line_data': {},
            'replicate_correlations': {},
            'cross_cell_correlations': self.calculate_cell_correlations(all_results, use_log=use_log, correlation_type=correlation_type),
            'gene_expression_summary': self.generate_gene_summary(all_results, use_log=use_log),
            'rna_type_groups': self._group_by_rna_type(all_results),
            'analysis_parameters': {
                'use_log': use_log,
                'correlation_type': correlation_type
            }
        }
        
        # Process each cell line
        for cell_line, results in all_results.items():
            # Calculate replicate correlation
            rep_corr = self.calculate_replicate_correlation(results)
            report['replicate_correlations'][cell_line] = rep_corr
            
            # Get metadata
            metadata = results['metadata']
            
            # Store cell line data
            report['cell_line_data'][cell_line] = {
                'has_replicates': len(results['replicates']) > 1,
                'rna_type': metadata['rna_type'],
                'experiment_id': metadata['experiment_id']
            }
        
        # Save report
        output_file = self.analysis_dir / "tsv_analysis_report.json"
        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        return report    

def print_report_summary(report):
    """Print a summary of the analysis report"""
    print("\nSummary Report")
    print("==============")
    
    # Replicate correlations
    print("\nReplicate Correlations:")
    for cell_line, corr_data in report['replicate_correlations'].items():
        if isinstance(corr_data, dict) and corr_data.get('correlation') is not None:
            print(f"  {cell_line}: r={corr_data['correlation']:.3f} ({corr_data['common_genes']}/{corr_data['total_genes']} genes, {corr_data['percent_common']:.1f}%)")
        elif isinstance(corr_data, dict):
            print(f"  {cell_line}: N/A (insufficient common genes: {corr_data.get('common_genes', 0)})")
        else:
            print(f"  {cell_line}: N/A (single replicate)")
    
    # RNA-seq type groups
    print("\nCell Lines by RNA-seq Type:")
    for rna_type, cells in report['rna_type_groups'].items():
        if cells:
            print(f"  {rna_type}: {', '.join(cells)}")
    
    # Cross-cell correlations
    print("\nCross-Cell Line Correlations (Same RNA Type Only):")
    for cell1, correlations in report['cross_cell_correlations'].items():
        for cell2, corr_data in correlations.items():
            if corr_data['same_rna_type'] and corr_data.get('correlation') is not None:
                print(f"  {cell1} vs {cell2}: r={corr_data['correlation']:.3f} ({corr_data['common_genes']} common genes, {corr_data['percent_common']:.1f}%) ({corr_data['cell1_rna_type']})")
    
    # Gene consistency
    print("\nHousekeeping Gene Consistency:")
    for gene, data in report['gene_expression_summary'].items():
        if data['consistency'] is not None:
            print(f"  {gene}: {data['consistency'].upper()} consistency")
            if data['variability'] is not None:
                cv = data['variability']['cv']
                if cv is not None:
                    print(f"    CV: {cv:.3f}")
    
    # Cell line recommendations
    print("\nCell Line Dataset Recommendations:")
    for cell_line, cell_data in report['cell_line_data'].items():
        has_replicates = cell_data['has_replicates']
        rep_corr_data = report['replicate_correlations'].get(cell_line)
        rna_type = cell_data['rna_type']
        
        # Extract correlation value if it's in dictionary format
        rep_corr = None
        if isinstance(rep_corr_data, dict):
            rep_corr = rep_corr_data.get('correlation')
        else:
            rep_corr = rep_corr_data
            
        if has_replicates and rep_corr is not None and rep_corr > 0.9:
            print(f"  {cell_line}: HIGH CONFIDENCE (r={rep_corr:.3f}, {rna_type})")
        elif has_replicates and rep_corr is not None and rep_corr > 0.7:
            print(f"  {cell_line}: MEDIUM CONFIDENCE (r={rep_corr:.3f}, {rna_type})")
        elif not has_replicates:
            print(f"  {cell_line}: LOW CONFIDENCE (single replicate, {rna_type})")
        else:
            print(f"  {cell_line}: LOW CONFIDENCE (r={rep_corr:.3f}, {rna_type})")
            
def generate_heatmap(report, output_file, correlation_type='pearson', use_log=False):
    """Generate correlation heatmap with RNA type info"""
    # Collect all cell lines
    cell_lines = list(report['replicate_correlations'].keys())
    rna_types = [report['cell_line_data'][cell]['rna_type'] for cell in cell_lines]
    
    # Create correlation matrix
    matrix_size = len(cell_lines)
    corr_matrix = np.ones((matrix_size, matrix_size))
    rna_same_matrix = np.ones((matrix_size, matrix_size), dtype=bool)
    
    # Fill in cross-correlations
    for i, cell1 in enumerate(cell_lines):
        for j, cell2 in enumerate(cell_lines):
            if i == j:
                # Diagonal: use replicate correlation or 1 if single replicate
                rep_corr_data = report['replicate_correlations'].get(cell1)
                
                # Check if rep_corr_data is a dictionary (new format) or float/None (old format)
                if isinstance(rep_corr_data, dict):
                    if correlation_type == 'spearman' and 'spearman' in rep_corr_data:
                        corr_matrix[i, j] = rep_corr_data['spearman'] if rep_corr_data['spearman'] is not None else 1.0
                    else:
                        corr_matrix[i, j] = rep_corr_data['correlation'] if rep_corr_data['correlation'] is not None else 1.0
                else:
                    corr_matrix[i, j] = rep_corr_data if rep_corr_data is not None else 1.0
                    
                rna_same_matrix[i, j] = True
            elif i < j:
                # Upper triangle
                corr_data = report['cross_cell_correlations'].get(cell1, {}).get(cell2)
                if corr_data:
                    if correlation_type == 'spearman' and 'spearman' in corr_data:
                        corr_matrix[i, j] = corr_data['spearman']
                        corr_matrix[j, i] = corr_data['spearman']
                    else:
                        corr_matrix[i, j] = corr_data['correlation']
                        corr_matrix[j, i] = corr_data['correlation']
                    
                    rna_same_matrix[i, j] = corr_data['same_rna_type']
                    rna_same_matrix[j, i] = corr_data['same_rna_type']
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    
    # First heatmap: All correlations
    im1 = ax1.imshow(corr_matrix, cmap='viridis', vmin=-1, vmax=1)
    
    # Add correlation values
    for i in range(matrix_size):
        for j in range(matrix_size):
            if not np.isnan(corr_matrix[i, j]):
                text_color = 'white' if abs(corr_matrix[i, j]) < 0.5 else 'black'
                ax1.text(j, i, f'{corr_matrix[i, j]:.2f}', 
                       ha='center', va='center', color=text_color)
    
    # Labels and title
    ax1.set_xticks(np.arange(matrix_size))
    ax1.set_yticks(np.arange(matrix_size))
    ax1.set_xticklabels(cell_lines, rotation=45, ha='right')
    ax1.set_yticklabels(cell_lines)
    
    # Add transformation info to title
    transform_text = "Log-transformed " if use_log else ""
    ax1.set_title(f'Cell Line Expression {transform_text}{correlation_type.capitalize()} Correlation Matrix (All)')
    
    # Add RNA type annotation
    for i, rna_type in enumerate(rna_types):
        ax1.text(-0.5, i, rna_type, ha='right', va='center', fontweight='bold')
    
    # Second heatmap: Only same RNA type
    masked_matrix = np.copy(corr_matrix)
    masked_matrix[~rna_same_matrix] = np.nan
    
    im2 = ax2.imshow(masked_matrix, cmap='viridis', vmin=-1, vmax=1)
    
    # Add correlation values
    for i in range(matrix_size):
        for j in range(matrix_size):
            if not np.isnan(masked_matrix[i, j]):
                text_color = 'white' if abs(masked_matrix[i, j]) < 0.5 else 'black'
                ax2.text(j, i, f'{masked_matrix[i, j]:.2f}', 
                       ha='center', va='center', color=text_color)
    
    # Labels and title
    ax2.set_xticks(np.arange(matrix_size))
    ax2.set_yticks(np.arange(matrix_size))
    ax2.set_xticklabels(cell_lines, rotation=45, ha='right')
    ax2.set_yticklabels(cell_lines)
    ax2.set_title(f'Cell Line Expression {transform_text}{correlation_type.capitalize()} Correlation Matrix (Same RNA Type Only)')
    
    # Add RNA type annotation
    for i, rna_type in enumerate(rna_types):
        ax2.text(-0.5, i, rna_type, ha='right', va='center', fontweight='bold')
    
    # Add colorbar
    fig.colorbar(im1, ax=ax1, label=f'{correlation_type.capitalize()} Correlation')
    fig.colorbar(im2, ax=ax2, label=f'{correlation_type.capitalize()} Correlation')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to {output_file}")
    
def generate_gene_barplot(report, output_file, use_log=False):
    """Generate bar plot of gene expression across cell lines grouped by RNA type"""
    # Prepare data
    genes = list(report['gene_expression_summary'].keys())
    
    # Check if we have any genes to plot
    if not genes:
        print(f"Warning: No genes found in the report for the bar plot")
        return
    
    # Get cell lines grouped by RNA type
    rna_types = report['rna_type_groups']
    
    # Create figure with subplots
    n_genes = len(genes)
    fig, axes = plt.subplots(n_genes, 1, figsize=(12, 4*n_genes))
    if n_genes == 1:
        axes = [axes]
    
    # Set colors for RNA types
    rna_colors = {
        'total': '#3498db',  # Blue
        'polyA+': '#e74c3c',  # Red
        'polyA-': '#2ecc71'   # Green
    }
    
    # Generate bar plot for each gene
    for i, gene in enumerate(genes):
        gene_data = report['gene_expression_summary'][gene]
        
        # Group data by RNA type
        grouped_data = []
        group_names = []
        colors = []
        
        for rna_type in ['total', 'polyA+', 'polyA-']:
            for cell in rna_types.get(rna_type, []):
                if cell in gene_data['expression_by_cell']:
                    value = gene_data['expression_by_cell'][cell]['mean']
                    grouped_data.append(value)
                    group_names.append(cell)
                    colors.append(rna_colors[rna_type])
        
        if use_log:
            # Add small value to avoid log(0)
            grouped_data = [np.log2(v + 0.1) for v in grouped_data]
            ylabel = 'log2(TPM + 0.1)'
        else:
            ylabel = 'TPM'  
                  
        # Create bars
        x = np.arange(len(group_names))
        bars = axes[i].bar(x, grouped_data, color=colors)
        axes[i].set_xticks(x)
        axes[i].set_xticklabels(group_names, rotation=45, ha='right')
        axes[i].set_title(f'{gene} Expression')
        axes[i].set_ylabel(ylabel)
        
        # Add consistency info
        if gene_data.get('consistency'):
            axes[i].text(0.02, 0.92, 
                       f"Consistency: {gene_data['consistency'].upper()}",
                       transform=axes[i].transAxes, fontweight='bold')
        
        # Add RNA type legend
        handles = [plt.Rectangle((0,0),1,1, color=color) for color in rna_colors.values()]
        labels = rna_colors.keys()
        axes[i].legend(handles, labels, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Gene expression plot saved to {output_file}")

def generate_cluster_heatmap(cluster_report, output_file):
    """Generate heatmap showing expression clusters"""
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
    
    # Create heatmap
    plt.figure(figsize=(12, 10))
    ax = plt.gca()
    
    # Plot heatmap
    im = ax.imshow(reordered_matrix, cmap='viridis', vmin=-1, vmax=1)
    
    # Add correlation values
    for i in range(len(ordered_cells)):
        for j in range(len(ordered_cells)):
            text_color = 'white' if abs(reordered_matrix[i, j]) < 0.5 else 'black'
            ax.text(j, i, f'{reordered_matrix[i, j]:.2f}', 
                   ha='center', va='center', color=text_color)
    
    # Add cluster boundaries
    current_cluster = ordered_clusters[0]
    boundaries = [0]
    
    for i, cluster in enumerate(ordered_clusters[1:], 1):
        if cluster != current_cluster:
            boundaries.append(i - 0.5)
            current_cluster = cluster
    
    # Draw lines for cluster boundaries
    for b in boundaries[1:]:
        ax.axhline(y=b, color='black', linestyle='-', linewidth=2)
        ax.axvline(x=b, color='black', linestyle='-', linewidth=2)
    
    # Labels and title
    ax.set_xticks(np.arange(len(ordered_cells)))
    ax.set_yticks(np.arange(len(ordered_cells)))
    ax.set_xticklabels(ordered_cells, rotation=45, ha='right')
    ax.set_yticklabels(ordered_cells)
    
    # Add cluster labels on right side
    for i, cluster in enumerate(ordered_clusters):
        ax.text(len(ordered_cells) + 0.5, i, cluster, 
               ha='left', va='center', fontweight='bold')
    
    plt.title('Cell Line Expression Clusters')
    plt.colorbar(im, label='Correlation')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Cluster heatmap saved to {output_file}")

import datetime
import shutil

# Modify your main function to include dated output directories
def main():
    parser = argparse.ArgumentParser(description='Analyze RNA-seq gene quantification TSV files across cell lines')
    parser.add_argument('--cell-lines', nargs='+', 
                      help='Specific cell lines to analyze (default: all)')
    parser.add_argument('--plots', action='store_true',
                      help='Generate visualization plots')
    parser.add_argument('--cluster-threshold', type=float, default=0.9,
                      help='Correlation threshold for clustering (default: 0.9)')
    parser.add_argument('--use_log', action='store_true',
                    help='Use log transformation for correlation calculations')
    parser.add_argument('--correlation_type', choices=['pearson', 'spearman'], default='pearson',
                        help='Type of correlation to use in visualizations (default: pearson)')
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = TSVAnalyzer()
    
    print("\n===== Gene Quantification TSV Analysis =====")
    print("This analysis uses gene-level TPM values directly from ENCODE's")
    print("quantification files, which is more suitable for WGS integration")
    print("than analyzing raw bigWig signal tracks.\n")
    
    try:
        # Create a dated folder for all outputs
        today = datetime.datetime.now().strftime("%Y%m%d")
        dated_analysis_dir = analyzer.analysis_dir / today
        dated_analysis_dir.mkdir(exist_ok=True)
        
        
        # Run analysis
        if args.cell_lines:
            all_results = {}
            for cell_line in args.cell_lines:
                if cell_line in analyzer.cell_lines:
                    results = analyzer.analyze_cell_line(cell_line)
                    if results and results['replicates']:
                        all_results[cell_line] = results
                else:
                    print(f"Warning: Unknown cell line '{cell_line}'")
        else:
            all_results = analyzer.analyze_all_cell_lines()
        
        # Generate report and save to dated directory
        report = analyzer.generate_report(all_results, use_log=args.use_log, correlation_type=args.correlation_type)
        dated_report_file = dated_analysis_dir / "tsv_analysis_report.json"
        with open(dated_report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Also save a copy to the main analysis directory for latest version
        main_report_file = analyzer.analysis_dir / "tsv_analysis_report.json"
        shutil.copy(dated_report_file, main_report_file)
        
        # Print summary
        print_report_summary(report)
        
        # Generate expression clusters
        cluster_report = analyzer.identify_expression_clusters(
            all_results, min_correlation=args.cluster_threshold)
        
        # Print cluster summary
        print("\nExpression Clusters:")
        for cluster in cluster_report['clusters']:
            print(f"  {cluster['id']}: {', '.join(cluster['cells'])}")
        
        if cluster_report['outliers']:
            print("\nOutliers:")
            for outlier in cluster_report['outliers']:
                print(f"  {outlier['cell']} ({outlier['rna_type']})")
        
        # Save cluster report to dated directory
        dated_cluster_file = dated_analysis_dir / "tsv_expression_clusters.json"
        with open(dated_cluster_file, 'w') as f:
            json.dump(cluster_report, f, indent=2)
            
        # Also save a copy to the main analysis directory
        main_cluster_file = analyzer.analysis_dir / "tsv_expression_clusters.json"
        shutil.copy(dated_cluster_file, main_cluster_file)
        
        # Generate plots if requested
        if args.plots:
            # Create plots directory within the dated directory
            dated_plots_dir = dated_analysis_dir / 'plots'
            dated_plots_dir.mkdir(exist_ok=True)
            
            # Also create plots directory in main analysis directory if it doesn't exist
            main_plots_dir = analyzer.analysis_dir / 'plots'
            main_plots_dir.mkdir(exist_ok=True)
            
            # Generate correlation heatmap
            dated_heatmap_file = dated_plots_dir / 'tsv_correlation_heatmap.png'
            generate_heatmap(report, dated_heatmap_file, 
                correlation_type=args.correlation_type, 
                use_log=args.use_log)
            
            
            main_heatmap_file = main_plots_dir / 'tsv_correlation_heatmap.png'
            shutil.copy(dated_heatmap_file, main_heatmap_file)

            # Generate gene expression bars
            dated_barplot_file = dated_plots_dir / 'tsv_gene_expression.png'
            generate_gene_barplot(report, dated_barplot_file, use_log=args.use_log)
            
            main_barplot_file = main_plots_dir / 'tsv_gene_expression.png'
            shutil.copy(dated_barplot_file, main_barplot_file)
            
            
            # Generate cluster heatmap
            dated_cluster_heatmap_file = dated_plots_dir / 'tsv_expression_clusters.png'
            cluster_report = analyzer.identify_expression_clusters(
                all_results, min_correlation=args.cluster_threshold, 
                use_log=args.use_log, correlation_type=args.correlation_type)
            generate_cluster_heatmap(cluster_report, dated_cluster_heatmap_file)
            
            main_cluster_heatmap_file = main_plots_dir / 'tsv_expression_clusters.png'
            shutil.copy(dated_cluster_heatmap_file, main_cluster_heatmap_file)
        
        print("\n=== Analysis Complete ===")
        print(f"Reports saved to:")
        print(f"  - Latest version: {analyzer.analysis_dir}/")
        print(f"  - Dated archive: {dated_analysis_dir}/")
        
        if args.plots:
            print(f"Plots saved to:")
            print(f"  - Latest version: {analyzer.analysis_dir}/plots/")
            print(f"  - Dated archive: {dated_plots_dir}/")
        
    except Exception as e:
        print(f"\n!!! Error during analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        
if __name__ == "__main__":
    main()
