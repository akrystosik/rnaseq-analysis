#!/usr/bin/env python3
"""
Basic Quality Control for RNA-seq Gene Quantification Data
---------------------------------------------------------

This module provides basic quality control functions for RNA-seq data:
- Checking housekeeping gene expression
- Replicating correlation assessment
- Identifying low-quality samples
"""

import json
import numpy as np
import pandas as pd
import os
from pathlib import Path

class BasicQC:
    def __init__(self, config_file=None):
        """Initialize the QC module with configuration"""
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
        
        # Define QC thresholds
        self.qc_thresholds = {
            'min_housekeeping_detected': 5,  # Minimum number of housekeeping genes detected
            'min_housekeeping_tpm': 1.0,     # Minimum TPM for a gene to be considered detected
            'min_replicate_correlation': 0.7, # Minimum correlation between replicates
            'min_housekeeping_correlation': 0.5 # Minimum correlation with expected housekeeping pattern
        }
    
    def check_housekeeping_genes(self, expression_data):
        """Check housekeeping gene expression in a sample"""
        detected = 0
        detected_genes = []
        expression_values = []
        
        for gene in self.housekeeping_genes:
            if gene in expression_data and expression_data[gene] >= self.qc_thresholds['min_housekeeping_tpm']:
                detected += 1
                detected_genes.append(gene)
                expression_values.append(expression_data[gene])
        
        # Calculate coefficient of variation for detected genes
        cv = None
        if detected > 0:
            cv = np.std(expression_values) / np.mean(expression_values)
        
        return {
            'detected': detected,
            'total': len(self.housekeeping_genes),
            'percent_detected': detected / len(self.housekeeping_genes) * 100 if self.housekeeping_genes else 0,
            'detected_genes': detected_genes,
            'coefficient_of_variation': cv,
            'passed': detected >= self.qc_thresholds['min_housekeeping_detected']
        }
    
    def check_replicate_correlation(self, replicates):
        """Check correlation between replicates"""
        if len(replicates) < 2:
            return {
                'has_replicates': False,
                'passed': True  # Single replicate passes by default
            }
        
        # Get first two replicates
        rep_keys = list(replicates.keys())
        rep1 = rep_keys[0]
        rep2 = rep_keys[1]
        
        # Get common genes
        common_genes = [gene for gene in self.housekeeping_genes 
                       if gene in replicates[rep1] and gene in replicates[rep2]]
        
        if len(common_genes) < 3:
            return {
                'has_replicates': True,
                'common_genes': len(common_genes),
                'passed': False,
                'reason': "Too few common genes between replicates"
            }
        
        # Calculate correlation
        values1 = [replicates[rep1][gene] for gene in common_genes]
        values2 = [replicates[rep2][gene] for gene in common_genes]
        
        corr = np.corrcoef(values1, values2)[0, 1]
        
        return {
            'has_replicates': True,
            'correlation': corr,
            'common_genes': len(common_genes),
            'total_genes': len(self.housekeeping_genes),
            'passed': corr >= self.qc_thresholds['min_replicate_correlation'],
            'reason': f"Replicate correlation ({corr:.3f}) below threshold" if corr < self.qc_thresholds['min_replicate_correlation'] else None
        }
    
    def assess_sample_quality(self, results):
        """Assess overall quality of a sample"""
        # Check if we have data
        if not results or 'replicates' not in results or not results['replicates']:
            return {
                'passed': False,
                'reason': "No data found"
            }
        
        # Extract average expression across replicates
        avg_expression = {}
        for gene in self.housekeeping_genes:
            values = []
            for rep in results['replicates'].values():
                if gene in rep:
                    values.append(rep[gene])
            
            if values:
                avg_expression[gene] = np.mean(values)
        
        # Check housekeeping genes
        hk_check = self.check_housekeeping_genes(avg_expression)
        
        # Check replicate correlation
        rep_check = self.check_replicate_correlation(results['replicates'])
        
        # Overall assessment
        passed = hk_check['passed'] and rep_check['passed']
        
        # Determine reason for failure
        reason = None
        if not passed:
            if not hk_check['passed']:
                reason = f"Only {hk_check['detected']}/{hk_check['total']} housekeeping genes detected"
            elif not rep_check['passed']:
                reason = rep_check['reason']
        
        return {
            'passed': passed,
            'reason': reason,
            'housekeeping_check': hk_check,
            'replicate_check': rep_check
        }
    
    def run_qc_analysis(self, all_results):
        """Run QC analysis on all results"""
        qc_results = {}
        
        for dataset_id, results in all_results.items():
            qc_results[dataset_id] = self.assess_sample_quality(results)
        
        # Generate summary
        passing = [dataset_id for dataset_id, result in qc_results.items() if result['passed']]
        failing = [dataset_id for dataset_id, result in qc_results.items() if not result['passed']]
        
        summary = {
            'total_datasets': len(qc_results),
            'passing_datasets': len(passing),
            'failing_datasets': len(failing),
            'pass_rate': len(passing) / len(qc_results) * 100 if qc_results else 0,
            'passing_datasets_list': passing,
            'failing_datasets_list': failing,
            'detailed_results': qc_results
        }
        
        # Save report
        qc_report_file = self.analysis_dir / "qc_report.json"
        with open(qc_report_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        return summary
    
    def filter_high_quality_results(self, all_results):
        """Filter results to keep only high-quality datasets"""
        qc_results = self.run_qc_analysis(all_results)
        
        high_quality_results = {}
        for dataset_id in qc_results['passing_datasets_list']:
            high_quality_results[dataset_id] = all_results[dataset_id]
        
        return high_quality_results