#!/usr/bin/env python3
"""
Main Analysis Script for RNA-seq Gene Quantification Data
--------------------------------------------------------

This script runs the complete analysis pipeline for RNA-seq gene quantification data,
including:
- Dataset management
- Gene identification
- Expression analysis
- Quality control
- Visualization
"""

import argparse
import json
import os
from pathlib import Path
import sys
import pandas as pd
import numpy as np
from datetime import datetime

# Add the project root to the path to import modules
script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
project_root = script_dir.parent
sys.path.append(str(project_root))

from core.dataset_manager import DatasetManager
from core.gene_identifier import GeneIdentifier
from core.expression_analyzer import ExpressionAnalyzer
from qc.basic_qc import BasicQC
from visualization.heatmap import HeatmapVisualizer

def setup_output_directories(base_dir):
    """Set up dated output directories"""
    today = datetime.now().strftime("%Y%m%d")
    analysis_dir = Path(base_dir) / "analysis"
    
    # Create dated directory
    dated_dir = analysis_dir / today
    dated_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    results_dir = dated_dir / "results"
    results_dir.mkdir(exist_ok=True)
    
    viz_dir = dated_dir / "visualizations"
    viz_dir.mkdir(exist_ok=True)
    
    return {
        'base': dated_dir,
        'results': results_dir,
        'visualizations': viz_dir
    }

def analyze_dataset(dataset_id, dataset_manager, gene_identifier):
    """Analyze a single dataset"""
    print(f"\nAnalyzing dataset: {dataset_id}")
    
    # Get dataset info
    dataset_info = dataset_manager.get_cell_line_info(dataset_id)
    if not dataset_info:
        print(f"  Error: Dataset {dataset_id} not found")
        return None
    
    experiment_id = dataset_info['experiment_id']
    print(f"  Experiment ID: {experiment_id}")
    print(f"  RNA type: {dataset_info['rna_type']}")
    
    # Get TSV files
    tsv_files = dataset_manager.get_tsv_files(experiment_id)
    if not tsv_files:
        print(f"  Error: No TSV files found for {dataset_id}")
        return None
    
    print(f"  Found {len(tsv_files)} TSV files")
    
    # Initialize results structure
    results = {
        'metadata': dataset_info,
        'replicates': {}
    }
    
    # Process each TSV file
    for file_info in tsv_files:
        accession = file_info['accession']
        url = file_info['url']
        rep = f"rep{file_info['rep']}"
        
        print(f"  Processing file: {accession}.tsv (Replicate {file_info['rep']})")
        
        # Download TSV if needed
        file_path = dataset_manager.download_tsv(accession, url, dataset_id)
        if not file_path:
            print(f"    Error: Could not download file")
            continue
        
        # Parse TSV
        df = gene_identifier.parse_tsv(file_path)
        if df is None:
            print(f"    Error: Could not parse file")
            continue
        
        print(f"    File parsed successfully ({len(df)} rows)")
        
        # Extract expression for housekeeping genes
        replicate_data = {}
        for gene in gene_identifier.housekeeping_genes:
            tpm = gene_identifier.extract_gene_data(df, gene)
            replicate_data[gene] = tpm
        
        results['replicates'][rep] = replicate_data
        print(f"    Extracted expression for {len(replicate_data)} genes")
    
    # Report success
    if results['replicates']:
        print(f"  Analysis complete: {len(results['replicates'])} replicates processed")
        return results
    else:
        print(f"  Analysis failed: No valid replicates found")
        return None

def main():
    """Run the main analysis"""
    parser = argparse.ArgumentParser(description='Analyze RNA-seq gene quantification data')
    
    # Dataset selection arguments
    parser.add_argument('--datasets', nargs='+', help='Specific datasets to analyze (default: all)')
    parser.add_argument('--include-entex', action='store_true', help='Include ENTEx tissue datasets')
    
    # Analysis options
    parser.add_argument('--cluster-threshold', type=float, default=0.9,
                      help='Correlation threshold for clustering (default: 0.9)')
    parser.add_argument('--use-log', action='store_true',
                      help='Use log transformation for correlation calculations')
    parser.add_argument('--correlation-type', choices=['pearson', 'spearman'], default='pearson',
                      help='Type of correlation to use in visualizations (default: pearson)')
    
    # Output options
    parser.add_argument('--output-dir', help='Custom output directory')
    parser.add_argument('--config', help='Path to custom config file')
    parser.add_argument('--skip-qc', action='store_true', help='Skip quality control filtering')
    parser.add_argument('--skip-plots', action='store_true', help='Skip generating plots')
    
    args = parser.parse_args()
    
    print("\n===== RNA-seq Gene Quantification Analysis =====")
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Set up components
    config_file = args.config
    dataset_manager = DatasetManager(config_file)
    gene_identifier = GeneIdentifier(config_file)
    expression_analyzer = ExpressionAnalyzer(config_file)
    qc_module = BasicQC(config_file)
    visualizer = HeatmapVisualizer(config_file)
    
    # Set up output directories
    if args.output_dir:
        output_dirs = {
            'base': Path(args.output_dir),
            'results': Path(args.output_dir) / "results",
            'visualizations': Path(args.output_dir) / "visualizations"
        }
        for dir_path in output_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
    else:
        output_dirs = setup_output_directories(dataset_manager.base_dir)
    
    print(f"Output directory: {output_dirs['base']}")
    
    # Add ENTEx datasets if requested
    if args.include_entex:
        print("\nAdding ENTEx tissue datasets...")
        added = dataset_manager.add_entex_datasets()
        print(f"  Added {added} ENTEx datasets")
    
    # Select datasets to analyze
    if args.datasets:
        selected_datasets = [d for d in args.datasets if d in dataset_manager.list_cell_lines()]
        print(f"\nSelected {len(selected_datasets)} datasets for analysis")
    else:
        selected_datasets = dataset_manager.list_cell_lines()
        print(f"\nAnalyzing all {len(selected_datasets)} available datasets")
    
    # Run analysis on each dataset
    all_results = {}
    for dataset_id in selected_datasets:
        result = analyze_dataset(dataset_id, dataset_manager, gene_identifier)
        if result and result['replicates']:
            all_results[dataset_id] = result
    
    print(f"\nSuccessfully analyzed {len(all_results)}/{len(selected_datasets)} datasets")
    
    # Save raw results
    raw_results_file = output_dirs['results'] / "raw_results.json"
    with open(raw_results_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"Raw results saved to {raw_results_file}")
    
    # Run quality control
    if not args.skip_qc:
        print("\nRunning quality control...")
        qc_summary = qc_module.run_qc_analysis(all_results)
        print(f"  {qc_summary['passing_datasets']} datasets passed QC")
        print(f"  {qc_summary['failing_datasets']} datasets failed QC")
        
        # Filter to high-quality datasets if needed
        if qc_summary['failing_datasets'] > 0:
            all_results = qc_module.filter_high_quality_results(all_results)
            print(f"  Filtered to {len(all_results)} high-quality datasets")
    
    # Generate expression matrix
    print("\nGenerating expression matrix...")
    expr_matrix = expression_analyzer.generate_expression_matrix(all_results)
    
    # Save expression matrix
    matrix_file = output_dirs['results'] / "expression_matrix.csv"
    expr_matrix.to_csv(matrix_file, index=False)
    print(f"Expression matrix saved to {matrix_file}")
    
    # Run comprehensive analysis
    print("\nRunning comprehensive expression analysis...")
    report = expression_analyzer.generate_comprehensive_report(
        all_results, 
        use_log=args.use_log,
        correlation_type=args.correlation_type
    )
    print("Expression analysis complete")
    
    # Identify expression clusters
    print("\nIdentifying expression clusters...")
    cluster_report = expression_analyzer.identify_expression_clusters(
        all_results,
        min_correlation=args.cluster_threshold,
        use_log=args.use_log,
        correlation_type=args.correlation_type
    )
    
    # Print cluster summary
    print("Expression clusters identified:")
    for cluster in cluster_report['clusters']:
        print(f"  {cluster['id']}: {', '.join(cluster['cells'])}")
    
    if cluster_report['outliers']:
        print("\nOutliers:")
        for outlier in cluster_report['outliers']:
            print(f"  {outlier['cell']} ({outlier.get('rna_type', 'unknown')})")
    
    # Save cluster report
    cluster_file = output_dirs['results'] / "expression_clusters.json"
    with open(cluster_file, 'w') as f:
        json.dump(cluster_report, f, indent=2)
    print(f"Cluster report saved to {cluster_file}")
    
    # Generate visualizations
    if not args.skip_plots:
        print("\nGenerating visualizations...")
        visualizer.generate_all_visualizations(
            report, 
            cluster_report, 
            expr_matrix,
            output_dir=output_dirs['visualizations']
        )
    
    print(f"\n===== Analysis Complete =====")
    print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Results saved to: {output_dirs['base']}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())