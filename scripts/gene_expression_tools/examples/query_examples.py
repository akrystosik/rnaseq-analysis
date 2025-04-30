#!/usr/bin/env python3
"""
Gene Expression Query Examples

This script demonstrates how to use the optimized gene expression functions
for various common use cases.
"""

import sys
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Add parent directory to path when running directly
if __name__ == "__main__":
    # We need to add the parent directory of gene_expression_tools to the path
    script_dir = os.path.dirname(os.path.abspath(__file__))  # examples directory
    package_dir = os.path.dirname(script_dir)  # gene_expression_tools
    scripts_dir = os.path.dirname(package_dir)  # scripts directory
    sys.path.append(scripts_dir)

# Import from the main package
from gene_expression_tools.optimized_gene_expression import (
    load_expression_data,
    get_gene_expression,
    get_gene_expression_matrix,
    get_available_tissues,
    get_available_donors,
    get_all_tissues_gene_expression
)

def example_1_basic_query(data_loader=None):
    """Basic gene expression query example."""
    print("\n=== Example 1: Basic Gene Expression Query ===")
    
    # Use provided data loader or initialize one
    if data_loader is None:
        # Initialize the data loader (only need to do this once)
        data_loader = load_expression_data()
    
    # Query TP53 expression in ENTEx dataset
    start_time = time.time()
    expr = get_gene_expression("entex", "ENSG00000141510", data_loader=data_loader)  # TP53 Ensembl ID
    query_time = time.time() - start_time
    
    print(f"TP53 expression in ENTEx dataset:")
    # Handle the case where 'mean_expression' might be 'N/A'
    mean_expr = expr.get('mean_expression', 'N/A')
    if mean_expr != 'N/A':
        print(f"  Mean TPM: {mean_expr:.4f}")
    else:
        print(f"  Mean TPM: {mean_expr}")
    
    # Handle the case where 'median_expression' might be 'N/A'
    median_expr = expr.get('median_expression', 'N/A')
    if median_expr != 'N/A':
        print(f"  Median TPM: {median_expr:.4f}")
    else:
        print(f"  Median TPM: {median_expr}")
        
    print(f"  Sample count: {expr.get('sample_count', 'N/A')}")
    print(f"  Query time: {query_time:.4f} seconds")
    
    # Query with tissue filter
    start_time = time.time()
    expr = get_gene_expression("entex", "ENSG00000141510", tissue="lung", data_loader=data_loader)
    query_time = time.time() - start_time
    
    print(f"\nTP53 expression in lung tissue (ENTEx):")
    # Handle the case where 'mean_expression' might be 'N/A'
    mean_expr = expr.get('mean_expression', 'N/A')
    if mean_expr != 'N/A':
        print(f"  Mean TPM: {mean_expr:.4f}")
    else:
        print(f"  Mean TPM: {mean_expr}")
    
    # Handle the case where 'median_expression' might be 'N/A'
    median_expr = expr.get('median_expression', 'N/A')
    if median_expr != 'N/A':
        print(f"  Median TPM: {median_expr:.4f}")
    else:
        print(f"  Median TPM: {median_expr}")
        
    print(f"  Sample count: {expr.get('sample_count', 'N/A')}")
    print(f"  Query time: {query_time:.4f} seconds")
    
    # Query with donor filter
    donors = get_available_donors("entex", data_loader=data_loader)
    if donors:
        start_time = time.time()
        expr = get_gene_expression("entex", "ENSG00000141510", donor=donors[0], data_loader=data_loader)
        query_time = time.time() - start_time
        
        print(f"\nTP53 expression for donor {donors[0]} (ENTEx):")
        # Handle the case where 'mean_expression' might be 'N/A'
        mean_expr = expr.get('mean_expression', 'N/A')
        if mean_expr != 'N/A':
            print(f"  Mean TPM: {mean_expr:.4f}")
        else:
            print(f"  Mean TPM: {mean_expr}")
        
        # Handle the case where 'median_expression' might be 'N/A'
        median_expr = expr.get('median_expression', 'N/A')
        if median_expr != 'N/A':
            print(f"  Median TPM: {median_expr:.4f}")
        else:
            print(f"  Median TPM: {median_expr}")
            
        print(f"  Sample count: {expr.get('sample_count', 'N/A')}")
        print(f"  Query time: {query_time:.4f} seconds")
    
    return data_loader

def example_2_expression_matrix(data_loader=None):
    """Expression matrix query example."""
    print("\n=== Example 2: Expression Matrix Query ===")
    
    # Use provided data loader or initialize one
    if data_loader is None:
        data_loader = load_expression_data()
    
    # Query for multiple genes (using Ensembl IDs for reliability)
    genes = ("ENSG00000141510", "ENSG00000012048", "ENSG00000139618")  # TP53, BRCA1, BRCA2
    
    # First query (cold cache)
    start_time = time.time()
    matrix, gene_ids, sample_ids = get_gene_expression_matrix("entex", genes, data_loader=data_loader)
    query_time = time.time() - start_time
    
    print(f"Expression matrix for {', '.join(genes)}:")
    print(f"  Matrix shape: {matrix.shape}")
    print(f"  Genes found: {', '.join(gene_ids)}")
    print(f"  Number of samples: {len(sample_ids)}")
    print(f"  Query time (cold cache): {query_time:.4f} seconds")
    
    # Second query (warm cache)
    start_time = time.time()
    matrix, gene_ids, sample_ids = get_gene_expression_matrix("entex", genes, data_loader=data_loader)
    query_time = time.time() - start_time
    
    print(f"\nQuery time (warm cache): {query_time:.4f} seconds")
    
    # Query with tissue filter
    tissues = get_available_tissues("entex", data_loader=data_loader)
    if tissues and len(tissues) >= 2:
        tissue_tuple = (tissues[0], tissues[1])
        
        start_time = time.time()
        matrix, gene_ids, sample_ids = get_gene_expression_matrix("entex", genes, tissues=tissue_tuple, data_loader=data_loader)
        query_time = time.time() - start_time
        
        print(f"\nExpression matrix for {', '.join(genes)} in {', '.join(tissue_tuple)}:")
        print(f"  Matrix shape: {matrix.shape}")
        print(f"  Number of samples: {len(sample_ids)}")
        print(f"  Query time: {query_time:.4f} seconds")
        
        # Calculate mean expression per gene
        if matrix.size > 0:
            means = np.mean(matrix, axis=0)
            for i, gene_id in enumerate(gene_ids):
                print(f"  Mean {gene_id} expression: {means[i]:.4f}")
    
    return data_loader

def example_3_cross_tissue_analysis(data_loader=None):
    """Cross-tissue gene expression analysis."""
    print("\n=== Example 3: Cross-Tissue Gene Expression Analysis ===")
    
    # Use provided data loader or initialize one
    if data_loader is None:
        data_loader = load_expression_data()
    
    # Get tissues
    tissues = get_available_tissues("entex", data_loader=data_loader)
    if not tissues:
        print("No tissues found in ENTEx dataset")
        return data_loader
    
    # Query gene across all tissues
    gene = "ENSG00000141510"  # TP53 Ensembl ID
    results = []
    
    print(f"Querying TP53 expression across {len(tissues)} tissues...")
    start_time = time.time()
    
    for tissue in tissues:
        expr = get_gene_expression("entex", gene, tissue=tissue, data_loader=data_loader)
        if expr and expr.get('sample_count', 0) > 0:
            results.append({
                'tissue': tissue,
                'mean_expression': expr['mean_expression'],
                'median_expression': expr['median_expression'],
                'sample_count': expr['sample_count']
            })
    
    query_time = time.time() - start_time
    print(f"Query time for all tissues: {query_time:.4f} seconds")
    
    # Alternative using bulk query function
    start_time = time.time()
    all_tissues = get_all_tissues_gene_expression(gene, ["entex"], data_loader=data_loader)
    bulk_query_time = time.time() - start_time
    
    print(f"Bulk query time: {bulk_query_time:.4f} seconds")
    
    # Sort by expression level
    if results:
        results.sort(key=lambda x: x['mean_expression'], reverse=True)
        
        print(f"\nTop 5 tissues by TP53 expression:")
        for result in results[:5]:
            print(f"  {result['tissue']}: {result['mean_expression']:.4f} TPM ({result['sample_count']} samples)")
        
        print(f"\nBottom 5 tissues by TP53 expression:")
        for result in results[-5:]:
            print(f"  {result['tissue']}: {result['mean_expression']:.4f} TPM ({result['sample_count']} samples)")
    
    return data_loader

def example_4_performance_benchmark(data_loader=None):
    """Performance benchmark for repeated queries."""
    print("\n=== Example 4: Performance Benchmark ===")
    
    # Use provided data loader or initialize one
    if data_loader is None:
        data_loader = load_expression_data()
    
    # Get all tissues and a list of genes
    tissues = get_available_tissues("entex", data_loader=data_loader)
    if not tissues:
        print("No tissues found in ENTEx dataset")
        return data_loader
    
    # List of genes to query (using Ensembl IDs)
    genes = [
        "ENSG00000141510",  # TP53
        "ENSG00000012048",  # BRCA1
        "ENSG00000139618",  # BRCA2
        "ENSG00000146648",  # EGFR
        "ENSG00000136997",  # MYC
        "ENSG00000171862",  # PTEN
        "ENSG00000139687",  # RB1
        "ENSG00000133703",  # KRAS
        "ENSG00000174775",  # HRAS
        "ENSG00000213281"   # NRAS
    ]
    
    # Reduce number of queries for testing
    num_queries = 20  # Reduced from 100 to make it faster for testing
    print(f"Running {num_queries} individual queries...")
    
    start_time = time.time()
    for _ in range(num_queries // len(genes)):
        for gene in genes:
            # Randomly select a tissue
            tissue = tissues[np.random.randint(0, len(tissues))]
            _ = get_gene_expression("entex", gene, tissue=tissue, data_loader=data_loader)
    
    individual_time = time.time() - start_time
    print(f"Individual queries time: {individual_time:.4f} seconds")
    print(f"Average time per query: {individual_time / num_queries:.6f} seconds")
    
    # Benchmark: Matrix queries
    num_matrix_queries = 2  # Reduced from 10 to make it faster for testing
    print(f"\nRunning {num_matrix_queries} matrix queries...")
    
    start_time = time.time()
    for _ in range(num_matrix_queries):
        # Randomly select 3 tissues
        tissue_indices = np.random.choice(len(tissues), 3, replace=False)
        tissue_tuple = tuple(tissues[i] for i in tissue_indices)
        
        # Query all genes at once
        _, _, _ = get_gene_expression_matrix("entex", tuple(genes), tissues=tissue_tuple, data_loader=data_loader)
    
    matrix_time = time.time() - start_time
    print(f"Matrix queries time: {matrix_time:.4f} seconds")
    print(f"Average time per matrix query: {matrix_time / num_matrix_queries:.6f} seconds")
    
    return data_loader

def example_5_visualization(data_loader=None):
    """Visualization example."""
    print("\n=== Example 5: Visualization ===")
    
    # Use provided data loader or initialize one
    if data_loader is None:
        data_loader = load_expression_data()
    
    # Get gene expression across tissues
    gene = "ENSG00000141510"  # TP53 Ensembl ID
    all_tissues = get_all_tissues_gene_expression(gene, ["entex"], data_loader=data_loader)
    
    if not all_tissues or "entex" not in all_tissues:
        print(f"No expression data found for TP53")
        return data_loader
    
    # Extract tissue expression data
    tissues = []
    expressions = []
    
    for tissue, expr in all_tissues["entex"].items():
        tissues.append(tissue)
        expressions.append(expr['mean_expression'])
    
    # Sort by expression level
    sorted_indices = np.argsort(expressions)[::-1]
    tissues = [tissues[i] for i in sorted_indices]
    expressions = [expressions[i] for i in sorted_indices]
    
    # Create bar chart
    plt.figure(figsize=(12, 8))
    plt.bar(range(len(tissues)), expressions)
    plt.xticks(range(len(tissues)), tissues, rotation=90)
    plt.xlabel('Tissue')
    plt.ylabel('Mean Expression (TPM)')
    plt.title('TP53 Expression Across Tissues (ENTEx)')
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(os.getcwd(), 'TP53_tissue_expression.png')
    plt.savefig(output_file)
    print(f"Saved visualization to {output_file}")
    
    # Print top 5 tissues
    print(f"\nTop 5 tissues by TP53 expression:")
    for i in range(min(5, len(tissues))):
        print(f"  {tissues[i]}: {expressions[i]:.4f} TPM")
    
    return data_loader

def example_6_direct_file_loading():
    """Example showing how to load data directly from a file."""
    print("\n=== Example 6: Direct File Loading ===")
    
    # Define specific files to load
    combined_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_all_genes_standardized.h5ad"
    
    print(f"Loading data directly from file: {combined_file}")
    
    # Load the file directly
    start_time = time.time()
    loader = load_expression_data(specific_files={"combined_direct": combined_file})
    load_time = time.time() - start_time
    
    print(f"Data loaded in {load_time:.2f} seconds")
    
    # Query gene expression
    gene = "ENSG00000141510"  # TP53 Ensembl ID
    expr = get_gene_expression("combined_direct", gene, data_loader=loader)
    
    print(f"\nTP53 expression in combined dataset:")
    # Handle the case where 'mean_expression' might be 'N/A'
    mean_expr = expr.get('mean_expression', 'N/A')
    if mean_expr != 'N/A':
        print(f"  Mean TPM: {mean_expr:.4f}")
    else:
        print(f"  Mean TPM: {mean_expr}")
    
    # Handle the case where 'median_expression' might be 'N/A'
    median_expr = expr.get('median_expression', 'N/A')
    if median_expr != 'N/A':
        print(f"  Median TPM: {median_expr:.4f}")
    else:
        print(f"  Median TPM: {median_expr}")
        
    print(f"  Sample count: {expr.get('sample_count', 'N/A')}")
    
    # Check dataset distribution if available
    if hasattr(loader.datasets["combined_direct"].obs, 'columns') and 'dataset' in loader.datasets["combined_direct"].obs.columns:
        dataset_counts = loader.datasets["combined_direct"].obs['dataset'].value_counts()
        print("\nDataset distribution in the combined file:")
        for dataset, count in dataset_counts.items():
            print(f"  {dataset}: {count} samples")
    
    # Get available tissues
    tissues = get_available_tissues("combined_direct", data_loader=loader)
    print(f"\nFound {len(tissues)} tissues")
    if tissues:
        print(f"First 5 tissues: {', '.join(tissues[:5])}")
    
    # Get available donors
    donors = get_available_donors("combined_direct", data_loader=loader)
    print(f"Found {len(donors)} donors")
    if donors:
        print(f"First 5 donors: {', '.join(donors[:5])}")
    
    return loader

def main():
    """Run all examples."""
    print("=== Gene Expression Query Examples ===")
    
    # Example with direct file loading
    print("\nRunning example with direct file loading...")
    loader = example_6_direct_file_loading()
    
    # Basic examples using the same loader to avoid reloading data
    loader = example_1_basic_query(loader)
    loader = example_2_expression_matrix(loader)
    loader = example_3_cross_tissue_analysis(loader)
    loader = example_4_performance_benchmark(loader)
    loader = example_5_visualization(loader)
    
    print("\nAll examples completed successfully!")

if __name__ == "__main__":
    main()