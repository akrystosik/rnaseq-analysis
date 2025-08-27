#!/usr/bin/env python3
"""
Generate visualizations directly from the raw results JSON file
"""
import json
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Load raw results
results_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis_new/analysis/20250304/results/raw_results.json"
output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis_new/analysis/20250304/visualizations"
os.makedirs(output_dir, exist_ok=True)

print(f"Loading results from {results_file}")
with open(results_file) as f:
    results = json.load(f)

print(f"Loaded data for {len(results)} datasets")

# Create a simple heatmap of gene expression
print("Generating gene expression heatmap...")
data = []
for dataset_id, dataset_data in results.items():
    # Calculate average expression across replicates
    avg_expr = {}
    for gene in ['ACTB', 'GAPDH', 'TUBB', 'B2M', 'PPIA', 'TBP', 'HPRT1', 'RPL13A']:
        values = []
        for rep, rep_data in dataset_data['replicates'].items():
            if gene in rep_data:
                values.append(rep_data[gene])
        if values:
            avg_expr[gene] = np.mean(values)
        else:
            avg_expr[gene] = 0
    
    # Add metadata
    row = {
        'dataset_id': dataset_id,
        'rna_type': dataset_data['metadata']['rna_type'],
        'is_tissue': 'ENTEx' in dataset_id
    }
    row.update(avg_expr)
    data.append(row)

# Create DataFrame
df = pd.DataFrame(data)
df.set_index('dataset_id', inplace=True)

# Split cell lines and tissues
cell_lines = df[~df['is_tissue']].drop(['is_tissue'], axis=1)
tissues = df[df['is_tissue']].drop(['is_tissue'], axis=1)

# Log transform for better visualization
rna_type = cell_lines.pop('rna_type')
log_data = np.log2(cell_lines + 0.1)
log_data['rna_type'] = rna_type

plt.figure(figsize=(12, 10))
sns.clustermap(log_data.drop(['rna_type'], axis=1), 
              cmap='viridis', 
              z_score=0,
              row_colors=pd.get_dummies(log_data['rna_type']).values,
              yticklabels=True)
plt.title('Gene Expression Across Cell Lines (Log2)')
plt.tight_layout()
plt.savefig(f"{output_dir}/gene_expression_heatmap.png", dpi=300)
print(f"Saved heatmap to {output_dir}/gene_expression_heatmap.png")

# Create tissues heatmap if we have tissue data
if not tissues.empty:
    tissues.pop('rna_type')
    log_tissues = np.log2(tissues + 0.1)
    plt.figure(figsize=(12, 8))
    sns.clustermap(log_tissues, 
                  cmap='viridis', 
                  z_score=0,
                  yticklabels=True)
    plt.title('Gene Expression Across Tissues (Log2)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/tissue_expression_heatmap.png", dpi=300)
    print(f"Saved tissue heatmap to {output_dir}/tissue_expression_heatmap.png")

print("Done!")
