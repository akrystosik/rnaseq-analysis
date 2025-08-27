
#!/usr/bin/env python3
"""
Simplified visualization generator for RNA-seq analysis
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

# Create a simple expression heatmap
print("Creating expression heatmap...")
plt.figure(figsize=(12, 16))

# Prepare data for heatmap
heatmap_data = df[['ACTB', 'GAPDH', 'TUBB', 'B2M', 'PPIA', 'TBP', 'HPRT1', 'RPL13A']]
log_data = np.log2(heatmap_data + 0.1)  # Log transform with small offset

# Add annotations for dataset types
is_tissue = df['is_tissue'].values
rna_types = df['rna_type'].values

# Create heatmap
ax = sns.heatmap(log_data, cmap='viridis', 
               yticklabels=df['dataset_id'].values)

# Add custom y-tick labels that include RNA type
ytick_labels = []
for i, dataset in enumerate(df['dataset_id']):
    prefix = "T: " if is_tissue[i] else "C: "
    ytick_labels.append(f"{prefix}{dataset} ({rna_types[i]})")

ax.set_yticklabels(ytick_labels)

plt.title('Log2 Gene Expression Across Datasets')
plt.tight_layout()
plt.savefig(f"{output_dir}/gene_expression_heatmap.png", dpi=300)
print(f"Saved heatmap to {output_dir}/gene_expression_heatmap.png")

# Create a correlation heatmap
print("Generating correlation heatmap...")
# Convert the expression data to a correlation matrix
expr_corr = log_data.T.corr()

plt.figure(figsize=(14, 12))
mask = np.triu(np.ones_like(expr_corr, dtype=bool))
cmap = sns.diverging_palette(220, 10, as_cmap=True)

# Plot correlation heatmap
sns.heatmap(expr_corr, mask=mask, cmap=cmap, vmin=-1, vmax=1,
           annot=True, fmt=".2f", square=True, linewidths=0.5,
           cbar_kws={"shrink": .5})

# Create custom tick labels
tick_labels = []
for i, dataset in enumerate(df['dataset_id']):
    prefix = "T: " if is_tissue[i] else "C: "
    tick_labels.append(f"{prefix}{dataset}")

plt.xticks(np.arange(len(tick_labels)) + 0.5, tick_labels, rotation=90)
plt.yticks(np.arange(len(tick_labels)) + 0.5, tick_labels)
plt.title('Gene Expression Correlation Between Datasets')
plt.tight_layout()
plt.savefig(f"{output_dir}/correlation_heatmap.png", dpi=300)
print(f"Saved correlation heatmap to {output_dir}/correlation_heatmap.png")

# Create a separate cell lines and tissues visualizations
cell_lines_df = df[~df['is_tissue']]
tissues_df = df[df['is_tissue']]

# Cell lines heatmap
if len(cell_lines_df) > 0:
    print("Generating cell lines heatmap...")
    plt.figure(figsize=(10, 8))
    cell_data = cell_lines_df[['ACTB', 'GAPDH', 'TUBB', 'B2M', 'PPIA', 'TBP', 'HPRT1', 'RPL13A']]
    log_cell_data = np.log2(cell_data + 0.1)
    
    ax = sns.clustermap(log_cell_data, cmap='viridis', 
                     figsize=(12, 10),
                     col_cluster=True, row_cluster=True,
                     yticklabels=cell_lines_df['dataset_id'])
    
    plt.title('Cell Lines Gene Expression Clusters')
    plt.savefig(f"{output_dir}/cell_lines_heatmap.png", dpi=300)
    print(f"Saved cell lines heatmap to {output_dir}/cell_lines_heatmap.png")

# Tissues heatmap
if len(tissues_df) > 0:
    print("Generating tissues heatmap...")
    plt.figure(figsize=(10, 8))
    tissue_data = tissues_df[['ACTB', 'GAPDH', 'TUBB', 'B2M', 'PPIA', 'TBP', 'HPRT1', 'RPL13A']]
    log_tissue_data = np.log2(tissue_data + 0.1)
    
    ax = sns.clustermap(log_tissue_data, cmap='viridis', 
                      figsize=(12, 10),
                      col_cluster=True, row_cluster=True,
                      yticklabels=tissues_df['dataset_id'])
    
    plt.title('Tissue Samples Gene Expression Clusters')
    plt.savefig(f"{output_dir}/tissues_heatmap.png", dpi=300)
    print(f"Saved tissues heatmap to {output_dir}/tissues_heatmap.png")

print("Done!")
