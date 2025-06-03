
# ADNI Gene Expression Clustering Analysis Report

## Dataset Overview
- **Samples**: 650
- **Genes**: 2,000 (highly variable genes used for clustering)
- **Date**: 2025-06-03 17:28:38

## Clustering Results

### Best Clustering Method: kmeans_k2
- **Silhouette Score**: 0.199
- **Number of Clusters**: 2

### All Methods Comparison:
- **kmeans_k2**: Silhouette = 0.199
- **kmeans_k3**: Silhouette = 0.124
- **kmeans_k4**: Silhouette = 0.114
- **kmeans_k5**: Silhouette = 0.101
- **hierarchical_k2**: Silhouette = 0.157
- **hierarchical_k3**: Silhouette = 0.091
- **hierarchical_k4**: Silhouette = 0.086
- **hierarchical_k5**: Silhouette = 0.094

## Cluster-Diagnosis Relationship
- **Adjusted Rand Index**: -0.003
  - 0.0 = random clustering relative to diagnosis
  - 1.0 = perfect agreement with diagnosis
  - Values > 0.3 suggest meaningful relationship

### Diagnosis Distribution by Cluster:

**Cluster 0:**
- Alzheimer's Disease: 33.4%
- Cognitively Normal: 25.1%
- Mild Cognitive Impairment: 41.5%

**Cluster 1:**
- Alzheimer's Disease: 32.2%
- Cognitively Normal: 22.0%
- Mild Cognitive Impairment: 45.9%

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
