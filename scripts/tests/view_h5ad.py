import scanpy as sc
import pandas as pd

# Load the combined dataset
adata = sc.read_h5ad('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_all_genes_standardized.h5ad')

# Basic information about the dataset
print(f"Dataset shape: {adata.shape}")
print(f"Number of samples: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")

# Look at the metadata fields
print("\nMetadata fields (adata.uns):")
for key in adata.uns.keys():
    print(f"  - {key}")

# Check observation (sample) annotations
print("\nSample annotation fields (adata.obs.columns):")
for col in adata.obs.columns:
    print(f"  - {col}")

# Show the first few samples
print("\nFirst 5 samples:")
print(adata.obs.head())

# Check if cell type information is included
has_cell_type = 'cell_type' in adata.obs.columns
has_cell_ont = 'cell_type_ontology' in adata.obs.columns
has_anat_ont = 'anatomical_entity_ontology' in adata.obs.columns

print(f"\nCell type information:")
print(f"  - cell_type column: {'Present' if has_cell_type else 'Missing'}")
print(f"  - cell_type_ontology column: {'Present' if has_cell_ont else 'Missing'}")
print(f"  - anatomical_entity_ontology column: {'Present' if has_anat_ont else 'Missing'}")

# Check ontology mapping coverage
if 'tissue_ontology' in adata.obs.columns:
    missing = adata.obs['tissue_ontology'].isna().sum()
    perc = (adata.n_obs - missing) / adata.n_obs * 100
    print(f"\nTissue ontology mapping: {perc:.1f}% complete ({adata.n_obs - missing}/{adata.n_obs})")

if has_cell_ont:
    missing = adata.obs['cell_type_ontology'].isna().sum()
    perc = (adata.n_obs - missing) / adata.n_obs * 100
    print(f"Cell type ontology mapping: {perc:.1f}% complete ({adata.n_obs - missing}/{adata.n_obs})")

# Show unique dataset sources
if 'dataset' in adata.obs.columns:
    print("\nDataset sources:")
    for ds, count in adata.obs['dataset'].value_counts().items():
        print(f"  - {ds}: {count} samples")