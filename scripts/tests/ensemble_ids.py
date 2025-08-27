import pandas as pd

# Load the Entrez to Ensembl mapping
mapping_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/entrez_to_ensembl_mapping.csv'
entrez_mapping = pd.read_csv(mapping_file, low_memory=False)

# Check the structure
print(f"Columns in mapping file: {entrez_mapping.columns.tolist()}")
print(f"Number of entries: {len(entrez_mapping)}")

# Check a few examples
print("\nSample mappings:")
print(entrez_mapping.head())

# Check distribution of Entrez IDs
print(f"\nNumber of unique Entrez IDs: {entrez_mapping['entrez_id'].nunique()}")