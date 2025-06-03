import pandas as pd

# Load the collected gene IDs
gene_ids_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/encode_original_gene_ids.csv'
gene_ids = pd.read_csv(gene_ids_file)

# Get sample of each ID type
ensembl_ids = gene_ids[gene_ids['is_ensembl']]['gene_id'].head(10).tolist()
numeric_ids = gene_ids[gene_ids['is_numeric']]['gene_id'].head(10).tolist()
other_ids = gene_ids[~gene_ids['is_ensembl'] & ~gene_ids['is_numeric']]['gene_id'].head(20).tolist()

print("Sample Ensembl IDs:")
for id in ensembl_ids:
    print(id)

print("\nSample numeric IDs:")
for id in numeric_ids:
    print(id)

print("\nSample 'other' IDs:")
for id in other_ids:
    print(id)

# Analyze patterns in the "other" IDs
other_patterns = []
for id in gene_ids[~gene_ids['is_ensembl'] & ~gene_ids['is_numeric']]['gene_id']:
    id_str = str(id)
    if ":" in id_str:
        pattern = id_str.split(":")[0] + ":"
        other_patterns.append(pattern)
    elif "_" in id_str:
        pattern = id_str.split("_")[0] + "_"
        other_patterns.append(pattern)
    else:
        first_chars = id_str[:3] if len(id_str) >= 3 else id_str
        other_patterns.append(first_chars)

# Count the most common patterns
from collections import Counter
pattern_counts = Counter(other_patterns).most_common(10)

print("\nMost common patterns in 'other' IDs:")
for pattern, count in pattern_counts:
    print(f"{pattern}: {count} ({count/len(other_patterns)*100:.2f}%)")