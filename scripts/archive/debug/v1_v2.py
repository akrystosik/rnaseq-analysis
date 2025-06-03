import scanpy as sc
import pandas as pd

# --- Load the V1 file ---
adata_v1_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250502_172003/encode_standardized_v1.h5ad'
try:
    adata_v1 = sc.read_h5ad(adata_v1_path)
    print("--- Stage 1 Output (encode_standardized_v1.h5ad) ---")
    print(f"Shape: {adata_v1.shape}")

    # Check var_names (index)
    print("\nSample var_names (index):")
    print(adata_v1.var_names[:10].tolist())
    # Check if index is numeric - important!
    try:
        _ = pd.to_numeric(adata_v1.var_names)
        is_numeric = True
    except ValueError:
        is_numeric = False
    print(f"Is index potentially numeric? {is_numeric}") 

    # Check var dataframe columns
    print("\nVar columns:")
    print(adata_v1.var.columns)

    # --- Check specific columns Stage 1 should have created ---
    if 'gene_id' in adata_v1.var.columns:
        print("\nSample 'gene_id' column values:")
        print(adata_v1.var['gene_id'].head())
    else:
        print("\n'gene_id' column MISSING in var!")

    if 'mapping_source' in adata_v1.var.columns:
        print("\nValue counts for 'mapping_source':")
        print(adata_v1.var['mapping_source'].value_counts())
    else:
        print("\n'mapping_source' column MISSING in var!")
        
    if 'original_ids' in adata_v1.var.columns:
         print("\nSample 'original_ids' column values:")
         print(adata_v1.var['original_ids'].head())
         # Count how many are potentially numeric (Entrez-like) vs Ensembl in the original_ids column
         numeric_orig = sum(1 for x in adata_v1.var['original_ids'] if isinstance(x, str) and all(s.isdigit() for s in x.split(';')))
         ensembl_orig = sum(1 for x in adata_v1.var['original_ids'] if isinstance(x, str) and 'ENSG' in x)
         print(f"Count of numeric-like original_ids: {numeric_orig}")
         print(f"Count of Ensembl-like original_ids: {ensembl_orig}")
    else:
        print("\n'original_ids' column MISSING in var!")

except Exception as e:
    print(f"Error loading or inspecting {adata_v1_path}: {e}")

# --- Load the V2 file ---
adata_v2_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250502_172003/encode_standardized_v2.h5ad'
try:
    adata_v2 = sc.read_h5ad(adata_v2_path)
    print("\n\n--- Stage 2 Output (encode_standardized_v2.h5ad) ---")
    print(f"Shape: {adata_v2.shape}")

    # Check var_names (index)
    print("\nSample var_names (index):")
    print(adata_v2.var_names[:10].tolist())
    try:
        _ = pd.to_numeric(adata_v2.var_names)
        is_numeric_v2 = True
    except ValueError:
        is_numeric_v2 = False
    print(f"Is index potentially numeric? {is_numeric_v2}")

    # Check var dataframe columns
    print("\nVar columns:")
    print(adata_v2.var.columns)

    # --- Check specific columns ---
    if 'gene_id' in adata_v2.var.columns:
        print("\nSample 'gene_id' column values:")
        print(adata_v2.var['gene_id'].head())
    else:
        print("\n'gene_id' column MISSING in var!")

    if 'mapping_source' in adata_v2.var.columns:
        print("\nValue counts for 'mapping_source':")
        print(adata_v2.var['mapping_source'].value_counts())
    else:
        print("\n'mapping_source' column MISSING in var!")
        
    if 'original_ids' in adata_v2.var.columns:
         print("\nSample 'original_ids' column values:")
         print(adata_v2.var['original_ids'].head())
    else:
        print("\n'original_ids' column MISSING in var!")

except Exception as e:
    print(f"Error loading or inspecting {adata_v2_path}: {e}")

