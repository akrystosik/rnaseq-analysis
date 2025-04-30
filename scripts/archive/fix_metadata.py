# First, let's check the structure of the GTEx dataset
def inspect_gtex_metadata():
    """Inspect the structure of the GTEx AnnData object to identify missing metadata."""
    import anndata as ad
    import pandas as pd
    
    gtex_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/gtex_standardized.h5ad"
    adata = ad.read_h5ad(gtex_path)
    
    # Check variable (gene) metadata
    print(f"GTEx variable (gene) metadata columns: {adata.var.columns.tolist()}")
    
    # Check observation (sample) metadata 
    print(f"GTEx observation (sample) metadata columns: {adata.obs.columns.tolist()}")
    
    # Check if any tissue-related columns exist but have a different name
    potential_tissue_cols = [col for col in adata.obs.columns if 'tissue' in col.lower() 
                           or 'site' in col.lower() or 'organ' in col.lower()]
    print(f"Potential tissue columns: {potential_tissue_cols}")
    
    # Look at uns (unstructured metadata)
    print(f"GTEx unstructured metadata keys: {adata.uns.keys()}")
    
    # If there's a potential tissue column, look at the first few entries
    if potential_tissue_cols:
        for col in potential_tissue_cols:
            print(f"\nSample values from {col}:")
            print(adata.obs[col].value_counts().head())
    
    return adata

def fix_empty_categories():
    """Fix empty category values in the datasets."""
    import anndata as ad
    import pandas as pd
    import numpy as np
    
    # Process each dataset
    for dataset in ["encode", "adni", "mage"]:
        path = f"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/{dataset}_standardized.h5ad"
        adata = ad.read_h5ad(path)
        
        # Fix tissue column
        tissue_col = None
        for col in adata.obs.columns:
            if "tissue" in col.lower():
                tissue_col = col
                break
        
        if tissue_col:
            # Replace nan and empty values with "unknown"
            mask = adata.obs[tissue_col].isna() | (adata.obs[tissue_col] == "")
            if mask.any():
                print(f"Replacing {mask.sum()} empty tissue values in {dataset}")
                adata.obs[tissue_col] = adata.obs[tissue_col].fillna("unknown")
                adata.obs.loc[adata.obs[tissue_col] == "", tissue_col] = "unknown"
                
                # Write back the fixed data
                adata.write_h5ad(path)
                print(f"Updated {dataset} dataset")
                
                
def fix_encode_subjects():
    """Fix empty subject labels in ENCODE dataset."""
    import anndata as ad
    import pandas as pd
    
    path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/encode_standardized.h5ad"
    adata = ad.read_h5ad(path)
    
    # Find subject column
    subject_col = None
    for col in adata.obs.columns:
        if "subject" in col.lower() or "donor" in col.lower():
            subject_col = col
            break
    
    if subject_col:
        # Check for empty values
        empty_mask = (adata.obs[subject_col] == "") | adata.obs[subject_col].isna()
        
        if empty_mask.any():
            print(f"Found {empty_mask.sum()} empty subject values in ENCODE")
            
            # Try to infer subject from cell line information
            if 'cell_line' in adata.obs.columns:
                for idx in adata.obs[empty_mask].index:
                    cell_line = adata.obs.loc[idx, 'cell_line']
                    if pd.notna(cell_line) and cell_line != "":
                        # Use the cell line info as a fallback
                        adata.obs.loc[idx, subject_col] = f"derived_from_{cell_line}"
            
            # For remaining empty values, use a placeholder
            still_empty = (adata.obs[subject_col] == "") | adata.obs[subject_col].isna()
            if still_empty.any():
                adata.obs.loc[still_empty, subject_col] = "unknown_subject"
            
            # Write back the fixed data
            adata.write_h5ad(path)
            print("Updated ENCODE dataset")                
            
            
if __name__ == '__main__':
    print("Running GTEx metadata inspection...")
    inspect_gtex_metadata()
    
    print("\nFixing empty category values...")
    fix_empty_categories()
    
    print("\nFixing ENCODE subject values...")
    fix_encode_subjects()
    
    print("\nMetadata fixes completed!")