#!/bin/bash
# Script to apply patches to fix AnnData string conversion issues

# Create a backup of the original file
cp /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py.bak

# Apply the patch to fix string conversion in create_standard_anndata
cat > /tmp/anndata_string_fix.patch << 'EOL'
--- standardize_datasets.py	2025-04-30 07:10:00.000000000 +0000
+++ standardize_datasets.py.fixed	2025-04-30 07:10:00.000000000 +0000
@@ -113,7 +113,22 @@
     # Ensure all uns values are serializable
     for key in dataset_info.keys():
-        if isinstance(dataset_info[key], (pd.Series, pd.DataFrame, np.ndarray)):
+        # Handle various types that need conversion
+        if isinstance(dataset_info[key], dict):
+            # Convert all dict values to strings
+            for subkey in dataset_info[key]:
+                if isinstance(dataset_info[key][subkey], (pd.Series, pd.DataFrame, np.ndarray, np.integer, np.floating)):
+                    dataset_info[key][subkey] = str(dataset_info[key][subkey])
+                elif pd.isna(dataset_info[key][subkey]):
+                    dataset_info[key][subkey] = ''
+        elif isinstance(dataset_info[key], (pd.Series, pd.DataFrame, np.ndarray)):
             logger.debug(
                 f"Converting uns[{key}] from {type(dataset_info[key])} to standard Python type"
             )
@@ -121,6 +136,9 @@
                 dataset_info[key].tolist()
                 if hasattr(dataset_info[key], "tolist")
                 else str(dataset_info[key])
             )
+        elif isinstance(dataset_info[key], (np.integer, np.floating)):
+            dataset_info[key] = dataset_info[key].item()
+        elif pd.isna(dataset_info[key]):
+            dataset_info[key] = ''
 
     # Validate metadata before preparing for AnnData
EOL

# Apply the patch
patch /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py /tmp/anndata_string_fix.patch

# Also add a fix to the save_anndata function 
cat > /tmp/save_anndata_fix.patch << 'EOL'
--- standardize_datasets.py	2025-04-30 07:10:00.000000000 +0000
+++ standardize_datasets.py.fixed	2025-04-30 07:10:00.000000000 +0000
@@ -208,6 +208,21 @@
             logger.warning(
                 f"Fixing index/column name conflict: Renaming column '{adata.var.index.name}' to 'original_{adata.var.index.name}'"
             )
             adata.var = adata.var.rename(
                 columns={adata.var.index.name: f"original_{adata.var.index.name}"}
             )
+
+        # Ensure all values in uns are serializable
+        for key in list(adata.uns.keys()):
+            if isinstance(adata.uns[key], dict):
+                # Convert all values in the dictionary to serializable types
+                for subkey in list(adata.uns[key].keys()):
+                    if isinstance(adata.uns[key][subkey], (np.integer, np.floating)):
+                        adata.uns[key][subkey] = adata.uns[key][subkey].item()
+                    elif pd.isna(adata.uns[key][subkey]):
+                        adata.uns[key][subkey] = ''
+            elif isinstance(adata.uns[key], (np.integer, np.floating)):
+                adata.uns[key] = adata.uns[key].item()
+            elif pd.isna(adata.uns[key]):
+                adata.uns[key] = ''
 
         # Create directory if it doesn't exist
         os.makedirs(os.path.dirname(file_path), exist_ok=True)
EOL

# Apply the save_anndata fix
patch /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py /tmp/save_anndata_fix.patch

echo "Patches applied successfully."