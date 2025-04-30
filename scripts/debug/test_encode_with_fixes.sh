#!/bin/bash
# Test script for ENCODE processing with cell lines only

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_ENCODE_cell_lines_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/encode_test_cell_lines_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname $LOG_FILE)"

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# Start testing
log_message "=== Testing ENCODE Cell Lines Processing Stage ==="
log_message "Output directory: $OUTPUT_DIR"

# Fix string conversion in save_anndata function
log_message "Fixing string conversion in save_anndata function"

cat > /tmp/fix_save_anndata.py << 'EOL'
import os
import sys
import re

# Path to the standardize_datasets.py file
file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py'

# Create a backup if it doesn't exist
backup_path = f"{file_path}.bak"
if not os.path.exists(backup_path):
    print(f"Creating backup at {backup_path}")
    with open(file_path, 'r') as f_in:
        with open(backup_path, 'w') as f_out:
            f_out.write(f_in.read())

# Read the file
with open(file_path, 'r') as f:
    content = f.read()

# Add string conversion to create_standard_anndata function
create_anndata_pattern = r'def create_standard_anndata\(data_df, obs_df, var_df, dataset_info\):'
create_anndata_fix = """def create_standard_anndata(data_df, obs_df, var_df, dataset_info):
    # Ensure all obs columns can be properly saved
    for col in obs_df.columns:
        if obs_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting obs column {col} to string")
            obs_df[col] = obs_df[col].astype(str)

    # Ensure all var columns can be properly saved
    for col in var_df.columns:
        if var_df[col].dtype.name not in ["category", "string", "object"]:
            logger.debug(f"Converting var column {col} to string")
            var_df[col] = var_df[col].astype(str)

    # Convert all dictionary values to strings
    safe_dataset_info = {}
    for key, value in dataset_info.items():
        if isinstance(value, dict):
            safe_dataset_info[key] = {k: str(v) if v is not None else "" for k, v in value.items()}
        elif value is None:
            safe_dataset_info[key] = ""
        else:
            safe_dataset_info[key] = str(value)"""

# Add string conversion to save_anndata function
save_anndata_pattern = r'def save_anndata\(adata, file_path\):'
save_anndata_fix = """def save_anndata(adata, file_path):
    \"\"\"Save AnnData object to file with validation.\"\"\"
    try:
        # Validate before saving
        if adata.var.shape[1] == 0:
            logger.error("Cannot save AnnData: var DataFrame is empty!")
            return False

        # Check for index/column name conflicts and fix them
        if adata.obs.index.name is not None and adata.obs.index.name in adata.obs.columns:
            logger.warning(
                f"Fixing index/column name conflict: Renaming column '{adata.obs.index.name}' to 'original_{adata.obs.index.name}'"
            )
            adata.obs = adata.obs.rename(
                columns={adata.obs.index.name: f"original_{adata.obs.index.name}"}
            )

        if adata.var.index.name is not None and adata.var.index.name in adata.var.columns:
            logger.warning(
                f"Fixing index/column name conflict: Renaming column '{adata.var.index.name}' to 'original_{adata.var.index.name}'"
            )
            adata.var = adata.var.rename(
                columns={adata.var.index.name: f"original_{adata.var.index.name}"}
            )

        # Convert all values in uns to serializable types
        safe_uns = {}
        for key, value in adata.uns.items():
            if isinstance(value, dict):
                safe_uns[key] = {k: str(v) if v is not None else "" for k, v in value.items()}
            elif value is None:
                safe_uns[key] = ""
            else:
                safe_uns[key] = str(value)

        # Store original uns
        original_uns = adata.uns.copy()
        
        # Replace with safe uns for saving
        adata.uns = safe_uns

        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        # Save the AnnData object
        adata.write_h5ad(file_path)
        
        # Restore original uns
        adata.uns = original_uns

        # Verify successful save by loading
        test_load = ad.read_h5ad(file_path)
        logger.info(
            f"Verification: Saved AnnData loaded with var shape={test_load.var.shape}, obs shape={test_load.obs.shape}"
        )

        logger.info(f"Saved AnnData to {file_path}")
        return True

    except Exception as e:
        logger.error(f"Error saving AnnData: {e}")
        return False"""

# Apply the fixes
content = re.sub(create_anndata_pattern, create_anndata_fix, content)
content = re.sub(save_anndata_pattern, save_anndata_fix, content)

# Write the fixed file
with open(file_path, 'w') as f:
    f.write(content)

print("String conversion fixes applied successfully.")
EOL

python /tmp/fix_save_anndata.py

# Run ENCODE processing with just cell_lines directory
log_message "Running: python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py --encode-dir \"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines\" --metadata-dir \"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json\" --output-dir \"$OUTPUT_DIR\""

python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py \
    --encode-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines" \
    --metadata-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json" \
    --output-dir "$OUTPUT_DIR" 2>&1 | tee -a "$LOG_FILE"

# Validate results
log_message "=== Validating Results ==="
log_message "Running validation script"

python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np

encode_file = '$OUTPUT_DIR/encode_standardized_v1.h5ad'
try:
    adata = sc.read_h5ad(encode_file)
    print(f'ENCODE dataset loaded successfully: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check for tissue information
    if 'tissue' in adata.obs.columns:
        # Check for true NaN values
        missing_tissue = pd.isna(adata.obs['tissue']).sum()
        print(f'Missing tissues (NaN): {missing_tissue}/{adata.n_obs} ({missing_tissue/adata.n_obs*100:.2f}%)')
        
        # Check for 'nan' string values
        nan_string = sum(adata.obs['tissue'].astype(str) == 'nan')
        print(f'String \'nan\' tissues: {nan_string}/{adata.n_obs} ({nan_string/adata.n_obs*100:.2f}%)')
        
        # Check for empty strings
        empty_string = sum(adata.obs['tissue'].astype(str) == '')
        print(f'Empty string tissues: {empty_string}/{adata.n_obs} ({empty_string/adata.n_obs*100:.2f}%)')
        
        # Show tissue distribution
        tissue_counts = adata.obs['tissue'].value_counts(dropna=False).head(10)
        print('\\nTissue distribution:')
        for tissue, count in tissue_counts.items():
            print(f'  {tissue}: {count}')
    
    # Check cell lines
    if 'cell_line' in adata.obs.columns:
        missing_cell_line = pd.isna(adata.obs['cell_line']).sum()
        print(f'\\nMissing cell lines (NaN): {missing_cell_line}/{adata.n_obs} ({missing_cell_line/adata.n_obs*100:.2f}%)')
        
        # Show cell line distribution
        cell_line_counts = adata.obs['cell_line'].value_counts(dropna=False).head(10)
        print('\\nCell line distribution:')
        for cell_line, count in cell_line_counts.items():
            print(f'  {cell_line}: {count}')
        
    # Check gene statistics
    print('\\nGene metadata:')
    if 'gene_type' in adata.var.columns:
        gene_type_counts = adata.var['gene_type'].value_counts(dropna=False).head(10)
        print('\\nGene type distribution:')
        for gene_type, count in gene_type_counts.items():
            print(f'  {gene_type}: {count}')
    
    if 'mapping_source' in adata.var.columns:
        mapping_counts = adata.var['mapping_source'].value_counts(dropna=False)
        print('\\nMapping source distribution:')
        for source, count in mapping_counts.items():
            print(f'  {source}: {count} ({count/adata.n_vars*100:.2f}%)')
except Exception as e:
    print(f'Error examining ENCODE file: {e}')
" 2>&1 | tee -a "$LOG_FILE"

# Update pipeline script to use cell_lines directory
log_message "Updating pipeline script to use cell_lines directory"

cat > /tmp/update_pipeline.py << 'EOL'
import os
import sys
import re

# Path to the run_rnaseq_pipeline.sh file
file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/run_rnaseq_pipeline.sh'

# Create a backup if it doesn't exist
backup_path = f"{file_path}.bak"
if not os.path.exists(backup_path):
    print(f"Creating backup at {backup_path}")
    with open(file_path, 'r') as f_in:
        with open(backup_path, 'w') as f_out:
            f_out.write(f_in.read())

# Read the file
with open(file_path, 'r') as f:
    content = f.read()

# Find and replace the ENCODE directory path
content = content.replace(
    '--encode-dir "${BASE_DIR}/encode/raw_data"',
    '--encode-dir "${BASE_DIR}/encode/raw_data/cell_lines"'
)

# Write the fixed file
with open(file_path, 'w') as f:
    f.write(content)

print("Pipeline ENCODE path updated successfully.")
EOL

python /tmp/update_pipeline.py

log_message "=== Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"