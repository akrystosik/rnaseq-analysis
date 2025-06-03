#!/bin/bash
# Integration of AnnData save wrapper into the pipeline

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1"
}

# Step 1: Copy save wrapper to scripts directory
log_message "Installing AnnData save wrapper to scripts directory"

WRAPPER_SRC="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/debug/anndata_save_wrapper.py"
WRAPPER_DST="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/anndata_save_wrapper.py"

cp "$WRAPPER_SRC" "$WRAPPER_DST"
chmod +x "$WRAPPER_DST"

log_message "Save wrapper installed at $WRAPPER_DST"

# Step 2: Update run_rnaseq_pipeline.sh to use the wrapper
log_message "Updating pipeline script to use save wrapper"

cat > /tmp/update_pipeline_save.py << 'EOL'
import os
import sys
import re

# Path to the run_rnaseq_pipeline.sh file
file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/run_rnaseq_pipeline.sh'

# Read the file
with open(file_path, 'r') as f:
    content = f.read()

# Define the save wrapper function
save_wrapper_function = """
# Function to safely save AnnData files
save_anndata_safely() {
    local input_file="$1"
    local output_file="$2"
    
    if [ -f "$input_file" ]; then
        log_message "Using save wrapper to safely save AnnData: $output_file"
        python "${SCRIPTS_DIR}/anndata_save_wrapper.py" "$input_file" "$output_file"
        return $?
    else
        log_message "Error: Input file not found: $input_file"
        return 1
    fi
}
"""

# Add the function after the run_command function
content = content.replace(
    "# Function to run a command with logging",
    "# Function to run a command with logging\n" + save_wrapper_function
)

# Find and replace standardize_datasets.py call to use temporary output and save wrapper
new_output_path = "${OUTPUT_DIR}/temp"
temp_output_pattern = r'--output-dir "\$OUTPUT_DIR"'
temp_output_replacement = f'--output-dir "{new_output_path}"'

# Replace the output directory in standardize_datasets.py call
content = content.replace(
    temp_output_pattern,
    temp_output_replacement
)

# Add mkdir commands to create temp directory
mkdir_command = "\n# Create temp directory for intermediate files\nmkdir -p \"${OUTPUT_DIR}/temp\"\n"
content = content.replace(
    "# Step 1: Initial Data Conversion (existing)",
    mkdir_command + "# Step 1: Initial Data Conversion (existing)"
)

# Add save wrapper calls after each standardize_datasets.py call
save_wrapper_calls = """
# Use save wrapper to properly save AnnData files
if [ -f "${OUTPUT_DIR}/temp/encode_standardized_v1.h5ad" ]; then
    save_anndata_safely "${OUTPUT_DIR}/temp/encode_standardized_v1.h5ad" "${OUTPUT_DIR}/encode_standardized_v1.h5ad"
fi
if [ -f "${OUTPUT_DIR}/temp/gtex_standardized_v1.h5ad" ]; then
    save_anndata_safely "${OUTPUT_DIR}/temp/gtex_standardized_v1.h5ad" "${OUTPUT_DIR}/gtex_standardized_v1.h5ad"
fi
if [ -f "${OUTPUT_DIR}/temp/mage_standardized_v1.h5ad" ]; then
    save_anndata_safely "${OUTPUT_DIR}/temp/mage_standardized_v1.h5ad" "${OUTPUT_DIR}/mage_standardized_v1.h5ad"
fi
if [ -f "${OUTPUT_DIR}/temp/adni_standardized_v1.h5ad" ]; then
    save_anndata_safely "${OUTPUT_DIR}/temp/adni_standardized_v1.h5ad" "${OUTPUT_DIR}/adni_standardized_v1.h5ad"
fi
"""

# Add save wrapper calls after standardize_datasets.py
content = re.sub(
    r'(python \${SCRIPTS_DIR}/standardize_datasets\.py .*?fi\n)',
    r'\1\n' + save_wrapper_calls,
    content,
    flags=re.DOTALL
)

# Write the updated file
with open(file_path, 'w') as f:
    f.write(content)

print("Pipeline script updated to use save wrapper")
EOL

python /tmp/update_pipeline_save.py

log_message "Pipeline updated to use save wrapper"

# Step 3: Create test script for full pipeline run
log_message "Creating test script for full pipeline run"

cat > /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/debug/test_full_pipeline.sh << 'EOL'
#!/bin/bash
# Test script for full pipeline run

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_full_pipeline_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/test_full_pipeline_${TIMESTAMP}.log"

# Create log directory
mkdir -p "$(dirname $LOG_FILE)"

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

log_message "=== Testing Full Pipeline with Fixes ==="
log_message "Output directory: $OUTPUT_DIR"

# Run the pipeline with the specific output directory
log_message "Running pipeline script"
/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/run_rnaseq_pipeline.sh "$OUTPUT_DIR" 2>&1 | tee -a "$LOG_FILE"

# Validate the results
log_message "=== Validating Results ==="

# Check ENCODE dataset
log_message "Validating ENCODE dataset"
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
        
        # Show tissue distribution
        tissue_counts = adata.obs['tissue'].value_counts(dropna=False).head(5)
        print('\\nTissue distribution (top 5):')
        for tissue, count in tissue_counts.items():
            print(f'  {tissue}: {count}')
    
    # Check cell lines
    if 'cell_line' in adata.obs.columns:
        missing_cell_line = pd.isna(adata.obs['cell_line']).sum()
        print(f'\\nMissing cell lines (NaN): {missing_cell_line}/{adata.n_obs} ({missing_cell_line/adata.n_obs*100:.2f}%)')
        
        # Show cell line distribution
        cell_line_counts = adata.obs['cell_line'].value_counts(dropna=False).head(5)
        print('\\nCell line distribution (top 5):')
        for cell_line, count in cell_line_counts.items():
            print(f'  {cell_line}: {count}')
    
    # Check gene mapping
    if 'mapping_source' in adata.var.columns:
        mapping_counts = adata.var['mapping_source'].value_counts(dropna=False)
        print('\\nMapping source distribution:')
        for source, count in mapping_counts.items():
            print(f'  {source}: {count} ({count/adata.n_vars*100:.2f}%)')
except Exception as e:
    print(f'Error examining ENCODE file: {e}')
" 2>&1 | tee -a "$LOG_FILE"

# Check combined dataset if it exists
combined_file="$OUTPUT_DIR/combined_all_genes_sparse_standardized.h5ad"
if [ -f "$combined_file" ]; then
    log_message "Validating combined dataset"
    python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np

combined_file = '$combined_file'
try:
    adata = sc.read_h5ad(combined_file)
    print(f'Combined dataset loaded successfully: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check dataset distribution
    if 'dataset' in adata.obs.columns:
        dataset_counts = adata.obs['dataset'].value_counts(dropna=False)
        print('\\nDataset distribution:')
        for dataset, count in dataset_counts.items():
            print(f'  {dataset}: {count} ({count/adata.n_obs*100:.2f}%)')
    
    # Check gene overlap
    if 'present_in_datasets' in adata.var.columns:
        # Count genes present in each dataset
        datasets = ['encode', 'gtex', 'mage', 'adni']
        for dataset in datasets:
            count = sum(1 for x in adata.var['present_in_datasets'] if dataset in str(x))
            percentage = count / adata.n_vars * 100
            print(f'\\nGenes present in {dataset}: {count}/{adata.n_vars} ({percentage:.2f}%)')
        
        # Count genes present in all datasets
        all_datasets = sum(1 for x in adata.var['present_in_datasets'] 
                           if all(dataset in str(x) for dataset in datasets))
        print(f'\\nGenes present in all datasets: {all_datasets}/{adata.n_vars} ({all_datasets/adata.n_vars*100:.2f}%)')
    
    # Check gene ID formats
    ensembl_count = sum(1 for x in adata.var_names if str(x).startswith('ENSG'))
    other_count = adata.n_vars - ensembl_count
    print(f'\\nGene ID format distribution:')
    print(f'  Ensembl IDs: {ensembl_count}/{adata.n_vars} ({ensembl_count/adata.n_vars*100:.2f}%)')
    print(f'  Other IDs: {other_count}/{adata.n_vars} ({other_count/adata.n_vars*100:.2f}%)')
except Exception as e:
    print(f'Error examining combined dataset: {e}')
" 2>&1 | tee -a "$LOG_FILE"
else
    log_message "Combined dataset not found at $combined_file"
fi

log_message "=== Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"
EOL

chmod +x /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/debug/test_full_pipeline.sh

log_message "Test script created at /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/debug/test_full_pipeline.sh"
log_message "Integration complete. You can now run the test script to validate the full pipeline."