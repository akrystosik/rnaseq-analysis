#!/bin/bash
# Test just the ENCODE processing stage

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR="${BASE_DIR}/scripts"
METADATA_DIR="${BASE_DIR}/metadata/json"

# Create a timestamp for output directory
TIMESTAMP="ENCODE_$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="${BASE_DIR}/standardized_data/test_${TIMESTAMP}"
LOG_FILE="${BASE_DIR}/logs/encode_test_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname "$LOG_FILE")"

# Log function
log() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# Run function with logging
run() {
    log "Running: $1"
    eval "$1" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

log "=== Testing ENCODE Processing Stage ==="
log "Output directory: $OUTPUT_DIR"

# Process ENCODE dataset
run "python ${SCRIPTS_DIR}/standardize_datasets.py --encode-dir \"${BASE_DIR}/encode/raw_data\" --metadata-dir \"${METADATA_DIR}\" --output-dir \"$OUTPUT_DIR\""

# Validate results
log "=== Validating Results ==="
run "python3 -c \"
import scanpy as sc
import pandas as pd
import numpy as np

encode_file = '${OUTPUT_DIR}/encode_standardized_v1.h5ad'
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
        print(f'String \\'nan\\' tissues: {nan_string}/{adata.n_obs} ({nan_string/adata.n_obs*100:.2f}%)')
        
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
except Exception as e:
    print(f'Error examining ENCODE file: {e}')
\""

log "=== Test Complete ==="
log "Results saved to $OUTPUT_DIR"
log "Log file: $LOG_FILE"
