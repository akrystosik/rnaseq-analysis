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

# Run ENCODE processing with just cell_lines directory
log_message "Running: python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py --encode-dir \"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines\" --metadata-dir \"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json\" --output-dir \"$OUTPUT_DIR\""

python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py \
    --encode-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines" \
    --metadata-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json" \
    --output-dir "$OUTPUT_DIR" 2>&1 | tee -a "$LOG_FILE"

# Validate results
log_message "=== Validating Results ==="
log_message "Running: python3 -c \"
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
\""

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
        print('\nTissue distribution:')
        for tissue, count in tissue_counts.items():
            print(f'  {tissue}: {count}')
    
    # Check cell lines
    if 'cell_line' in adata.obs.columns:
        missing_cell_line = pd.isna(adata.obs['cell_line']).sum()
        print(f'\nMissing cell lines (NaN): {missing_cell_line}/{adata.n_obs} ({missing_cell_line/adata.n_obs*100:.2f}%)')
        
        # Show cell line distribution
        cell_line_counts = adata.obs['cell_line'].value_counts(dropna=False).head(10)
        print('\nCell line distribution:')
        for cell_line, count in cell_line_counts.items():
            print(f'  {cell_line}: {count}')
        
    # Check gene statistics
    print('\nGene metadata:')
    if 'gene_type' in adata.var.columns:
        gene_type_counts = adata.var['gene_type'].value_counts(dropna=False).head(10)
        print('\nGene type distribution:')
        for gene_type, count in gene_type_counts.items():
            print(f'  {gene_type}: {count}')
    
    if 'mapping_source' in adata.var.columns:
        mapping_counts = adata.var['mapping_source'].value_counts(dropna=False)
        print('\nMapping source distribution:')
        for source, count in mapping_counts.items():
            print(f'  {source}: {count} ({count/adata.n_vars*100:.2f}%)')
except Exception as e:
    print(f'Error examining ENCODE file: {e}')
" 2>&1 | tee -a "$LOG_FILE"

log_message "=== Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"