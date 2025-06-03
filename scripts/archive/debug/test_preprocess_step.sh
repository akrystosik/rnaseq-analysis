#!/bin/bash
# Test script for preprocessing step to fix placeholder gene IDs

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_preprocess_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/preprocess_test_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname $LOG_FILE)"

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# Start testing
log_message "=== Testing Gene ID Preprocessing ==="
log_message "Output directory: $OUTPUT_DIR"

# Copy our fixed ENCODE dataset to the input directory
mkdir -p "$OUTPUT_DIR/input"
log_message "Copying fixed ENCODE dataset as input"
cp "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_ENCODE_with_save_fix_20250430_070735/encode_standardized_v1.h5ad" "$OUTPUT_DIR/input/"

# Run the preprocessing script on the ENCODE dataset
log_message "Running: python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/preprocess_dataset_gene_ids.py --data-dir \"$OUTPUT_DIR/input\" --reference-mapping \"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_id_reference_mapping.csv\" --output-dir \"$OUTPUT_DIR/output\" --datasets encode"

python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/preprocess_dataset_gene_ids.py \
    --data-dir "$OUTPUT_DIR/input" \
    --reference-mapping "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_id_reference_mapping.csv" \
    --output-dir "$OUTPUT_DIR/output" \
    --datasets encode 2>&1 | tee -a "$LOG_FILE"

# Validate results
log_message "=== Validating Results ==="
log_message "Running validation script"

python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np

encode_file = '$OUTPUT_DIR/output/encode_standardized_preprocessed.h5ad'
try:
    adata = sc.read_h5ad(encode_file)
    print(f'Preprocessed ENCODE dataset loaded successfully: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check for placeholder IDs
    if 'ensembl_id' in adata.var.columns:
        placeholder_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('PLACEHOLDER_'))
        placeholder_percent = placeholder_count / len(adata.var) * 100
        print(f'Placeholder IDs found: {placeholder_count}/{len(adata.var)} ({placeholder_percent:.2f}%)')
        
        # Show examples of placeholder IDs
        if placeholder_count > 0:
            placeholder_examples = [id for id in adata.var['ensembl_id'] if str(id).startswith('PLACEHOLDER_')][:5]
            print('\\nExamples of placeholder IDs:')
            for id in placeholder_examples:
                print(f'  {id}')
    
    # Check mapping source distribution
    if 'mapping_source' in adata.var.columns:
        mapping_counts = adata.var['mapping_source'].value_counts(dropna=False)
        print('\\nMapping source distribution:')
        for source, count in mapping_counts.items():
            print(f'  {source}: {count} ({count/adata.n_vars*100:.2f}%)')
    
    # Check ensembl_id distribution
    ensembl_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('ENSG'))
    empty_count = sum(1 for id in adata.var['ensembl_id'] if pd.isna(id) or str(id) == '')
    other_count = len(adata.var) - ensembl_count - empty_count - placeholder_count
    
    print('\\nEnsembl ID distribution:')
    print(f'  Valid Ensembl IDs: {ensembl_count}/{len(adata.var)} ({ensembl_count/len(adata.var)*100:.2f}%)')
    print(f'  Placeholder IDs: {placeholder_count}/{len(adata.var)} ({placeholder_count/len(adata.var)*100:.2f}%)')
    print(f'  Empty IDs: {empty_count}/{len(adata.var)} ({empty_count/len(adata.var)*100:.2f}%)')
    print(f'  Other IDs: {other_count}/{len(adata.var)} ({other_count/len(adata.var)*100:.2f}%)')
except Exception as e:
    print(f'Error examining preprocessed ENCODE file: {e}')
" 2>&1 | tee -a "$LOG_FILE"

log_message "=== Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"