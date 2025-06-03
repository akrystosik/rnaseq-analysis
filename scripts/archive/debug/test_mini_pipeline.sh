#!/bin/bash
# Test script for running a small subset of the pipeline

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_mini_pipeline_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/test_mini_pipeline_${TIMESTAMP}.log"

# Create log directory
mkdir -p "$(dirname $LOG_FILE)"

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

log_message "=== Running Mini Pipeline Test ==="
log_message "Output directory: $OUTPUT_DIR"

# Create subdirectories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/preprocessed_data"

# Step 1: Run standardize_datasets.py for ENCODE only
log_message "Step 1: Running standardize_datasets.py for ENCODE only"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py \
    --encode-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines" \
    --metadata-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json" \
    --output-dir "$OUTPUT_DIR" 2>&1 | tee -a "$LOG_FILE"

# Use save wrapper to properly save AnnData object
log_message "Using save wrapper for ENCODE dataset"
ENCODE_TEMP="$OUTPUT_DIR/encode_standardized_v1.h5ad"
if [ -f "$ENCODE_TEMP" ]; then
    python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/anndata_save_wrapper.py \
        "$ENCODE_TEMP" "$ENCODE_TEMP.fixed" 2>&1 | tee -a "$LOG_FILE"
    mv "$ENCODE_TEMP.fixed" "$ENCODE_TEMP"
else
    log_message "ERROR: ENCODE file not found"
    exit 1
fi

# Step 2: Run preprocess_dataset_gene_ids.py
log_message "Step 2: Running preprocess_dataset_gene_ids.py"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/preprocess_dataset_gene_ids.py \
    --data-dir "$OUTPUT_DIR" \
    --reference-mapping "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_id_reference_mapping.csv" \
    --output-dir "$OUTPUT_DIR/preprocessed_data" \
    --datasets encode 2>&1 | tee -a "$LOG_FILE"

# Step 3: Fix placeholder IDs
log_message "Step 3: Fixing placeholder IDs"
PREPROCESSED_FILE="$OUTPUT_DIR/preprocessed_data/encode_standardized_preprocessed.h5ad"
if [ -f "$PREPROCESSED_FILE" ]; then
    python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/fix_placeholder_ids.py \
        "$PREPROCESSED_FILE" "$PREPROCESSED_FILE.fixed" 2>&1 | tee -a "$LOG_FILE"
    mv "$PREPROCESSED_FILE.fixed" "$PREPROCESSED_FILE"
else
    log_message "ERROR: Preprocessed ENCODE file not found"
    exit 1
fi

# Step 4: Analyze results
log_message "Step 4: Analyzing results"

# Check original dataset
log_message "Analyzing original ENCODE dataset"
python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

original_file = '$OUTPUT_DIR/encode_standardized_v1.h5ad'
try:
    adata = sc.read_h5ad(original_file)
    print(f'Original ENCODE dataset loaded: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check tissue and cell line distribution
    print('\\nTissue distribution:')
    print(adata.obs['tissue'].value_counts())
    
    print('\\nCell line distribution:')
    print(adata.obs['cell_line'].value_counts())
    
    # Check gene mapping
    if 'mapping_source' in adata.var.columns:
        mapping_counts = adata.var['mapping_source'].value_counts(dropna=False)
        print('\\nMapping source distribution:')
        for source, count in mapping_counts.items():
            print(f'  {source}: {count} ({count/adata.n_vars*100:.2f}%)')
            
    # Check columns in obs and var
    print('\\nSample metadata columns:')
    for col in adata.obs.columns:
        print(f'  {col}')
        
    print('\\nGene metadata columns:')
    for col in adata.var.columns:
        print(f'  {col}')
    
except Exception as e:
    print(f'Error examining original file: {e}')
" 2>&1 | tee -a "$LOG_FILE"

# Check preprocessed dataset
log_message "Analyzing preprocessed ENCODE dataset"
python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

preprocessed_file = '$OUTPUT_DIR/preprocessed_data/encode_standardized_preprocessed.h5ad'
try:
    adata = sc.read_h5ad(preprocessed_file)
    print(f'Preprocessed ENCODE dataset loaded: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check for placeholder IDs
    if 'ensembl_id' in adata.var.columns:
        placeholder_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('PLACEHOLDER_'))
        print(f'Placeholder IDs found: {placeholder_count}/{len(adata.var)} ({placeholder_count/len(adata.var)*100:.2f}%)')
        
        entrez_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('ENTREZ:'))
        print(f'ENTREZ: format IDs found: {entrez_count}/{len(adata.var)} ({entrez_count/len(adata.var)*100:.2f}%)')
    
    # Check mapping source distribution
    if 'mapping_source' in adata.var.columns:
        mapping_counts = adata.var['mapping_source'].value_counts(dropna=False)
        print('\\nMapping source distribution:')
        for source, count in mapping_counts.items():
            print(f'  {source}: {count} ({count/adata.n_vars*100:.2f}%)')
    
    # Check ensembl_id distribution
    ensembl_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('ENSG'))
    empty_count = sum(1 for id in adata.var['ensembl_id'] if pd.isna(id) or str(id) == '')
    other_count = len(adata.var) - ensembl_count - empty_count - entrez_count
    
    print('\\nEnsembl ID distribution:')
    print(f'  Valid Ensembl IDs: {ensembl_count}/{len(adata.var)} ({ensembl_count/len(adata.var)*100:.2f}%)')
    print(f'  ENTREZ: format IDs: {entrez_count}/{len(adata.var)} ({entrez_count/len(adata.var)*100:.2f}%)')
    print(f'  Empty IDs: {empty_count}/{len(adata.var)} ({empty_count/len(adata.var)*100:.2f}%)')
    print(f'  Other IDs: {other_count}/{len(adata.var)} ({other_count/len(adata.var)*100:.2f}%)')
    
    # Sample data from obs and var
    print('\\nSample metadata preview:')
    print(adata.obs.head())
    
    print('\\nGene metadata preview:')
    print(adata.var.head())
except Exception as e:
    print(f'Error examining preprocessed file: {e}')
" 2>&1 | tee -a "$LOG_FILE"

log_message "=== Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"