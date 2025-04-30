#!/bin/bash
# Direct fix for AnnData string conversion issue

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_ENCODE_with_save_fix_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/encode_save_fix_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname $LOG_FILE)"

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# Start testing
log_message "=== Testing ENCODE Processing with Save Fix ==="
log_message "Output directory: $OUTPUT_DIR"

# Create a simple wrapper script that will fix the string conversion issue
log_message "Creating AnnData save wrapper script"

cat > /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/debug/anndata_save_wrapper.py << 'EOL'
#!/usr/bin/env python3
"""
Wrapper script to save AnnData objects with proper string conversion.
Usage: python anndata_save_wrapper.py input_file output_file
"""

import sys
import os
import scanpy as sc
import pandas as pd
import numpy as np
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('anndata_save_wrapper')

def convert_to_serializable(obj):
    """Convert dict values to serializable types."""
    if isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif pd.isna(obj):
        return ""
    elif isinstance(obj, (pd.Series, pd.DataFrame, np.ndarray)):
        if hasattr(obj, "tolist"):
            return obj.tolist()
        else:
            return str(obj)
    else:
        return str(obj)

def save_adata_safely(input_file, output_file):
    """Load AnnData and save it with proper string conversion."""
    try:
        logger.info(f"Loading AnnData from {input_file}")
        adata = sc.read_h5ad(input_file)
        
        logger.info(f"Loaded AnnData with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Make a copy of the original uns
        original_uns = adata.uns.copy()
        
        # Convert all uns values to serializable types
        logger.info("Converting uns dictionary to serializable types")
        adata.uns = convert_to_serializable(original_uns)
        
        # Save the AnnData object
        logger.info(f"Saving AnnData to {output_file}")
        adata.write(output_file)
        
        # Verify the save worked
        logger.info("Verifying saved file")
        test_adata = sc.read_h5ad(output_file)
        logger.info(f"Verification successful: {test_adata.n_obs} samples, {test_adata.n_vars} genes")
        
        return True
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python anndata_save_wrapper.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    success = save_adata_safely(input_file, output_file)
    
    if success:
        print(f"Successfully saved AnnData to {output_file}")
        sys.exit(0)
    else:
        print("Failed to save AnnData")
        sys.exit(1)
EOL

chmod +x /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/debug/anndata_save_wrapper.py

# Run ENCODE processing with just cell_lines directory
log_message "Running: python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py --encode-dir \"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines\" --metadata-dir \"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json\" --output-dir \"$OUTPUT_DIR/temp\""

# Create temp directory
mkdir -p "$OUTPUT_DIR/temp"

python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py \
    --encode-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines" \
    --metadata-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json" \
    --output-dir "$OUTPUT_DIR/temp" 2>&1 | tee -a "$LOG_FILE"

# Now use our wrapper to save it properly
ENCODE_TEMP="$OUTPUT_DIR/temp/encode_standardized_v1.h5ad"
ENCODE_FIXED="$OUTPUT_DIR/encode_standardized_v1.h5ad"

if [ -f "$ENCODE_TEMP" ]; then
    log_message "Using save wrapper to properly save the AnnData object"
    python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/debug/anndata_save_wrapper.py "$ENCODE_TEMP" "$ENCODE_FIXED" 2>&1 | tee -a "$LOG_FILE"
else
    log_message "ERROR: Temporary file $ENCODE_TEMP not found"
    exit 1
fi

# Validate results
log_message "=== Validating Results ==="
log_message "Running validation script"

python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np

encode_file = '$ENCODE_FIXED'
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

log_message "=== Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"