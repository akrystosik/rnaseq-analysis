#!/bin/bash
# Fix for placeholder gene IDs in preprocessing step

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1"
}

log_message "Implementing fix for placeholder gene IDs"

# Create backup of original file
PREPROCESS_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/preprocess_dataset_gene_ids.py"
cp "$PREPROCESS_FILE" "${PREPROCESS_FILE}.bak"

# Create Python script to modify the file
cat > /tmp/fix_placeholder_ids.py << 'EOL'
import os
import re

# Path to the file
file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/preprocess_dataset_gene_ids.py'

# Read the file
with open(file_path, 'r') as f:
    content = f.read()

# Find the problematic section in preprocess_encode_dataset function
pattern = r"""        # Try reference mapping as last resort
        elif str\(original_id\)\.isdigit\(\) and str\(original_id\) in numeric_to_ensembl:
            ensembl_id = numeric_to_ensembl\[str\(original_id\)\]
            if not ensembl_id\.startswith\('PLACEHOLDER_'\):
                mapping_source = 'reference_mapping'
                mapped_from_reference \+= 1
            else:
                ensembl_id = ''  # Don't use placeholders
                mapping_source = 'unmapped'
                unmapped_count \+= 1"""

# Replacement that simply skips placeholders
replacement = """        # Try reference mapping as last resort
        elif str(original_id).isdigit() and str(original_id) in numeric_to_ensembl:
            ensembl_id = numeric_to_ensembl[str(original_id)]
            # Skip placeholder IDs and keep original numeric ID
            if ensembl_id.startswith('PLACEHOLDER_'):
                # Use original ID as Ensembl ID to avoid empty values
                ensembl_id = f"ENTREZ:{str(original_id)}"
                mapping_source = 'entrez_id'
                mapped_from_reference += 1
            else:
                mapping_source = 'reference_mapping'
                mapped_from_reference += 1"""

# Apply the replacement
updated_content = re.sub(pattern, replacement, content)

# Write the updated content back to the file
with open(file_path, 'w') as f:
    f.write(updated_content)

print("Placeholder ID fix applied to preprocess_dataset_gene_ids.py")
EOL

python /tmp/fix_placeholder_ids.py

log_message "Testing fix with the preprocessing script"

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_preprocess_fix_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/preprocess_fix_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR/input"
mkdir -p "$OUTPUT_DIR/output"
mkdir -p "$(dirname $LOG_FILE)"

# Copy the fixed ENCODE dataset to input
cp "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_ENCODE_with_save_fix_20250430_070735/encode_standardized_v1.h5ad" "$OUTPUT_DIR/input/"

# Run the preprocessing script with our fix
log_message "Running: python $PREPROCESS_FILE --data-dir \"$OUTPUT_DIR/input\" --reference-mapping \"/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_id_reference_mapping.csv\" --output-dir \"$OUTPUT_DIR/output\" --datasets encode"

python $PREPROCESS_FILE \
    --data-dir "$OUTPUT_DIR/input" \
    --reference-mapping "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_id_reference_mapping.csv" \
    --output-dir "$OUTPUT_DIR/output" \
    --datasets encode 2>&1 | tee -a "$LOG_FILE"

# Validate the results
log_message "Validating results"

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
except Exception as e:
    print(f'Error examining preprocessed ENCODE file: {e}')
" 2>&1 | tee -a "$LOG_FILE"

log_message "Fix testing complete. Check the log file for results: $LOG_FILE"
