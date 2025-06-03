#!/bin/bash
# Direct fix for placeholder gene IDs

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1"
}

log_message "Creating direct fix for placeholder gene IDs"

# Create a post-processing script
cat > /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/fix_placeholder_ids.py << 'EOL'
#!/usr/bin/env python3
"""
Fix Placeholder Gene IDs

This script fixes placeholder gene IDs in preprocessed datasets by converting
them to proper Entrez identifiers.

Usage:
    python fix_placeholder_ids.py input_file output_file
"""

import sys
import scanpy as sc
import pandas as pd
import numpy as np
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('placeholder_id_fixer')

def fix_placeholder_ids(input_file, output_file):
    """Fix placeholder gene IDs in a preprocessed dataset."""
    try:
        logger.info(f"Loading dataset from {input_file}")
        adata = sc.read_h5ad(input_file)
        
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Count placeholder IDs before fix
        placeholder_count_before = sum(1 for id in adata.var['ensembl_id'] 
                                      if str(id).startswith('PLACEHOLDER_'))
        
        if placeholder_count_before > 0:
            logger.info(f"Found {placeholder_count_before} placeholder IDs to fix")
            
            # Create a mapping from placeholder IDs to Entrez IDs
            placeholder_to_entrez = {}
            for idx, row in adata.var.iterrows():
                ensembl_id = row['ensembl_id']
                if str(ensembl_id).startswith('PLACEHOLDER_'):
                    # Extract the numeric part from the placeholder ID
                    placeholder_num = str(ensembl_id).replace('PLACEHOLDER_', '')
                    # Use this as the Entrez ID with ENTREZ: prefix
                    placeholder_to_entrez[ensembl_id] = f"ENTREZ:{placeholder_num}"
            
            # Fix the placeholder IDs
            for idx, row in adata.var.iterrows():
                ensembl_id = row['ensembl_id']
                if ensembl_id in placeholder_to_entrez:
                    adata.var.at[idx, 'ensembl_id'] = placeholder_to_entrez[ensembl_id]
                    adata.var.at[idx, 'mapping_source'] = 'entrez_id'
            
            # Count placeholder IDs after fix
            placeholder_count_after = sum(1 for id in adata.var['ensembl_id'] 
                                         if str(id).startswith('PLACEHOLDER_'))
            
            entrez_count = sum(1 for id in adata.var['ensembl_id'] 
                              if str(id).startswith('ENTREZ:'))
            
            logger.info(f"Fix results: {placeholder_count_before} placeholder IDs converted to Entrez IDs")
            logger.info(f"Remaining placeholder IDs: {placeholder_count_after}")
            logger.info(f"Total Entrez IDs: {entrez_count}")
            
            # Save the fixed dataset
            logger.info(f"Saving fixed dataset to {output_file}")
            adata.write(output_file)
            
            logger.info("Fix completed successfully")
            return True
        else:
            logger.info("No placeholder IDs found, no fix needed")
            # Save the dataset as is
            adata.write(output_file)
            return True
    
    except Exception as e:
        logger.error(f"Error fixing placeholder IDs: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_placeholder_ids.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    success = fix_placeholder_ids(input_file, output_file)
    
    if success:
        print(f"Successfully fixed placeholder IDs and saved to {output_file}")
        sys.exit(0)
    else:
        print("Failed to fix placeholder IDs")
        sys.exit(1)
EOL

chmod +x /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/fix_placeholder_ids.py

log_message "Testing placeholder ID fix"

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_placeholder_fix_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/placeholder_fix_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname $LOG_FILE)"

# Copy the preprocessed ENCODE dataset
PREPROCESSED_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_preprocess_20250430_071147/output/encode_standardized_preprocessed.h5ad"
FIXED_FILE="$OUTPUT_DIR/encode_standardized_preprocessed_fixed.h5ad"

log_message "Running fix on preprocessed ENCODE dataset"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/fix_placeholder_ids.py "$PREPROCESSED_FILE" "$FIXED_FILE" 2>&1 | tee -a "$LOG_FILE"

# Validate the fixed dataset
log_message "Validating fixed dataset"

python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np

encode_file = '$FIXED_FILE'
try:
    adata = sc.read_h5ad(encode_file)
    print(f'Fixed ENCODE dataset loaded successfully: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
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
    print(f'Error examining fixed ENCODE file: {e}')
" 2>&1 | tee -a "$LOG_FILE"

# Integrate into pipeline
log_message "Integrating placeholder ID fix into pipeline"

cat > /tmp/update_pipeline_placeholder_fix.py << 'EOL'
import os
import sys
import re

# Path to the run_rnaseq_pipeline.sh file
file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/run_rnaseq_pipeline.sh'

# Read the file
with open(file_path, 'r') as f:
    content = f.read()

# Add placeholder fix commands after preprocessing step
fix_commands = """
# Fix placeholder IDs in preprocessed datasets
log_message "=== Fixing placeholder IDs in preprocessed datasets ==="
for dataset in encode gtex mage adni; do
    preprocessed_file="${PREPROCESSED_DIR}/${dataset}_standardized_preprocessed.h5ad"
    if [ -f "$preprocessed_file" ]; then
        fixed_file="${PREPROCESSED_DIR}/${dataset}_standardized_preprocessed_fixed.h5ad"
        log_message "Fixing placeholder IDs in ${dataset} dataset"
        run_command "python ${SCRIPTS_DIR}/fix_placeholder_ids.py $preprocessed_file $fixed_file"
        if [ $? -eq 0 ]; then
            # Replace the original file with the fixed file
            mv "$fixed_file" "$preprocessed_file"
            log_message "Placeholder IDs fixed in ${dataset} dataset"
        else
            log_message "Warning: Failed to fix placeholder IDs in ${dataset} dataset"
        fi
    fi
done
"""

# Find the appropriate place to add the fix commands (after Step 2.5)
preprocessed_pattern = r'(# Step 2\.5: Preprocess Datasets for Consistent Gene IDs.*?exit 1\n})'
replacement = r'\1\n\n' + fix_commands

# Apply the replacement
updated_content = re.sub(preprocessed_pattern, replacement, content, flags=re.DOTALL)

# Write the updated file
with open(file_path, 'w') as f:
    f.write(updated_content)

print("Pipeline updated with placeholder ID fix")
EOL

python /tmp/update_pipeline_placeholder_fix.py

log_message "Placeholder ID fix successfully integrated into pipeline"
log_message "Fix testing complete. Check the log file for results: $LOG_FILE"