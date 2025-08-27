#!/bin/bash
# Fix for categorical columns when handling placeholder IDs

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1"
}

log_message "Creating improved fix for placeholder gene IDs (handling categorical data)"

# Create an improved post-processing script
cat > /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/fix_placeholder_ids.py << 'EOL'
#!/usr/bin/env python3
"""
Fix Placeholder Gene IDs

This script fixes placeholder gene IDs in preprocessed datasets by converting
them to proper Entrez identifiers, properly handling categorical data types.

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
            
            # First handle categorical columns by converting to string
            var_df = adata.var.copy()
            
            # Convert categorical columns to strings
            for col in var_df.columns:
                if pd.api.types.is_categorical_dtype(var_df[col]):
                    logger.info(f"Converting categorical column {col} to string")
                    var_df[col] = var_df[col].astype(str)
            
            # Create a mapping from placeholder IDs to Entrez IDs
            placeholder_to_entrez = {}
            placeholder_indices = []
            
            # Identify placeholder IDs and create mapping
            for idx, row in var_df.iterrows():
                ensembl_id = row['ensembl_id']
                if str(ensembl_id).startswith('PLACEHOLDER_'):
                    # Extract the numeric part from the placeholder ID
                    placeholder_num = str(ensembl_id).replace('PLACEHOLDER_', '')
                    # Use this as the Entrez ID with ENTREZ: prefix
                    entrez_id = f"ENTREZ:{placeholder_num}"
                    placeholder_to_entrez[ensembl_id] = entrez_id
                    placeholder_indices.append(idx)
            
            # Fix the placeholder IDs
            logger.info(f"Replacing {len(placeholder_indices)} placeholder IDs with Entrez IDs")
            for idx in placeholder_indices:
                old_id = var_df.at[idx, 'ensembl_id']
                var_df.at[idx, 'ensembl_id'] = placeholder_to_entrez[old_id]
                var_df.at[idx, 'mapping_source'] = 'entrez_id'
            
            # Replace adata.var with the fixed DataFrame
            adata.var = var_df
            
            # Count IDs after fix
            entrez_count = sum(1 for id in adata.var['ensembl_id'] 
                              if str(id).startswith('ENTREZ:'))
            
            logger.info(f"Fix results: {placeholder_count_before} placeholder IDs converted to Entrez IDs")
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
        import traceback
        logger.error(traceback.format_exc())
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

log_message "Testing improved placeholder ID fix"

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

log_message "Running improved fix on preprocessed ENCODE dataset"
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

log_message "Fix testing complete. Check the log file for results: $LOG_FILE"