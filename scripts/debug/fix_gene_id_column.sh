#!/bin/bash
# Fix script for gene_id column to use proper Ensembl IDs

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1"
}

log_message "Creating fix for gene_id column in preprocessed datasets"

# Create a Python script to fix the preprocess_dataset_gene_ids.py file
cat > /tmp/fix_gene_id_column.py << 'EOL'
import re

# Path to the preprocess_dataset_gene_ids.py file
file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/preprocess_dataset_gene_ids.py'

# Read the file
with open(file_path, 'r') as f:
    content = f.read()

# Find the section where the var DataFrame is created in preprocess_encode_dataset
# We want to modify how gene_id is populated, using ensembl_id instead of the index
pattern = r"""    # Create new var DataFrame
    new_var = pd.DataFrame\(var_columns\)"""

# Replace with modified code that sets the index to ensembl_id and uses it for gene_id
replacement = """    # Create new var DataFrame
    new_var = pd.DataFrame(var_columns)
    
    # Fix gene_id column to use ensembl_id (or ENTREZ: id) instead of sequential index
    for i, row in new_var.iterrows():
        if row['ensembl_id']:
            # Use the ensembl_id as the gene_id
            new_var.at[i, 'gene_id'] = row['ensembl_id']
        elif row['original_gene_id'].startswith('gSpikein'):
            # For spike-in controls, use the original ID
            new_var.at[i, 'gene_id'] = row['original_gene_id']"""

# Apply the replacement
content = re.sub(pattern, replacement, content)

# Write the updated file
with open(file_path, 'w') as f:
    f.write(content)

print("gene_id column fix applied to preprocess_dataset_gene_ids.py")
EOL

python /tmp/fix_gene_id_column.py

log_message "Testing the gene_id column fix"

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_gene_id_fix_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/gene_id_fix_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR/input"
mkdir -p "$OUTPUT_DIR/output"
mkdir -p "$(dirname $LOG_FILE)"

# Copy the original dataset for testing
cp "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_mini_pipeline_20250430_073015/encode_standardized_v1.h5ad" "$OUTPUT_DIR/input/"

# Run the preprocessing script with our fix
log_message "Running preprocess_dataset_gene_ids.py with the fix"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/preprocess_dataset_gene_ids.py \
    --data-dir "$OUTPUT_DIR/input" \
    --reference-mapping "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_id_reference_mapping.csv" \
    --output-dir "$OUTPUT_DIR/output" \
    --datasets encode 2>&1 | tee -a "$LOG_FILE"

# Analyze the results
log_message "Analyzing the results with fixed gene_id column"
cat > /tmp/analyze_gene_id_fix.py << 'EOL'
#!/usr/bin/env python3
"""
Analyze gene_id fix in preprocessed dataset.
"""
import scanpy as sc
import pandas as pd
import numpy as np
import sys

# Get the output directory from command line
output_dir = sys.argv[1]

# Path to preprocessed file
preprocessed_file = f"{output_dir}/output/encode_standardized_preprocessed.h5ad"

# Load the dataset
print(f"Loading preprocessed dataset from {preprocessed_file}")
try:
    adata = sc.read_h5ad(preprocessed_file)
    print(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
    
    # Examine gene_id and ensembl_id columns
    print("\nComparing gene_id and ensembl_id columns:")
    matching_ids = 0
    for i, (idx, row) in enumerate(adata.var.iterrows()):
        if i < 5:  # Print first 5 for inspection
            print(f"Index: {idx}, gene_id: {row.get('gene_id', 'N/A')}, ensembl_id: {row.get('ensembl_id', 'N/A')}")
        
        # Count matching values
        if str(row.get('gene_id', '')) == str(row.get('ensembl_id', '')):
            matching_ids += 1
    
    print(f"\nMatching gene_id and ensembl_id values: {matching_ids}/{len(adata.var)} ({matching_ids/len(adata.var)*100:.2f}%)")
    
    # Analyze gene_id values
    ensembl_count = sum(1 for id in adata.var['gene_id'] if str(id).startswith('ENSG'))
    entrez_count = sum(1 for id in adata.var['gene_id'] if str(id).startswith('ENTREZ:'))
    spike_count = sum(1 for id in adata.var['gene_id'] if str(id).startswith('gSpikein'))
    numeric_count = sum(1 for id in adata.var['gene_id'] if str(id).isdigit())
    
    print("\nGene ID distribution in gene_id column:")
    print(f"  Ensembl IDs: {ensembl_count}/{len(adata.var)} ({ensembl_count/len(adata.var)*100:.2f}%)")
    print(f"  Entrez IDs: {entrez_count}/{len(adata.var)} ({entrez_count/len(adata.var)*100:.2f}%)")
    print(f"  Spike-in IDs: {spike_count}/{len(adata.var)} ({spike_count/len(adata.var)*100:.2f}%)")
    print(f"  Numeric IDs: {numeric_count}/{len(adata.var)} ({numeric_count/len(adata.var)*100:.2f}%)")
    
    print("\nFix successfully applied!" if numeric_count == 0 else "\nFix not fully applied, still have numeric IDs")
except Exception as e:
    print(f"Error examining dataset: {e}")
    import traceback
    print(traceback.format_exc())
EOL

chmod +x /tmp/analyze_gene_id_fix.py
python /tmp/analyze_gene_id_fix.py "$OUTPUT_DIR" 2>&1 | tee -a "$LOG_FILE"

log_message "Fix testing complete. Check the log file for results: $LOG_FILE"