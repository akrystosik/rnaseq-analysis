#!/bin/bash
# Test pipeline script with fixed components

# Set base directory and output directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR="${BASE_DIR}/scripts"
METADATA_DIR="${BASE_DIR}/metadata/json"

# Create a timestamp for output directory and log file
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="${BASE_DIR}/standardized_data/test_run_${TIMESTAMP}"
LOG_FILE="${BASE_DIR}/logs/test_pipeline_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname $LOG_FILE)"

# Function to log messages
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# Function to run a command
run_command() {
    log_message "Running: $1"
    eval "$1" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

# Log the start of the pipeline
log_message "=== Starting RNA-seq Pipeline Test ==="
log_message "Output directory: $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"

# Step 1: Process ADNI dataset (previously failing)
log_message "=== Step 1: Processing ADNI dataset ==="
run_command "python ${SCRIPTS_DIR}/standardize_datasets.py \\
    --adni-dir \"${BASE_DIR}/adni_microarray\" \\
    --metadata-dir \"${METADATA_DIR}\" \\
    --output-dir \"$OUTPUT_DIR\""

# Step 2: Fix categorical columns in all datasets
log_message "=== Step 2: Fixing categorical columns ==="
run_command "python ${SCRIPTS_DIR}/fix_categorical_columns.py \\
    --input \"$OUTPUT_DIR\" \\
    --output \"$OUTPUT_DIR\" \\
    --pattern \"*.h5ad\""

# Step 3: Generate gene ID reference mapping
log_message "=== Step 3: Generating gene ID reference mapping ==="
run_command "python ${SCRIPTS_DIR}/gene_id_mapping_reference.py \\
    --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
    --entex-dir \"${BASE_DIR}/encode/entex\" \\
    --entrez-mapping \"${METADATA_DIR}/entrez_to_ensembl_mapping.csv\" \\
    --output \"${OUTPUT_DIR}/gene_id_reference_mapping.csv\" \\
    --force"

# Step 4: Preprocess dataset gene IDs
log_message "=== Step 4: Preprocessing dataset gene IDs ==="
run_command "python ${SCRIPTS_DIR}/preprocess_dataset_gene_ids.py \\
    --data-dir \"$OUTPUT_DIR\" \\
    --reference-mapping \"${OUTPUT_DIR}/gene_id_reference_mapping.csv\" \\
    --output-dir \"${OUTPUT_DIR}/preprocessed_data\" \\
    --datasets \"adni\" \\
    --force"

# Step 5: Fix placeholder IDs
log_message "=== Step 5: Fixing placeholder IDs ==="
preprocessed_file="${OUTPUT_DIR}/preprocessed_data/adni_standardized_preprocessed.h5ad"
if [ -f "$preprocessed_file" ]; then
    run_command "python ${SCRIPTS_DIR}/fix_placeholder_ids.py \\
        \"$preprocessed_file\" \\
        \"${preprocessed_file}.fixed\""
    
    if [ $? -eq 0 ]; then
        mv "${preprocessed_file}.fixed" "$preprocessed_file"
    fi
fi

# Step 6: Apply metadata from JSON files
log_message "=== Step 6: Applying metadata from JSON files ==="
run_command "python ${SCRIPTS_DIR}/apply_dataset_metadata.py \\
    --input \"${OUTPUT_DIR}/preprocessed_data\" \\
    --output \"${OUTPUT_DIR}/preprocessed_data\" \\
    --metadata-dir \"${METADATA_DIR}\" \\
    --pattern \"*.h5ad\""

# Step 7: Analyze results
log_message "=== Step 7: Analyzing results ==="
run_command "python - <<EOF
import scanpy as sc
import pandas as pd
import numpy as np

print('\\n===== TESTING ADNI PREPROCESSED DATASET =====')
try:
    adata = sc.read_h5ad('${OUTPUT_DIR}/preprocessed_data/adni_standardized_preprocessed.h5ad')
    print(f'Dataset loaded: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Basic dataset info
    print('\\nBasic information:')
    print(f'  Data type: {adata.uns.get(\"dataset_info\", {}).get(\"data_type\", \"Unknown\")}')
    print(f'  Expression unit: {adata.uns.get(\"dataset_info\", {}).get(\"expression_unit\", \"Unknown\")}')
    print(f'  GENCODE version: {adata.uns.get(\"gencode_version\", \"Unknown\")}')
    print(f'  Reference genome: {adata.uns.get(\"reference_genome\", \"Unknown\")}')
    
    # Check tissue information
    if 'tissue' in adata.obs.columns:
        missing_tissue = pd.isna(adata.obs['tissue']).sum()
        nan_string = sum(adata.obs['tissue'].astype(str) == 'nan')
        print(f'\\nTissue information:')
        print(f'  Missing tissues (NaN): {missing_tissue}/{adata.n_obs} ({missing_tissue/adata.n_obs*100:.2f}%)')
        print(f'  String \"nan\" tissues: {nan_string}/{adata.n_obs} ({nan_string/adata.n_obs*100:.2f}%)')
        
        # Show top tissues
        if adata.n_obs > 0:
            tissue_counts = adata.obs['tissue'].value_counts().head(5)
            print('\\n  Top 5 tissues:')
            for tissue, count in tissue_counts.items():
                print(f'    {tissue}: {count}')
    
    # Check gene ID information
    if 'gene_id' in adata.var.columns and 'ensembl_id' in adata.var.columns:
        # Check ID formats
        ensembl_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('ENSG'))
        entrez_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('ENTREZ:'))
        placeholder_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('PLACEHOLDER_'))
        
        print('\\nGene ID analysis:')
        print(f'  Ensembl IDs: {ensembl_count}/{adata.n_vars} ({ensembl_count/adata.n_vars*100:.2f}%)')
        print(f'  Entrez IDs: {entrez_count}/{adata.n_vars} ({entrez_count/adata.n_vars*100:.2f}%)')
        print(f'  Placeholder IDs: {placeholder_count}/{adata.n_vars} ({placeholder_count/adata.n_vars*100:.2f}%)')
        
        # Check if gene_id matches ensembl_id (convert both to string first)
        matching_ids = sum(adata.var['gene_id'].astype(str) == adata.var['ensembl_id'].astype(str))
        print(f'  gene_id matches ensembl_id: {matching_ids}/{adata.n_vars} ({matching_ids/adata.n_vars*100:.2f}%)')
        
        # Check original gene IDs with version
        if 'original_gene_id' in adata.var.columns:
            version_pattern = r'ENSG\\d+\\.\\d+'
            has_version = adata.var['original_gene_id'].astype(str).str.match(version_pattern)
            original_versions = sum(has_version)
            print(f'  Original IDs with version numbers: {original_versions}/{adata.n_vars} ({original_versions/adata.n_vars*100:.2f}%)')
    
    # Check required fields
    required_obs = ['dataset', 'sample_id', 'data_type', 'expression_unit', 'tissue', 'assay_ontology', 'species', 'species_ontology']
    required_var = ['gene_id', 'ensembl_id', 'gene_name', 'gene_type', 'chromosome', 'mapping_source']
    
    print('\\nRequired fields check:')
    for field in required_obs:
        if field in adata.obs.columns:
            print(f'  obs[{field}]: Present')
        else:
            print(f'  obs[{field}]: MISSING')
    
    for field in required_var:
        if field in adata.var.columns:
            print(f'  var[{field}]: Present')
        else:
            print(f'  var[{field}]: MISSING')
    
    # Check metadata source
    print('\\nMetadata sources:')
    if 'dataset_info' in adata.uns:
        for key, value in adata.uns['dataset_info'].items():
            print(f'  {key}: {value}')
except Exception as e:
    print(f'Error analyzing dataset: {e}')
    import traceback
    print(traceback.format_exc())
EOF"

log_message "=== Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"
