#!/bin/bash
# Quick test version of RNA-seq pipeline that excludes GTEx data

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR="${BASE_DIR}/scripts"
METADATA_DIR="${BASE_DIR}/metadata/json"

# Create a timestamp for output directory and log file
TIMESTAMP="QUICK_$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="${BASE_DIR}/standardized_data/quick_test_${TIMESTAMP}"
PREPROCESSED_DIR="${BASE_DIR}/preprocessed_data/quick_test_${TIMESTAMP}"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/quick_test_${TIMESTAMP}.log"

# Create directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$PREPROCESSED_DIR"
mkdir -p "$METADATA_DIR"

# Function to log messages
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# Function to run a command with logging
run_command() {
    log_message "Running: $1"
    eval "$1" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

# Log the start of the pipeline with versioned directories
log_message "=== Starting Quick Test RNA-seq Pipeline (No GTEx) ==="
log_message "Output directory: $OUTPUT_DIR"
log_message "Preprocessed directory: $PREPROCESSED_DIR"
log_message "Log file: $LOG_FILE"

# Step 0: Generate Entrez to Ensembl mapping
log_message "=== Stage 0: Generating Entrez to Ensembl Mapping ==="
ENTREZ_MAPPING="${METADATA_DIR}/entrez_to_ensembl_mapping.csv"

# Check if --force flag was passed
FORCE_FLAG=""
if [ "$1" == "--force" ] || [ "$1" == "--force-mapping" ]; then
    FORCE_FLAG="--force"
    log_message "Force flag detected. Will regenerate mapping files."
fi

run_command "python ${SCRIPTS_DIR}/entrez-to-ensembl-mapping.py --output ${ENTREZ_MAPPING} --species human ${FORCE_FLAG}"

# Step 1: Process only ENCODE, MAGE, and ADNI
log_message "=== Stage 1: Processing Selected Datasets Only ==="

# Process ENCODE
log_message "Processing ENCODE data..."
run_command "python ${SCRIPTS_DIR}/standardize_datasets.py \\
    --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
    --metadata-dir \"${METADATA_DIR}\" \\
    --output-dir \"$OUTPUT_DIR\""

# Process MAGE
log_message "Processing MAGE data..."
run_command "python ${SCRIPTS_DIR}/standardize_datasets.py \\
    --mage-dir \"${BASE_DIR}/mage\" \\
    --metadata-dir \"${METADATA_DIR}\" \\
    --output-dir \"$OUTPUT_DIR\""

# Process ADNI
log_message "Processing ADNI data..."
run_command "python ${SCRIPTS_DIR}/standardize_datasets.py \\
    --adni-dir \"${BASE_DIR}/adni_microarray\" \\
    --metadata-dir \"${METADATA_DIR}\" \\
    --output-dir \"$OUTPUT_DIR\""

# Step 2: Generate Gene ID Reference Mapping
log_message "=== Stage 2: Generating Gene ID Reference Mapping ==="
run_command "python ${SCRIPTS_DIR}/gene_id_mapping_reference.py \\
    --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
    --entex-dir \"${BASE_DIR}/encode/entex\" \\
    --entrez-mapping \"${ENTREZ_MAPPING}\" \\
    --output \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
    ${FORCE_FLAG}"

# Step 3: Generate ENCODE ID to Ensembl mapping
log_message "=== Stage 3: Generating ENCODE ID to Ensembl Mapping ==="
run_command "python ${SCRIPTS_DIR}/generate_encode_mapping.py \\
    --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
    --output-dir \"${METADATA_DIR}/gene_mapping\" \\
    ${FORCE_FLAG}"

# Step 4: Preprocess Datasets for Consistent Gene IDs
log_message "=== Stage 4: Preprocessing Datasets for Consistent Gene IDs ==="
run_command "python ${SCRIPTS_DIR}/preprocess_dataset_gene_ids.py \\
    --data-dir \"$OUTPUT_DIR\" \\
    --reference-mapping \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
    --output-dir \"$PREPROCESSED_DIR\" \\
    --datasets \"encode,mage,adni\" \\
    ${FORCE_FLAG}"

# Step 5: Create combined dataset with ALL genes (new approach with improved gene ID mapping)
log_message "=== Stage 5: Creating Combined Dataset with All Genes ==="
run_command "python ${SCRIPTS_DIR}/create_combined_dataset_all_genes_sparse.py \\
    --input-dir \"$PREPROCESSED_DIR\" \\
    --reference-mapping \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
    --output-file \"${OUTPUT_DIR}/combined_all_genes_sparse_standardized.h5ad\" \\
    --include-datasets \"encode,mage,adni\" \\
    ${FORCE_FLAG}"

# Step 6: Run validation on standardized datasets
log_message "=== Stage 6: Validating Standardized Datasets ==="
VALIDATION_OUTPUT="${OUTPUT_DIR}/validation_report_${TIMESTAMP}.json"
run_command "python ${SCRIPTS_DIR}/validate_standardized_datasets.py \\
    --input-dir \"$OUTPUT_DIR\" \\
    --output-file \"$VALIDATION_OUTPUT\" \\
    --file-pattern \"*_standardized*.h5ad\""

# Also create a symbolic link to the latest quick test
log_message "=== Creating symbolic link to latest quick test ==="
LATEST_LINK="${BASE_DIR}/standardized_data/latest_quick_test"
run_command "ln -sfn \"$OUTPUT_DIR\" \"$LATEST_LINK\""

log_message "Quick test pipeline executed successfully!"
log_message "Results saved to $OUTPUT_DIR"
log_message "Symbolic link to latest quick test: $LATEST_LINK"
log_message "Log file: $LOG_FILE"
