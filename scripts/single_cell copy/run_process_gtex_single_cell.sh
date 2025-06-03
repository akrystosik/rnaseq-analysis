#!/bin/bash
# Script to process GTEx single-cell RNA-seq data for pseudobulking.

# --- Configuration ---
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_PIPELINE_DIR="${BASE_DIR}/scripts/pipeline"
SCRIPTS_SINGLE_CELL_DIR="${BASE_DIR}/scripts/single_cell"

# Input GTEx single-cell H5AD file
INPUT_H5AD_SC="${BASE_DIR}/gtex/raw_data/snRNAseq_atlas/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad" 

METADATA_JSON_DIR="${BASE_DIR}/metadata/json" # Ensure gtex_scrnaseq_metadata.json is here

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR_BASE="${BASE_DIR}/standardized_data/gtex_single_cell_runs"
OUTPUT_DIR="${OUTPUT_DIR_BASE}/run_${TIMESTAMP}"
LOG_DIR="${BASE_DIR}/logs/single_cell"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/process_gtex_sc_${TIMESTAMP}.log"

# --- Helper Functions ---
log_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

run_command() {
    log_message "Running: $1"
    # Set PYTHONPATH to include the pipeline scripts directory for imports
    PYTHONPATH="${SCRIPTS_PIPELINE_DIR}:${PYTHONPATH}" eval "$1" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

# --- Pre-flight Checks ---
log_message "=== Starting GTEx Single-Cell Processing ==="
log_message "Output directory: $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"

if [ ! -f "$INPUT_H5AD_SC" ]; then
    log_message "ERROR: Input GTEx single-cell H5AD file not found at $INPUT_H5AD_SC"
    log_message "Please download it from https://storage.googleapis.com/adult-gtex/single-cell/v9/snrna-seq-data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad and update INPUT_H5AD_SC path in this script."
    exit 1
fi

if [ ! -d "$METADATA_JSON_DIR" ]; then
    log_message "ERROR: Metadata JSON directory not found at $METADATA_JSON_DIR"
    exit 1
fi

if [ ! -f "${METADATA_JSON_DIR}/gtex_scrnaseq_metadata.json" ]; then
    log_message "ERROR: gtex_scrnaseq_metadata.json not found in $METADATA_JSON_DIR"
    log_message "Please create this file (see previous instructions)."
    exit 1
fi

if [ ! -f "${SCRIPTS_SINGLE_CELL_DIR}/process_gtex_single_cell.py" ]; then
    log_message "ERROR: Python script process_gtex_single_cell.py not found in $SCRIPTS_SINGLE_CELL_DIR"
    exit 1
fi

if [ ! -d "$SCRIPTS_PIPELINE_DIR" ]; then
    log_message "ERROR: Pipeline scripts directory for utils ($SCRIPTS_PIPELINE_DIR) not found."
    exit 1
fi

# --- Main Processing Step ---
OUTPUT_H5AD_PB="${OUTPUT_DIR}/gtex_scrnaseq_pseudobulk_standardized.h5ad"

COMMAND="python \"${SCRIPTS_SINGLE_CELL_DIR}/process_gtex_single_cell.py\" \
    --input-h5ad \"$INPUT_H5AD_SC\" \
    --output-h5ad \"$OUTPUT_H5AD_PB\" \
    --metadata-json-dir \"$METADATA_JSON_DIR\""

run_command "$COMMAND"

if [ $? -eq 0 ]; then
    log_message "GTEx single-cell processing completed successfully!"
    log_message "Pseudobulk output saved to: $OUTPUT_H5AD_PB"
else
    log_message "GTEx single-cell processing failed. Check log for details: $LOG_FILE"
    exit 1
fi

# Create a symbolic link to the latest run
LATEST_LINK="${OUTPUT_DIR_BASE}/latest"
log_message "Creating symbolic link to latest run: $LATEST_LINK"
ln -sfn "$OUTPUT_DIR" "$LATEST_LINK"

log_message "=== GTEx Single-Cell Processing Finished ==="