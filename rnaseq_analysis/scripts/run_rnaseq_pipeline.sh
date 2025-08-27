#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/run_rnaseq_pipeline.sh
#!/bin/bash
# Run the complete RNA-seq standardization pipeline

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR="${BASE_DIR}/scripts"
METADATA_DIR="${BASE_DIR}/metadata/json"
OUTPUT_DIR="${BASE_DIR}/standardized_data"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Process datasets with Stage 1 (standardize_datasets.py)
echo "=== Stage 1: Initial Data Conversion ==="
python "${SCRIPTS_DIR}/standardize_datasets.py" \
    --encode-dir "${BASE_DIR}/encode/raw_data" \
    --encode-entex-dir "${BASE_DIR}/encode/entex" \
    --entex-metadata-file "${BASE_DIR}/encode/metadata/entex_metadata.json" \
    --gtex-file "${BASE_DIR}/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz" \
    --mage-dir "${BASE_DIR}/mage" \
    --adni-dir "${BASE_DIR}/adni_microarray" \
    --metadata-dir "${METADATA_DIR}" \
    --output-dir "$OUTPUT_DIR"

# Check if Stage 1 ran successfully
if [ $? -eq 0 ]; then
    echo "Stage 1 completed successfully!"
else
    echo "Stage 1 failed. Check the log file for details."
    exit 1
fi

# Process datasets with Stage 2 (standardize_metadata.py)
echo "=== Stage 2: Enhanced Metadata Standardization ==="
python "${SCRIPTS_DIR}/standardize_metadata.py" \
    --data-dir "$OUTPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --metadata-dir "${METADATA_DIR}"

# Check if Stage 2 ran successfully
if [ $? -eq 0 ]; then
    echo "Stage 2 completed successfully!"
else
    echo "Stage 2 failed. Check the log file for details."
    exit 1
fi

echo "Complete pipeline executed successfully!"
echo "Results saved to $OUTPUT_DIR"
echo "See the log files for details."