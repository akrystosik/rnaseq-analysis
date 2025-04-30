#!/bin/bash
# Test script to run just the combined dataset creation and validation steps

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR="${BASE_DIR}/scripts"
OUTPUT_DIR="${BASE_DIR}/standardized_data"

echo "=== TEST: Creating Combined Dataset with ALL Genes ==="
# Create a timestamped log file
LOG_FILE="${OUTPUT_DIR}/combined_all_genes_$(date '+%Y%m%d_%H%M%S').log"
echo "Starting test at $(date)" > $LOG_FILE
echo "Running create_combined_dataset_all_genes.py" >> $LOG_FILE

# Run the script and capture output to the log file
python "${SCRIPTS_DIR}/create_combined_dataset_all_genes.py" 2>&1 | tee -a $LOG_FILE

# Check if all-genes combined dataset creation was successful
if [ $? -eq 0 ]; then
    echo "Combined dataset (all genes) created successfully!" | tee -a $LOG_FILE
else
    echo "Combined dataset (all genes) creation failed. Check the log file for details." | tee -a $LOG_FILE
    # Continue execution even if this fails
fi

# Validate the new combined dataset
echo "=== TEST: Validating the Combined Dataset with ALL Genes ===" | tee -a $LOG_FILE
VALIDATION_OUTPUT="${OUTPUT_DIR}/validation_combined_all_genes_$(date '+%Y%m%d_%H%M%S').json"
python "${SCRIPTS_DIR}/validate_standardized_datasets.py" \
    --input-dir "$OUTPUT_DIR" \
    --output-file "$VALIDATION_OUTPUT" \
    --file-pattern "combined_all_genes_standardized.h5ad" 2>&1 | tee -a $LOG_FILE

# Check if validation ran successfully
if [ $? -eq 0 ]; then
    echo "Validation completed successfully!" | tee -a $LOG_FILE
    echo "Validation report saved to $VALIDATION_OUTPUT" | tee -a $LOG_FILE
else
    echo "Validation failed. Check the log file for details." | tee -a $LOG_FILE
fi

echo "Test completed at $(date)" | tee -a $LOG_FILE
echo "Log file saved to: $LOG_FILE"