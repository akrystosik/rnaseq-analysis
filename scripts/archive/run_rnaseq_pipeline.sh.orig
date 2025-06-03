#!/bin/bash
# Run the complete RNA-seq standardization pipeline with improved gene ID mapping
# Modified to add date suffix to output directory for version control

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR="${BASE_DIR}/scripts"
METADATA_DIR="${BASE_DIR}/metadata/json"

# Create a timestamp for output directory and log file
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="${BASE_DIR}/standardized_data/run_${TIMESTAMP}"
PREPROCESSED_DIR="${BASE_DIR}/preprocessed_data/run_${TIMESTAMP}"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/pipeline_${TIMESTAMP}.log"

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
log_message "=== Starting RNA-seq Pipeline with Versioned Output ==="
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

# Check if mapping generation was successful
if [ $? -eq 0 ]; then
    log_message "Entrez to Ensembl mapping generated or verified successfully!"
else
    log_message "Warning: Entrez to Ensembl mapping generation failed. Check the log file for details."
    log_message "Will continue without this mapping, which may affect ENCODE gene ID mapping quality."
fi

# Step 1: Initial Data Conversion (existing)
log_message "=== Stage 1: Initial Data Conversion ==="
run_command "python ${SCRIPTS_DIR}/standardize_datasets.py \\
    --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
    --encode-entex-dir \"${BASE_DIR}/encode/entex\" \\
    --entex-metadata-file \"${BASE_DIR}/encode/metadata/entex_metadata.json\" \\
    --gtex-file \"${BASE_DIR}/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz\" \\
    --mage-dir \"${BASE_DIR}/mage\" \\
    --adni-dir \"${BASE_DIR}/adni_microarray\" \\
    --metadata-dir \"${METADATA_DIR}\" \\
    --output-dir \"$OUTPUT_DIR\""

# Check if Stage 1 ran successfully
if [ $? -eq 0 ]; then
    log_message "Stage 1 completed successfully!"
else
    log_message "Stage 1 failed. Check the log file for details."
    exit 1
fi

# Step 1.5: Generate Gene ID Reference Mapping with Entrez mapping
log_message "=== Stage 1.5: Generating Gene ID Reference Mapping ==="

# Check if Entrez mapping exists and add it to the command if it does
if [ -f "$ENTREZ_MAPPING" ]; then
    log_message "Using Entrez to Ensembl mapping for enhanced gene ID mapping"
    run_command "python ${SCRIPTS_DIR}/gene_id_mapping_reference.py \\
        --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
        --entex-dir \"${BASE_DIR}/encode/entex\" \\
        --entrez-mapping \"${ENTREZ_MAPPING}\" \\
        --output \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
        ${FORCE_FLAG}"
else
    log_message "Warning: Entrez to Ensembl mapping not found. Proceeding without it."
    run_command "python ${SCRIPTS_DIR}/gene_id_mapping_reference.py \\
        --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
        --entex-dir \"${BASE_DIR}/encode/entex\" \\
        --output \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
        ${FORCE_FLAG}"
fi

# Check if Gene ID Reference Mapping ran successfully
if [ $? -eq 0 ]; then
    log_message "Gene ID Reference Mapping completed successfully!"
else
    log_message "Gene ID Reference Mapping failed. Check the log file for details."
    exit 1
fi

# Step 1.6: Generate ENCODE ID to Ensembl mapping
log_message "=== Stage 1.6: Generating ENCODE ID to Ensembl Mapping ==="
run_command "python ${SCRIPTS_DIR}/generate_encode_mapping.py \\
    --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
    --output-dir \"${METADATA_DIR}/gene_mapping\" \\
    ${FORCE_FLAG}"

# Check if ENCODE ID mapping generation was successful
if [ $? -eq 0 ]; then
    log_message "ENCODE ID mapping generated successfully!"
else
    log_message "Warning: ENCODE ID mapping generation failed. Check the log file for details."
    log_message "Will continue without this mapping, which may affect ENCODE gene ID mapping quality."
fi

# Step 2: Enhanced Metadata Standardization
log_message "=== Stage 2: Enhanced Metadata Standardization ==="
run_command "python ${SCRIPTS_DIR}/standardize_metadata.py \\
    --data-dir \"$OUTPUT_DIR\" \\
    --output-dir \"$OUTPUT_DIR\" \\
    --metadata-dir \"${METADATA_DIR}\""

# Check if Stage 2 ran successfully
if [ $? -eq 0 ]; then
    log_message "Stage 2 completed successfully!"
else
    log_message "Stage 2 failed. Check the log file for details."
    exit 1
fi

# Step 2.5: Preprocess Datasets for Consistent Gene IDs
log_message "=== Stage 2.5: Preprocessing Datasets for Consistent Gene IDs ==="
run_command "python ${SCRIPTS_DIR}/preprocess_dataset_gene_ids.py \\
    --data-dir \"$OUTPUT_DIR\" \\
    --reference-mapping \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
    --output-dir \"$PREPROCESSED_DIR\" \\
    ${FORCE_FLAG}"

# Check if Preprocessing ran successfully
if [ $? -eq 0 ]; then
    log_message "Dataset Preprocessing completed successfully!"
else
    log_message "Dataset Preprocessing failed. Check the log file for details."
    exit 1
fi

# Step 2.6: Analyze ENCODE mapping quality
log_message "=== Stage 2.6: Analyzing ENCODE Gene Mapping Quality ==="
run_command "python - <<EOF
import scanpy as sc
import pandas as pd
import numpy as np

# Load the ENCODE preprocessed dataset
encode_path = '${PREPROCESSED_DIR}/encode_standardized_preprocessed.h5ad'
try:
    adata = sc.read_h5ad(encode_path)
    print(f\"ENCODE dataset shape: {adata.shape}\")

    # Check ensembl_id column
    if 'ensembl_id' in adata.var.columns:
        # Count non-empty ensembl_ids
        non_empty = sum(1 for x in adata.var['ensembl_id'] if x and str(x).strip() != '')
        percentage = non_empty / len(adata.var) * 100
        print(f\"ENCODE genes with mapped Ensembl IDs: {non_empty}/{len(adata.var)} ({percentage:.2f}%)\")
        
        # Sample of mapped genes
        mapped_ids = adata.var.loc[adata.var['ensembl_id'] != '', 'ensembl_id'].iloc[:5].tolist()
        print(f\"Sample mapped Ensembl IDs: {mapped_ids}\")
        
        # Mapping source distribution
        if 'mapping_source' in adata.var.columns:
            source_counts = adata.var['mapping_source'].value_counts()
            for source, count in source_counts.items():
                source_percentage = count / len(adata.var) * 100
                print(f\"Mapping source '{source}': {count} ({source_percentage:.2f}%)\")
    else:
        print(\"Error: ENCODE dataset does not have an 'ensembl_id' column!\")
    
    # Save mapping stats to file
    mapping_stats = {
        'total_genes': len(adata.var),
        'mapped_genes': non_empty,
        'mapping_percentage': percentage,
        'mapping_sources': {k: int(v) for k, v in source_counts.items()} if 'source_counts' in locals() else {}
    }
    
    import json
    with open('${OUTPUT_DIR}/encode_mapping_stats.json', 'w') as f:
        json.dump(mapping_stats, f, indent=2)
    
    print(f\"ENCODE mapping stats saved to ${OUTPUT_DIR}/encode_mapping_stats.json\")
    
except Exception as e:
    print(f\"Error analyzing ENCODE dataset: {e}\")
    print(\"Continuing with pipeline execution...\")
EOF"

# Step 3a: Create combined dataset with common genes (original approach)
log_message "=== Stage 3a: Creating Combined Dataset with Common Genes ==="
run_command "python ${SCRIPTS_DIR}/create_combined_dataset.py \\
    --output-file \"${OUTPUT_DIR}/combined_common_genes_standardized.h5ad\""

# Check if combined dataset creation was successful
if [ $? -eq 0 ]; then
    log_message "Combined dataset (common genes) created successfully!"
else
    log_message "Combined dataset (common genes) creation failed. Check the log file for details."
    # Continue execution even if this fails
fi

# Step 3b: Create combined dataset with ALL genes (new approach with improved gene ID mapping)
log_message "=== Stage 3b: Creating Sparse Combined Dataset with ALL Genes ==="
run_command "python ${SCRIPTS_DIR}/create_combined_dataset_all_genes_sparse.py \\
    --input-dir \"$PREPROCESSED_DIR\" \\
    --reference-mapping \"${METADATA_DIR}/gene_id_reference_mapping.csv\" \\
    --output-file \"${OUTPUT_DIR}/combined_all_genes_sparse_standardized.h5ad\" \\
    --include-datasets \"encode,gtex,mage,adni,entex\" \\
    ${FORCE_FLAG}"

# Check if all-genes combined dataset creation was successful
if [ $? -eq 0 ]; then
    log_message "Combined dataset (all genes) created successfully!"
else
    log_message "Combined dataset (all genes) creation failed. Check the log file for details."
    # Continue execution even if this fails
fi

# Step 3c: Analyze the combined dataset
log_message "=== Stage 3c: Analyzing Combined Dataset ==="
run_command "python - <<EOF
import scanpy as sc
import pandas as pd
import numpy as np

# Load the combined dataset
combined_path = '${OUTPUT_DIR}/combined_all_genes_sparse_standardized.h5ad'
try:
    adata = sc.read_h5ad(combined_path)
    print(f\"Combined dataset shape: {adata.shape}\")

    # Check dataset sample distribution
    print(\"\\nSample distribution by dataset:\")
    dataset_counts = adata.obs['dataset'].value_counts()
    for dataset, count in dataset_counts.items():
        percentage = count / len(adata.obs) * 100
        print(f\"  - {dataset}: {count} ({percentage:.2f}%)\")

    # Check gene presence by dataset
    print(\"\\nGene presence by dataset:\")
    dataset_gene_counts = {}
    for dataset in ['adni', 'encode', 'entex', 'gtex', 'mage']:
        count = sum(1 for x in adata.var['present_in_datasets'] if dataset in str(x))
        percentage = count / len(adata.var) * 100
        print(f\"  - {dataset}: {count}/{len(adata.var)} ({percentage:.2f}%)\")
        dataset_gene_counts[dataset] = count

    # Check ensembl ID vs. other ID format distribution
    ensembl_count = sum(1 for x in adata.var_names if str(x).startswith('ENSG'))
    spike_in_count = sum(1 for x in adata.var_names if str(x).startswith('gSpikein'))
    other_count = len(adata.var) - ensembl_count - spike_in_count

    print(\"\\nGene ID format distribution:\")
    print(f\"  - Ensembl IDs: {ensembl_count}/{len(adata.var)} ({ensembl_count/len(adata.var)*100:.2f}%)\")
    print(f\"  - Spike-in controls: {spike_in_count}/{len(adata.var)} ({spike_in_count/len(adata.var)*100:.2f}%)\")
    print(f\"  - Other IDs: {other_count}/{len(adata.var)} ({other_count/len(adata.var)*100:.2f}%)\")
    
    # Save stats to file
    combined_stats = {
        'total_samples': len(adata.obs),
        'total_genes': len(adata.var),
        'dataset_distribution': {k: int(v) for k, v in dataset_counts.items()},
        'gene_presence': dataset_gene_counts,
        'gene_id_distribution': {
            'ensembl': int(ensembl_count),
            'spike_in': int(spike_in_count),
            'other': int(other_count)
        }
    }
    
    import json
    with open('${OUTPUT_DIR}/combined_dataset_stats.json', 'w') as f:
        json.dump(combined_stats, f, indent=2)
    
    print(f\"Combined dataset stats saved to ${OUTPUT_DIR}/combined_dataset_stats.json\")
    
except Exception as e:
    print(f\"Error analyzing combined dataset: {e}\")
    print(\"Continuing with pipeline execution...\")
EOF"

# Step 4: Run validation on all standardized datasets
log_message "=== Stage 4: Validating Standardized Datasets ==="
VALIDATION_OUTPUT="${OUTPUT_DIR}/validation_report_${TIMESTAMP}.json"
run_command "python ${SCRIPTS_DIR}/validate_standardized_datasets.py \\
    --input-dir \"$OUTPUT_DIR\" \\
    --output-file \"$VALIDATION_OUTPUT\" \\
    --file-pattern \"*_standardized*.h5ad\""

# Check if validation ran successfully
if [ $? -eq 0 ]; then
    log_message "Validation completed successfully!"
    log_message "Validation report saved to $VALIDATION_OUTPUT"
else
    log_message "Validation failed. Check the log file for details."
    # Continue execution even if this fails
fi

# Also create a symbolic link to the latest run
log_message "=== Creating symbolic link to latest run ==="
LATEST_LINK="${BASE_DIR}/standardized_data/latest"
run_command "ln -sfn \"$OUTPUT_DIR\" \"$LATEST_LINK\""

log_message "Complete pipeline executed successfully!"
log_message "Results saved to $OUTPUT_DIR"
log_message "Symbolic link to latest run: $LATEST_LINK"
log_message "Log file: $LOG_FILE"

# Print a summary of the ENCODE mapping improvement
log_message "=== ENCODE Mapping Summary ==="
if [ -f "${OUTPUT_DIR}/encode_mapping_stats.json" ]; then
    # Extract mapping percentage from the stats file
    MAPPING_PERCENTAGE=$(grep -o '"mapping_percentage":[^,}]*' "${OUTPUT_DIR}/encode_mapping_stats.json" | cut -d':' -f2)
    log_message "ENCODE gene mapping percentage: ${MAPPING_PERCENTAGE}%"
fi

if [ -f "${OUTPUT_DIR}/combined_dataset_stats.json" ]; then
    # Extract ENCODE gene presence
    ENCODE_GENES=$(grep -o '"encode":[^,}]*' "${OUTPUT_DIR}/combined_dataset_stats.json" | head -1 | cut -d':' -f2)
    TOTAL_GENES=$(grep -o '"total_genes":[^,}]*' "${OUTPUT_DIR}/combined_dataset_stats.json" | head -1 | cut -d':' -f2)
    
    if [ ! -z "$ENCODE_GENES" ] && [ ! -z "$TOTAL_GENES" ]; then
        ENCODE_PERCENTAGE=$(echo "scale=2; 100 * $ENCODE_GENES / $TOTAL_GENES" | bc)
        log_message "ENCODE genes in combined dataset: ${ENCODE_GENES}/${TOTAL_GENES} (${ENCODE_PERCENTAGE}%)"
    fi
fi

log_message "Run the pipeline again with --force flag to regenerate all files"