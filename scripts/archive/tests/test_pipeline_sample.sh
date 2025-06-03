#!/bin/bash
# test_pipeline_sample.sh
# Fast test for RNA-seq standardization pipeline using data samples

# Set base directory and create timestamp for output
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="${BASE_DIR}/test_results/sample_${TIMESTAMP}"
LOG_FILE="${BASE_DIR}/logs/test_sample_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname $LOG_FILE)"

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

log_message "=== Starting Sample Pipeline Test ==="
log_message "Output directory: $OUTPUT_DIR"

# Step 1: Create sample directories for each dataset
SAMPLE_DIR="${OUTPUT_DIR}/samples"
mkdir -p "${SAMPLE_DIR}"

# Create a function to sample each dataset
sample_dataset() {
    local dataset=$1
    local source_dir=$2
    local sample_size=$3
    local pattern=$4
    
    log_message "Sampling ${dataset} dataset (${sample_size} files)"
    mkdir -p "${SAMPLE_DIR}/${dataset}"
    
    # Find files matching pattern and select a sample
    find "${source_dir}" -name "${pattern}" | head -n ${sample_size} > "${SAMPLE_DIR}/${dataset}_files.txt"
    
    # Create symbolic links to sample files
    while read file; do
        ln -s "$file" "${SAMPLE_DIR}/${dataset}/$(basename $file)"
    done < "${SAMPLE_DIR}/${dataset}_files.txt"
    
    log_message "Created sample for ${dataset} with $(ls ${SAMPLE_DIR}/${dataset} | wc -l) files"
}

# Sample each dataset (adjust sample sizes as needed)
log_message "Step 1: Creating dataset samples"
sample_dataset "encode" "${BASE_DIR}/encode/raw_data/cell_lines" 3 "*.tsv"
sample_dataset "mage" "${BASE_DIR}/mage" 5 "*.csv"
sample_dataset "adni" "${BASE_DIR}/adni_microarray" 5 "*.csv"

# Step 2: Use existing entrez mapping to save time
log_message "Step 2: Using existing Entrez to Ensembl mapping"
cp "${BASE_DIR}/metadata/json/entrez_to_ensembl_mapping.csv" "${OUTPUT_DIR}/"

# Step 3: Run standardize_datasets.py with sample data
log_message "Step 3: Running standardize_datasets.py on samples"
python "${BASE_DIR}/scripts/standardize_datasets.py" \
    --encode-dir "${SAMPLE_DIR}/encode" \
    --mage-dir "${SAMPLE_DIR}/mage" \
    --adni-dir "${SAMPLE_DIR}/adni" \
    --metadata-dir "${BASE_DIR}/metadata/json" \
    --output-dir "${OUTPUT_DIR}" 2>&1 | tee -a "$LOG_FILE"

# Step 4: Run preprocess_dataset_gene_ids.py on standardized data
log_message "Step 4: Running preprocess_dataset_gene_ids.py"
mkdir -p "${OUTPUT_DIR}/preprocessed"
python "${BASE_DIR}/scripts/preprocess_dataset_gene_ids.py" \
    --data-dir "${OUTPUT_DIR}" \
    --reference-mapping "${BASE_DIR}/metadata/json/gene_id_reference_mapping.csv" \
    --output-dir "${OUTPUT_DIR}/preprocessed" \
    --datasets "encode,mage,adni" 2>&1 | tee -a "$LOG_FILE"

# Step 5: Validate the sample output
log_message "Step 5: Validating sample output"
python - <<EOF | tee -a "$LOG_FILE"
import scanpy as sc
import pandas as pd
import os
import glob

print("\n=== VALIDATION RESULTS ===")

# Check for existence of expected files
output_dir = "${OUTPUT_DIR}"
preprocessed_dir = "${OUTPUT_DIR}/preprocessed"

print("\nStandardized datasets:")
for dataset in ['encode', 'mage', 'adni']:
    file_path = os.path.join(output_dir, f"{dataset}_standardized_v1.h5ad")
    status = "✅ Found" if os.path.exists(file_path) else "❌ Missing"
    print(f"  {dataset}: {status}")

print("\nPreprocessed datasets:")
for dataset in ['encode', 'mage', 'adni']:
    file_path = os.path.join(preprocessed_dir, f"{dataset}_standardized_preprocessed.h5ad")
    status = "✅ Found" if os.path.exists(file_path) else "❌ Missing"
    print(f"  {dataset}: {status}")

# Analyze metadata and structure
print("\nValidating dataset structure:")
for dataset in ['encode', 'mage', 'adni']:
    try:
        file_path = os.path.join(preprocessed_dir, f"{dataset}_standardized_preprocessed.h5ad")
        if os.path.exists(file_path):
            adata = sc.read_h5ad(file_path)
            print(f"\n{dataset.upper()}:")
            print(f"  Samples: {adata.n_obs}, Genes: {adata.n_vars}")
            
            # Check required metadata fields
            required_obs = ['sample_id', 'subject_id', 'dataset', 'tissue']
            required_var = ['gene_id', 'original_gene_id', 'ensembl_id', 'gene_name', 'gene_type']
            
            print("  Observation metadata:")
            for field in required_obs:
                status = "✅" if field in adata.obs.columns else "❌"
                print(f"    {field}: {status}")
            
            print("  Variable metadata:")
            for field in required_var:
                status = "✅" if field in adata.var.columns else "❌"
                print(f"    {field}: {status}")
            
            # Check gene_id format
            if 'ensembl_id' in adata.var.columns:
                ensembl_count = sum(str(id).startswith('ENSG') for id in adata.var['ensembl_id'])
                print(f"  Gene ID format:")
                print(f"    Ensembl IDs: {ensembl_count}/{adata.n_vars} ({ensembl_count/adata.n_vars*100:.1f}%)")
            
            # Check for proper preservation of ensemble version in original_gene_id
            if 'original_gene_id' in adata.var.columns:
                version_pattern = sum(1 for id in adata.var['original_gene_id'] if str(id).startswith('ENSG') and '.' in str(id))
                if version_pattern > 0:
                    print(f"  Original IDs with version numbers: {version_pattern}/{adata.n_vars} ({version_pattern/adata.n_vars*100:.1f}%)")
                    print(f"    Example: {str(adata.var['original_gene_id'].iloc[0])}")
    except Exception as e:
        print(f"\n{dataset.upper()} validation error: {str(e)}")
EOF

log_message "=== Sample Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"