#!/bin/bash
# Comprehensive test for the RNA-seq standardization pipeline across all datasets

# Create timestamp for output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/test_full_pipeline_${TIMESTAMP}"
LOG_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs/test_full_pipeline_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/preprocessed_data"
mkdir -p "$(dirname $LOG_FILE)"

# Log function
log_message() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

log_message "=== Starting Full Pipeline Test Across All Datasets ==="
log_message "Output directory: $OUTPUT_DIR"

# Step 1: Generate Entrez to Ensembl Mapping
log_message "Step 1: Generating Entrez to Ensembl Mapping"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/entrez-to-ensembl-mapping.py \
    --output "$OUTPUT_DIR/entrez_to_ensembl_mapping.csv" \
    --species human 2>&1 | tee -a "$LOG_FILE"

# Step 2: Run standardize_datasets.py for all datasets
log_message "Step 2: Running standardize_datasets.py for all datasets"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py \
    --encode-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines" \
    --gtex-file "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz" \
    --mage-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage" \
    --adni-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni_microarray" \
    --metadata-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json" \
    --output-dir "$OUTPUT_DIR" 2>&1 | tee -a "$LOG_FILE"

# Step 3: Use save wrapper to properly save AnnData objects
log_message "Step 3: Using save wrapper for all datasets"
for dataset in encode gtex mage adni; do
    dataset_file="$OUTPUT_DIR/${dataset}_standardized_v1.h5ad"
    if [ -f "$dataset_file" ]; then
        log_message "Using save wrapper for $dataset dataset"
        python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/anndata_save_wrapper.py \
            "$dataset_file" "${dataset_file}.fixed" 2>&1 | tee -a "$LOG_FILE"
        mv "${dataset_file}.fixed" "$dataset_file"
    else
        log_message "WARNING: $dataset file not found"
    fi
done

# Step 4: Generate gene ID reference mapping
log_message "Step 4: Generating gene ID reference mapping"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/gene_id_mapping_reference.py \
    --encode-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data" \
    --entex-dir "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/entex" \
    --entrez-mapping "$OUTPUT_DIR/entrez_to_ensembl_mapping.csv" \
    --output "$OUTPUT_DIR/gene_id_reference_mapping.csv" \
    --force 2>&1 | tee -a "$LOG_FILE"

# Step 5: Run preprocess_dataset_gene_ids.py
log_message "Step 5: Running preprocess_dataset_gene_ids.py for all datasets"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/preprocess_dataset_gene_ids.py \
    --data-dir "$OUTPUT_DIR" \
    --reference-mapping "$OUTPUT_DIR/gene_id_reference_mapping.csv" \
    --output-dir "$OUTPUT_DIR/preprocessed_data" \
    --datasets "encode,gtex,mage,adni" 2>&1 | tee -a "$LOG_FILE"

# Step 6: Fix placeholder IDs
log_message "Step 6: Fixing placeholder IDs in all datasets"
for dataset in encode gtex mage adni; do
    preprocessed_file="$OUTPUT_DIR/preprocessed_data/${dataset}_standardized_preprocessed.h5ad"
    if [ -f "$preprocessed_file" ]; then
        log_message "Fixing placeholder IDs in $dataset dataset"
        python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/fix_placeholder_ids.py \
            "$preprocessed_file" "${preprocessed_file}.fixed" 2>&1 | tee -a "$LOG_FILE"
        mv "${preprocessed_file}.fixed" "$preprocessed_file"
    else
        log_message "WARNING: Preprocessed $dataset file not found"
    fi
done

# Step 7: Create a small combined dataset with common genes
log_message "Step 7: Creating small combined dataset with common genes"
python /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/create_combined_dataset.py \
    --output-file "$OUTPUT_DIR/combined_common_genes_standardized.h5ad" 2>&1 | tee -a "$LOG_FILE"

# Step 8: Analyze results for each dataset
log_message "Step 8: Analyzing results for all datasets"

# Function to analyze a dataset
analyze_dataset() {
    local dataset=$1
    local stage=$2  # "original" or "preprocessed"
    
    if [ "$stage" == "original" ]; then
        file_path="$OUTPUT_DIR/${dataset}_standardized_v1.h5ad"
    else
        file_path="$OUTPUT_DIR/preprocessed_data/${dataset}_standardized_preprocessed.h5ad"
    fi
    
    log_message "Analyzing $stage $dataset dataset: $file_path"
    
    if [ -f "$file_path" ]; then
        python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np

print(f'\\n===== {\"$dataset\".upper()} {\"$stage\".upper()} DATASET ANALYSIS =====')
try:
    adata = sc.read_h5ad('$file_path')
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
    
    # Check gene ID format in preprocessed datasets
    if '${stage}' == 'preprocessed':
        if 'gene_id' in adata.var.columns and 'ensembl_id' in adata.var.columns:
            # Count placeholder IDs
            placeholder_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('PLACEHOLDER_'))
            entrez_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('ENTREZ:'))
            ensembl_count = sum(1 for id in adata.var['ensembl_id'] if str(id).startswith('ENSG'))
            
            print('\\nGene ID analysis:')
            print(f'  Placeholder IDs: {placeholder_count}/{adata.n_vars} ({placeholder_count/adata.n_vars*100:.2f}%)')
            print(f'  Entrez IDs: {entrez_count}/{adata.n_vars} ({entrez_count/adata.n_vars*100:.2f}%)')
            print(f'  Ensembl IDs: {ensembl_count}/{adata.n_vars} ({ensembl_count/adata.n_vars*100:.2f}%)')
            
            # Check if gene_id matches ensembl_id
            matching_ids = sum(adata.var['gene_id'] == adata.var['ensembl_id'])
            print(f'  gene_id matches ensembl_id: {matching_ids}/{adata.n_vars} ({matching_ids/adata.n_vars*100:.2f}%)')
            
            # Check version numbers
            if 'original_gene_id' in adata.var.columns:
                version_pattern = r'ENSG\\d+\\.\\d+'
                original_versions = sum(adata.var['original_gene_id'].astype(str).str.match(version_pattern))
                print(f'  Original IDs with version numbers: {original_versions}/{adata.n_vars} ({original_versions/adata.n_vars*100:.2f}%)')
        
        # Check mapping source distribution
        if 'mapping_source' in adata.var.columns:
            mapping_counts = adata.var['mapping_source'].value_counts()
            print('\\nMapping source distribution:')
            for source, count in mapping_counts.items():
                print(f'  {source}: {count} ({count/adata.n_vars*100:.2f}%)')
    
    # Check for common columns that should be present
    expected_obs_cols = ['dataset', 'sample_id', 'data_type', 'expression_unit']
    expected_var_cols = []
    if '${stage}' == 'preprocessed':
        expected_var_cols = ['gene_id', 'ensembl_id', 'gene_name', 'gene_type', 'chromosome', 'mapping_source']
    
    print('\\nMetadata completeness:')
    for col in expected_obs_cols:
        if col in adata.obs.columns:
            print(f'  obs[{col}]: Present')
        else:
            print(f'  obs[{col}]: MISSING')
    
    for col in expected_var_cols:
        if col in adata.var.columns:
            print(f'  var[{col}]: Present')
        else:
            print(f'  var[{col}]: MISSING')
    
    # Sample metadata
    if adata.n_obs > 0:
        print('\\nSample metadata preview:')
        preview_cols = min(5, len(adata.obs.columns))
        print(adata.obs.iloc[0:1, 0:preview_cols])

except Exception as e:
    print(f'Error analyzing dataset: {e}')
    import traceback
    print(traceback.format_exc())
" 2>&1 | tee -a "$LOG_FILE"
    else
        log_message "ERROR: $stage $dataset file not found: $file_path"
    fi
}

# Analyze original and preprocessed versions of each dataset
for dataset in encode gtex mage adni; do
    analyze_dataset "$dataset" "original"
    analyze_dataset "$dataset" "preprocessed"
done

# Analyze combined dataset
if [ -f "$OUTPUT_DIR/combined_common_genes_standardized.h5ad" ]; then
    log_message "Analyzing combined dataset"
    python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np

print('\\n===== COMBINED DATASET ANALYSIS =====')
try:
    adata = sc.read_h5ad('$OUTPUT_DIR/combined_common_genes_standardized.h5ad')
    print(f'Combined dataset loaded: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Dataset distribution
    if 'dataset' in adata.obs.columns:
        dataset_counts = adata.obs['dataset'].value_counts()
        print('\\nSamples per dataset:')
        for dataset, count in dataset_counts.items():
            print(f'  {dataset}: {count} ({count/adata.n_obs*100:.2f}%)')
    
    # Check gene ID format
    if adata.n_vars > 0:
        ensembl_count = sum(1 for id in adata.var_names if str(id).startswith('ENSG'))
        other_count = adata.n_vars - ensembl_count
        
        print('\\nGene ID format:')
        print(f'  Ensembl IDs: {ensembl_count}/{adata.n_vars} ({ensembl_count/adata.n_vars*100:.2f}%)')
        print(f'  Other IDs: {other_count}/{adata.n_vars} ({other_count/adata.n_vars*100:.2f}%)')
    
    # Print dataset information
    if 'dataset_info' in adata.uns:
        print('\\nDataset Info:')
        for key, value in adata.uns['dataset_info'].items():
            if key != 'sample_counts':  # Skip detailed sample counts
                print(f'  {key}: {value}')
    
except Exception as e:
    print(f'Error analyzing combined dataset: {e}')
    import traceback
    print(traceback.format_exc())
" 2>&1 | tee -a "$LOG_FILE"
else
    log_message "ERROR: Combined dataset file not found"
fi

log_message "=== Test Complete ==="
log_message "Results saved to $OUTPUT_DIR"
log_message "Log file: $LOG_FILE"