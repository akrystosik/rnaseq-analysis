#!/bin/bash
# Run each stage of the pipeline individually with validation

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR="${BASE_DIR}/scripts"
METADATA_DIR="${BASE_DIR}/metadata/json"

# Create a timestamp for output directory
TIMESTAMP="STAGES_$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="${BASE_DIR}/standardized_data/test_${TIMESTAMP}"
PREPROCESSED_DIR="${BASE_DIR}/preprocessed_data/test_${TIMESTAMP}"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/stages_${TIMESTAMP}.log"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$PREPROCESSED_DIR"

# Log function
log() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# Run function with logging
run() {
    log "Running: $1"
    eval "$1" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

# Validation function
validate() {
    log "Validating stage $1..."
    python3 -c "$2" 2>&1 | tee -a "$LOG_FILE"
}

log "=== Starting Stage-by-Stage Testing ==="
log "Output directory: $OUTPUT_DIR"
log "Preprocessed directory: $PREPROCESSED_DIR"

# Stage 0: Entrez to Ensembl Mapping
log "=== Stage 0: Generating Entrez to Ensembl Mapping ==="
run "python ${SCRIPTS_DIR}/entrez-to-ensembl-mapping.py --output ${METADATA_DIR}/entrez_to_ensembl_mapping.csv --species human --force"

validate "0" "
import pandas as pd
mapping_file = '${METADATA_DIR}/entrez_to_ensembl_mapping.csv'
try:
    df = pd.read_csv(mapping_file)
    print(f'Entrez mapping loaded successfully: {len(df)} entries')
    
    # Check columns
    print(f'Columns: {df.columns.tolist()}')
    
    # Check for any issues
    ensembl_prefix = sum(str(x).startswith('ENSG') for x in df['ensembl_id'])
    print(f'Entries with ENSG prefix: {ensembl_prefix}/{len(df)} ({ensembl_prefix/len(df)*100:.2f}%)')
except Exception as e:
    print(f'Error validating Entrez mapping: {e}')
"

# Stage 1: Process ENCODE dataset
log "=== Stage 1: Processing ENCODE Dataset ==="
run "python ${SCRIPTS_DIR}/standardize_datasets.py --encode-dir \"${BASE_DIR}/encode/raw_data\" --metadata-dir \"${METADATA_DIR}\" --output-dir \"$OUTPUT_DIR\""

validate "1" "
import scanpy as sc
import pandas as pd
encode_file = '${OUTPUT_DIR}/encode_standardized_v1.h5ad'
try:
    adata = sc.read_h5ad(encode_file)
    print(f'ENCODE dataset loaded successfully: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check for tissue information
    if 'tissue' in adata.obs.columns:
        missing_tissue = pd.isna(adata.obs['tissue']).sum()
        print(f'Missing tissues: {missing_tissue}/{adata.n_obs} ({missing_tissue/adata.n_obs*100:.2f}%)')
        
        # Check for 'nan' string values
        nan_string = sum(adata.obs['tissue'].astype(str) == 'nan')
        print(f'String \"nan\" tissues: {nan_string}/{adata.n_obs} ({nan_string/adata.n_obs*100:.2f}%)')
        
        # Show tissue distribution
        tissue_counts = adata.obs['tissue'].value_counts(dropna=False).head(5)
        print('Top tissues:')
        for tissue, count in tissue_counts.items():
            print(f'  {tissue}: {count}')
except Exception as e:
    print(f'Error validating ENCODE dataset: {e}')
"

# Stage 2: Process MAGE dataset
log "=== Stage 2: Processing MAGE Dataset ==="
run "python ${SCRIPTS_DIR}/standardize_datasets.py --mage-dir \"${BASE_DIR}/mage\" --metadata-dir \"${METADATA_DIR}\" --output-dir \"$OUTPUT_DIR\""

validate "2" "
import scanpy as sc
import pandas as pd
mage_file = '${OUTPUT_DIR}/mage_standardized_v1.h5ad'
try:
    adata = sc.read_h5ad(mage_file)
    print(f'MAGE dataset loaded successfully: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check for tissue information
    if 'tissue' in adata.obs.columns:
        missing_tissue = pd.isna(adata.obs['tissue']).sum()
        print(f'Missing tissues: {missing_tissue}/{adata.n_obs} ({missing_tissue/adata.n_obs*100:.2f}%)')
        
        # Show tissue distribution
        tissue_counts = adata.obs['tissue'].value_counts(dropna=False).head(5)
        print('Top tissues:')
        for tissue, count in tissue_counts.items():
            print(f'  {tissue}: {count}')
except Exception as e:
    print(f'Error validating MAGE dataset: {e}')
"

# Stage 3: Process ADNI dataset
log "=== Stage 3: Processing ADNI Dataset ==="
run "python ${SCRIPTS_DIR}/standardize_datasets.py --adni-dir \"${BASE_DIR}/adni_microarray\" --metadata-dir \"${METADATA_DIR}\" --output-dir \"$OUTPUT_DIR\""

validate "3" "
import scanpy as sc
import pandas as pd
adni_file = '${OUTPUT_DIR}/adni_standardized_v1.h5ad'
try:
    adata = sc.read_h5ad(adni_file)
    print(f'ADNI dataset loaded successfully: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check for tissue information
    if 'tissue' in adata.obs.columns:
        missing_tissue = pd.isna(adata.obs['tissue']).sum()
        print(f'Missing tissues: {missing_tissue}/{adata.n_obs} ({missing_tissue/adata.n_obs*100:.2f}%)')
        
        # Show tissue distribution
        tissue_counts = adata.obs['tissue'].value_counts(dropna=False).head(5)
        print('Top tissues:')
        for tissue, count in tissue_counts.items():
            print(f'  {tissue}: {count}')
except Exception as e:
    print(f'Error validating ADNI dataset: {e}')
"

# Stage 4: Generate Gene ID Reference Mapping
log "=== Stage 4: Generating Gene ID Reference Mapping ==="
run "python ${SCRIPTS_DIR}/gene_id_mapping_reference.py --encode-dir \"${BASE_DIR}/encode/raw_data\" --entex-dir \"${BASE_DIR}/encode/entex\" --entrez-mapping \"${METADATA_DIR}/entrez_to_ensembl_mapping.csv\" --output \"${METADATA_DIR}/gene_id_reference_mapping.csv\" --force"

validate "4" "
import pandas as pd
mapping_file = '${METADATA_DIR}/gene_id_reference_mapping.csv'
try:
    df = pd.read_csv(mapping_file)
    print(f'Gene ID reference mapping loaded successfully: {len(df)} entries')
    
    # Check columns
    print(f'Columns: {df.columns.tolist()}')
    
    # Check for placeholders
    placeholder_count = sum(str(x).startswith('PLACEHOLDER_') for x in df['gene_id'])
    print(f'Entries with PLACEHOLDER prefix: {placeholder_count}/{len(df)} ({placeholder_count/len(df)*100:.2f}%)')
except Exception as e:
    print(f'Error validating gene ID reference mapping: {e}')
"

# Stage 5: Generate ENCODE ID Mapping
log "=== Stage 5: Generating ENCODE ID Mapping ==="
run "python ${SCRIPTS_DIR}/generate_encode_mapping.py --encode-dir \"${BASE_DIR}/encode/raw_data\" --output-dir \"${METADATA_DIR}/gene_mapping\" --force"

validate "5" "
import pandas as pd
import os
mapping_file = '${METADATA_DIR}/gene_mapping/encode_id_to_ensembl_mapping.csv'
try:
    if os.path.exists(mapping_file):
        df = pd.read_csv(mapping_file)
        print(f'ENCODE ID mapping loaded successfully: {len(df)} entries')
        
        # Check columns
        print(f'Columns: {df.columns.tolist()}')
        
        # Check mapping stats
        ensembl_prefix = sum(str(x).startswith('ENSG') for x in df['ensembl_id'])
        print(f'Entries with ENSG prefix: {ensembl_prefix}/{len(df)} ({ensembl_prefix/len(df)*100:.2f}%)')
    else:
        print(f'ENCODE ID mapping file not found: {mapping_file}')
except Exception as e:
    print(f'Error validating ENCODE ID mapping: {e}')
"

# Stage 6: Preprocess Datasets
log "=== Stage 6: Preprocessing Datasets ==="
run "python ${SCRIPTS_DIR}/preprocess_dataset_gene_ids.py --data-dir \"$OUTPUT_DIR\" --reference-mapping \"${METADATA_DIR}/gene_id_reference_mapping.csv\" --output-dir \"$PREPROCESSED_DIR\" --datasets \"encode,mage,adni\" --force"

validate "6" "
import scanpy as sc
import pandas as pd
import os
import glob

# Check all preprocessed files
preprocessed_files = glob.glob('${PREPROCESSED_DIR}/*_standardized_preprocessed.h5ad')
for file_path in preprocessed_files:
    dataset_name = os.path.basename(file_path).split('_')[0]
    try:
        adata = sc.read_h5ad(file_path)
        print(f'\\n{dataset_name.upper()} preprocessed dataset: {adata.shape[0]} samples, {adata.shape[1]} genes')
        
        # Check for placeholder IDs
        if 'ensembl_id' in adata.var.columns:
            placeholder_count = sum(isinstance(x, str) and x.startswith('PLACEHOLDER_') for x in adata.var['ensembl_id'])
            print(f'Placeholder gene IDs: {placeholder_count}/{adata.n_vars} ({placeholder_count/adata.n_vars*100:.2f}%)')
            
            # Show mapping sources if available
            if 'mapping_source' in adata.var.columns:
                source_counts = adata.var['mapping_source'].value_counts()
                print('Mapping sources:')
                for source, count in source_counts.items():
                    print(f'  {source}: {count} ({count/adata.n_vars*100:.2f}%)')
    except Exception as e:
        print(f'Error validating {dataset_name} preprocessed dataset: {e}')
"

# Stage 7: Create Combined Dataset
log "=== Stage 7: Creating Combined Dataset ==="
run "python ${SCRIPTS_DIR}/create_combined_dataset_all_genes_sparse.py --input-dir \"$PREPROCESSED_DIR\" --reference-mapping \"${METADATA_DIR}/gene_id_reference_mapping.csv\" --output-file \"${OUTPUT_DIR}/combined_all_genes_sparse_standardized.h5ad\" --include-datasets \"encode,mage,adni\" --force"

validate "7" "
import scanpy as sc
import pandas as pd
combined_file = '${OUTPUT_DIR}/combined_all_genes_sparse_standardized.h5ad'
try:
    adata = sc.read_h5ad(combined_file)
    print(f'Combined dataset loaded successfully: {adata.shape[0]} samples, {adata.shape[1]} genes')
    
    # Check dataset distribution
    if 'dataset' in adata.obs.columns:
        print('Dataset distribution:')
        for dataset, count in adata.obs['dataset'].value_counts().items():
            print(f'  {dataset}: {count} ({count/adata.n_obs*100:.2f}%)')
        
        # Check for 'combined' label
        combined_count = sum(adata.obs['dataset'] == 'combined')
        if combined_count > 0:
            print(f'WARNING: Found {combined_count} samples labeled as \"combined\"')
        else:
            print('No samples labeled as \"combined\" - issue fixed!')
except Exception as e:
    print(f'Error validating combined dataset: {e}')
"

log "=== Stage-by-Stage Testing Complete ==="
log "Results saved to $OUTPUT_DIR"
log "Log file: $LOG_FILE"
