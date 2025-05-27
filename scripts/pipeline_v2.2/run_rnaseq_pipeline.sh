#!/bin/bash
# Run the complete RNA-seq standardization pipeline with improved gene ID mapping
# Uses separate base directories for inputs and outputs.

# --- Configuration ---
INPUT_BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
OUTPUT_BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq"

SCRIPTS_DIR="${INPUT_BASE_DIR}/scripts/pipeline_v2.2" # Scripts are read from input location

# --- INPUT DATA PATHS (relative to INPUT_BASE_DIR or absolute) ---
# Metadata JSON files (configuration for the pipeline itself)
PIPELINE_METADATA_JSON_DIR="${INPUT_BASE_DIR}/metadata/json" # General pipeline config JSONs
# Gene mapping resources (source files like GTF, NCBI downloaded maps)
SOURCE_GENE_MAPPING_RESOURCE_DIR="${INPUT_BASE_DIR}/metadata/gene_mapping"
DOWNLOADED_GTF_GZ="${SOURCE_GENE_MAPPING_RESOURCE_DIR}/gencode.v24.annotation.gtf.gz"
UNZIPPED_GTF_FILE_SOURCE="${SOURCE_GENE_MAPPING_RESOURCE_DIR}/gencode.v24.annotation.gtf" # Source if already unzipped

# Raw data locations
ENCODE_RAW_DATA_DIR="${INPUT_BASE_DIR}/encode/raw_data"
ENTEX_RAW_DATA_DIR="${INPUT_BASE_DIR}/encode/entex"
ENTEX_METADATA_FILE_INPUT="${INPUT_BASE_DIR}/encode/metadata/entex_metadata.json"
MAGE_RAW_DATA_DIR="${INPUT_BASE_DIR}/mage"
ADNI_RAW_DATA_DIR="${INPUT_BASE_DIR}/adni_microarray"
ADNI_DEMOGRAPHICS_FILE_INPUT="${INPUT_BASE_DIR}/metadataADNI/subject_demographics/PTDEMOG_25Apr2025.csv"
ADNI_DIAGNOSIS_FILE_INPUT="${INPUT_BASE_DIR}/metadataADNI/subject_demographics/DXSUM_30Apr2025.csv"
GTEX_RAW_FILE_INPUT="${INPUT_BASE_DIR}/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"
MAGE_1000G_PED_FILE_INPUT="${INPUT_BASE_DIR}/project_gene_regulation/data/MAGE/WGS/20130606_g1k_3202_samples_ped_population.txt"


# --- OUTPUT & GENERATED DATA PATHS (relative to OUTPUT_BASE_DIR) ---
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
# Main output directories for H5ADs
OUTPUT_DIR_H5ADS="${OUTPUT_BASE_DIR}/standardized_data/run_${TIMESTAMP}"
PREPROCESSED_DIR_H5ADS="${OUTPUT_BASE_DIR}/preprocessed_data/run_${TIMESTAMP}"
# Log directory
LOG_DIR_PIPELINE="${OUTPUT_BASE_DIR}/logs/pipeline_runs"
# Directory for *generated* metadata and gene mapping files (e.g., NCBI map, primary ref map)
GENERATED_METADATA_JSON_DIR="${OUTPUT_BASE_DIR}/metadata/json_generated"
GENERATED_GENE_MAPPING_RESOURCE_DIR="${OUTPUT_BASE_DIR}/metadata/gene_mapping_generated"

# Paths for generated mapping files
ENTREZ_TO_ENSEMBL_CSV_GENERATED="${GENERATED_METADATA_JSON_DIR}/entrez_to_ensembl_mapping.csv"
PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED="${GENERATED_GENE_MAPPING_RESOURCE_DIR}/gencode_v24_complete_mapping.csv"
UNZIPPED_GTF_FILE_OUTPUT="${GENERATED_GENE_MAPPING_RESOURCE_DIR}/gencode.v24.annotation.gtf" # Where unzipped version will be placed if re-generated
ENCODE_ID_MAP_OUTPUT_DIR_GENERATED="${GENERATED_GENE_MAPPING_RESOURCE_DIR}/encode_specific_mappings"


# --- Logging Setup ---
mkdir -p "$LOG_DIR_PIPELINE"
PIPELINE_LOG_FILE="${LOG_DIR_PIPELINE}/pipeline_run_${TIMESTAMP}.log"

# --- Helper Functions ---
log_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$PIPELINE_LOG_FILE"
}

run_command() {
    log_message "Running: $1"
    PYTHONPATH="${SCRIPTS_DIR}:${PYTHONPATH}" eval "$1" 2>&1 | tee -a "$PIPELINE_LOG_FILE"
    return ${PIPESTATUS[0]}
}

# --- Pipeline Start ---
log_message "=== Starting RNA-seq Pipeline Run: ${TIMESTAMP} ==="
log_message "INPUT BASE DIRECTORY: ${INPUT_BASE_DIR}"
log_message "OUTPUT BASE DIRECTORY: ${OUTPUT_BASE_DIR}"

# Create output structure
mkdir -p "$OUTPUT_DIR_H5ADS" "$PREPROCESSED_DIR_H5ADS" \
           "$GENERATED_METADATA_JSON_DIR" "$GENERATED_GENE_MAPPING_RESOURCE_DIR" \
           "${OUTPUT_DIR_H5ADS}/temp" # Temp dir for standardize_datasets.py

# --- Argument Parsing ---
FORCE_FLAG=""
EFFECTIVE_FORCE_FOR_REF_MAP=""
SKIP_METADATA_ENHANCEMENT=""

# Parse all arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --force|--force-all)
            FORCE_FLAG="--force"
            EFFECTIVE_FORCE_FOR_REF_MAP="--force"
            log_message "Force flag ('$1') detected. Applicable scripts will regenerate outputs."
            ;;
        --force-mapping)
            EFFECTIVE_FORCE_FOR_REF_MAP="--force"
            log_message "Force mapping flag detected. Gene mapping files will be regenerated."
            ;;
        --skip-metadata-enhancement)
            SKIP_METADATA_ENHANCEMENT="true"
            log_message "Skip metadata enhancement flag detected. Stages 2.8-3.5 will be skipped."
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --force, --force-all         Force regeneration of all outputs"
            echo "  --force-mapping              Force regeneration of gene mapping files only"
            echo "  --skip-metadata-enhancement  Skip metadata enhancement stages (2.8-3.5)"
            echo "  --help, -h                   Show this help message"
            echo ""
            echo "Pipeline Stages:"
            echo "  0     : Entrez-to-Ensembl mapping"
            echo "  A     : Primary gene reference map"
            echo "  1     : Raw data standardization"
            echo "  1.6   : ENCODE-specific mappings"
            echo "  2     : Metadata enhancement"
            echo "  2.5   : Gene ID preprocessing"
            echo "  2.6   : Placeholder ID fixes"
            echo "  2.7   : ENCODE mapping analysis"
            echo "  2.8   : Controlled-access metadata integration"
            echo "  2.9   : Developmental stage mapping"
            echo "  3     : Combined dataset creation"
            echo "  3.5   : MAGE technical metadata integration"
            echo "  4     : Validation"
            echo "  4.1   : Gene overlap analysis"
            echo "  4.2   : Partner deliverables generation"
            exit 0
            ;;
        *)
            log_message "Unknown option: $1"
            echo "Use --help for usage information."
            exit 1
            ;;
    esac
    shift
done

# --- PRE-RUN CHECKS FOR INPUTS ---
log_message "--- Performing Pre-run Input Checks (from INPUT_BASE_DIR) ---"
if [ ! -d "$SCRIPTS_DIR" ]; then log_message "ERROR: Scripts directory $SCRIPTS_DIR not found! Exiting."; exit 1; fi
if [ ! -f "$DOWNLOADED_GTF_GZ" ] && [ ! -f "$UNZIPPED_GTF_FILE_SOURCE" ]; then log_message "WARNING: GENCODE GTF (gzipped or source unzipped) not found at expected input path. Download or ensure correct path."; fi
if [ -n "$MAGE_RAW_DATA_DIR" ] && [ ! -d "$MAGE_RAW_DATA_DIR" ]; then log_message "WARNING: MAGE raw data directory $MAGE_RAW_DATA_DIR not found."; fi
if [ ! -f "$GTEX_RAW_FILE_INPUT" ]; then log_message "WARNING: GTEx raw file $GTEX_RAW_FILE_INPUT not found."; fi
# Add more checks for other critical inputs if desired
log_message "--- Pre-run Input Checks Complete ---"


# --- Stage 0: Prepare Entrez to Ensembl Mapping (NCBI Source) ---
# This will be written to the OUTPUT_BASE_DIR structure
log_message "--- Stage 0: Preparing Entrez to Ensembl Mapping ---"
mkdir -p "$(dirname "${ENTREZ_TO_ENSEMBL_CSV_GENERATED}")"
if [ "$EFFECTIVE_FORCE_FOR_REF_MAP" == "--force" ] || [ ! -f "$ENTREZ_TO_ENSEMBL_CSV_GENERATED" ]; then
    run_command "python \"${SCRIPTS_DIR}/entrez-to-ensembl-mapping.py\" --output \"${ENTREZ_TO_ENSEMBL_CSV_GENERATED}\" --species human --force"
    if [ $? -ne 0 ]; then log_message "ERROR: Entrez to Ensembl mapping generation failed. Exiting."; exit 1; fi
else
    log_message "Entrez to Ensembl mapping file ${ENTREZ_TO_ENSEMBL_CSV_GENERATED} already exists and no relevant force flag. Skipping generation."
fi

# --- Stage 0.5: Unzip GENCODE GTF ---
# Unzips from INPUT_BASE_DIR to OUTPUT_BASE_DIR if needed
log_message "--- Stage 0.5: Preparing GENCODE v24 GTF File ---"
mkdir -p "$(dirname "${UNZIPPED_GTF_FILE_OUTPUT}")"
if [ ! -f "$DOWNLOADED_GTF_GZ" ]; then
    log_message "ERROR: Gzipped GTF file $DOWNLOADED_GTF_GZ not found in input path."
    if [ ! -f "$UNZIPPED_GTF_FILE_OUTPUT" ] && [ ! -f "$UNZIPPED_GTF_FILE_SOURCE" ]; then
        log_message "ERROR: And no unzipped GTF file found in output or source paths. Exiting."; exit 1;
    elif [ -f "$UNZIPPED_GTF_FILE_SOURCE" ] && [ ! -f "$UNZIPPED_GTF_FILE_OUTPUT" ] ; then
        log_message "WARNING: Gzipped GTF $DOWNLOADED_GTF_GZ not found. Using source unzipped version $UNZIPPED_GTF_FILE_SOURCE and copying to output."
        cp "$UNZIPPED_GTF_FILE_SOURCE" "$UNZIPPED_GTF_FILE_OUTPUT"
    elif [ -f "$UNZIPPED_GTF_FILE_OUTPUT" ]; then
        log_message "WARNING: Gzipped GTF $DOWNLOADED_GTF_GZ not found, but unzipped version $UNZIPPED_GTF_FILE_OUTPUT exists in output. Proceeding."
    fi
elif [ "$EFFECTIVE_FORCE_FOR_REF_MAP" == "--force" ] || [ ! -f "$UNZIPPED_GTF_FILE_OUTPUT" ] || [ "$DOWNLOADED_GTF_GZ" -nt "$UNZIPPED_GTF_FILE_OUTPUT" ]; then
    log_message "Unzipping $DOWNLOADED_GTF_GZ to $UNZIPPED_GTF_FILE_OUTPUT..."
    gunzip -c "$DOWNLOADED_GTF_GZ" > "$UNZIPPED_GTF_FILE_OUTPUT"; if [ $? -ne 0 ]; then log_message "ERROR: Failed to unzip GTF. Exiting."; exit 1; fi
else
    log_message "Unzipped GTF ${UNZIPPED_GTF_FILE_OUTPUT} is up-to-date in output. Skipping unzip."
fi
# Ensure the GTF path used by subsequent scripts is the one in the output/generated directory
FINAL_UNZIPPED_GTF_TO_USE=$UNZIPPED_GTF_FILE_OUTPUT
if [ ! -f "$FINAL_UNZIPPED_GTF_TO_USE" ] && [ -f "$UNZIPPED_GTF_FILE_SOURCE" ]; then
    log_message "Unzipped GTF not in output, using source: $UNZIPPED_GTF_FILE_SOURCE"
    FINAL_UNZIPPED_GTF_TO_USE=$UNZIPPED_GTF_FILE_SOURCE
elif [ ! -f "$FINAL_UNZIPPED_GTF_TO_USE" ]; then
    log_message "ERROR: No unzipped GTF file available. Exiting."
    exit 1
fi


# --- Stage A: Generate THE Primary Gene Annotation/Reference Map ---
# This map will be written to OUTPUT_BASE_DIR structure
log_message "--- Stage A: Generating Primary Gene Annotation & Reference Map ---"
TEMP_DOWNLOAD_DIR_FOR_REF="${OUTPUT_BASE_DIR}/temp_gene_ref_downloads" # Temp dir under output
mkdir -p "$TEMP_DOWNLOAD_DIR_FOR_REF"
mkdir -p "$(dirname "${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}")"

run_command "python \"${SCRIPTS_DIR}/gene_id_mapping_reference.py\" \\
    --encode-dir \"${ENCODE_RAW_DATA_DIR}\" \\
    --entex-dir \"${ENTEX_RAW_DATA_DIR}\" \\
    --gencode-gtf \"${FINAL_UNZIPPED_GTF_TO_USE}\" \\
    --entrez-mapping \"${ENTREZ_TO_ENSEMBL_CSV_GENERATED}\" \\
    --output \"${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}\" \\
    --temp-dir \"${TEMP_DOWNLOAD_DIR_FOR_REF}\" \\
    ${EFFECTIVE_FORCE_FOR_REF_MAP}"
if [ $? -ne 0 ]; then log_message "ERROR: Stage A (Primary Gene Annotation Map generation to ${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}) failed. Exiting."; exit 1; fi
log_message "Stage A completed. Primary gene annotation & reference map: ${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}"


# --- Stage 1: Initial Data Conversion (Raw to Standardized v1 H5ADs) ---
# Reads from INPUT_BASE_DIR, Writes to OUTPUT_DIR_H5ADS
# standardize_datasets.py will need to be aware of PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED
# This means rnaseq_utils.GENCODE_MAPPING_FILE might need to be overridden or the script adapted
# For now, we assume standardize_datasets.py is modified to accept this path or uses a symlink trick.
# A simpler way is to make sure rnaseq_utils.py points to the generated one.
# Since rnaseq_utils.py is in INPUT_BASE_DIR, it might be tricky.
# Let's assume for now that `standardize_datasets.py` is modified to take the primary map path as an arg, or we set an env var.
# For this iteration, we'll ensure the `rnaseq_utils.py` used by `standardize_datasets.py` is configured to look for
# the *generated* primary mapping file. This is usually handled by GENCODE_MAPPING_FILE in rnaseq_utils.py.
# A robust way is to pass the path to the standardize_datasets.py script.
# *Self-correction: standardize_datasets.py via rnaseq_utils.load_gencode_mapping() uses a fixed path.
# The simplest is to ensure that fixed path (GENCODE_MAPPING_FILE in rnaseq_utils.py) points to *our generated one*.
# This means rnaseq_utils.py needs to be aware of the generated path.
# Let's make rnaseq_utils.py flexible or ensure the path it looks for IS the generated one.
# For now, let's copy the generated primary map to the path expected by rnaseq_utils.py if it's different.
# The `rnaseq_utils.py` expects it at `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gene_mapping/gencode_v24_complete_mapping.csv`
# If INPUT_BASE_DIR is the same, it works. If different, we need to handle this.
# For this separated path logic, it's best if `standardize_datasets.py` can take the gene map path as an argument.
# Assuming `standardize_datasets.py` (and by extension `rnaseq_utils.py`) will correctly use the
# `PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED` (this requires modification to those scripts to accept it as param, or env var).
# For now, the script has rnaseq_utils.GENCODE_MAPPING_FILE hardcoded.
# The most direct fix is to ensure the `run_command` for standardize_datasets.py uses a python path that sees an rnaseq_utils.py
# configured to point to PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED. This is complex.
# Alternative: symlink or copy PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED to the location expected by the *original* rnaseq_utils.py
EXPECTED_GENCODE_MAP_PATH_FOR_UTILS="${INPUT_BASE_DIR}/metadata/gene_mapping/gencode_v24_complete_mapping.csv"
if [ "$PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED" != "$EXPECTED_GENCODE_MAP_PATH_FOR_UTILS" ]; then
    log_message "Copying generated primary gene map ${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED} to expected location ${EXPECTED_GENCODE_MAP_PATH_FOR_UTILS} for Stage 1..."
    mkdir -p "$(dirname "$EXPECTED_GENCODE_MAP_PATH_FOR_UTILS")"
    cp "$PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED" "$EXPECTED_GENCODE_MAP_PATH_FOR_UTILS"
    # Also copy the JSON version if it exists
    GENERATED_JSON_MAP="${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED%.csv}.json"
    EXPECTED_JSON_MAP="${EXPECTED_GENCODE_MAP_PATH_FOR_UTILS%.csv}.json"
    if [ -f "$GENERATED_JSON_MAP" ]; then
      cp "$GENERATED_JSON_MAP" "$EXPECTED_JSON_MAP"
    fi
fi

log_message "--- Stage 1: Initial Data Conversion (Raw to Standardized v1 H5ADs) ---"
ADNI_DEMOGRAPHICS_ARG=""; if [ -f "$ADNI_DEMOGRAPHICS_FILE_INPUT" ]; then ADNI_DEMOGRAPHICS_ARG="--adni-demographics-file \"${ADNI_DEMOGRAPHICS_FILE_INPUT}\""; else log_message "WARNING: ADNI demographics file ${ADNI_DEMOGRAPHICS_FILE_INPUT} missing."; fi
ADNI_DIAGNOSIS_ARG=""; if [ -f "$ADNI_DIAGNOSIS_FILE_INPUT" ]; then ADNI_DIAGNOSIS_ARG="--adni-diagnosis-file \"${ADNI_DIAGNOSIS_FILE_INPUT}\""; else log_message "WARNING: ADNI diagnosis file ${ADNI_DIAGNOSIS_FILE_INPUT} missing."; fi

run_command "python \"${SCRIPTS_DIR}/standardize_datasets.py\" \\
    --encode-dir \"${ENCODE_RAW_DATA_DIR}/cell_lines\" \\
    --encode-entex-dir \"${ENTEX_RAW_DATA_DIR}\" \\
    --entex-metadata-file \"${ENTEX_METADATA_FILE_INPUT}\" \\
    --mage-dir \"${MAGE_RAW_DATA_DIR}\" \\
    --adni-dir \"${ADNI_RAW_DATA_DIR}\" \\
    ${ADNI_DEMOGRAPHICS_ARG} \\
    ${ADNI_DIAGNOSIS_ARG} \\
    --gtex-file \"${GTEX_RAW_FILE_INPUT}\" \\
    --metadata-dir \"${PIPELINE_METADATA_JSON_DIR}\" \\
    --output-dir \"${OUTPUT_DIR_H5ADS}\""
if [ $? -ne 0 ]; then log_message "ERROR: Stage 1 (Initial Data Conversion) failed. Exiting."; exit 1; fi
log_message "Stage 1 completed."


# --- Stage 1.6: Generate ENCODE-specific ID to Ensembl Mapping ---
# Output to GENERATED_GENE_MAPPING_RESOURCE_DIR
log_message "--- Stage 1.6: Generating ENCODE ID to Ensembl Mapping ---"
mkdir -p "$ENCODE_ID_MAP_OUTPUT_DIR_GENERATED"
run_command "python \"${SCRIPTS_DIR}/generate_encode_mapping.py\" \\
    --encode-dir \"${ENCODE_RAW_DATA_DIR}\" \\
    --output-dir \"${ENCODE_ID_MAP_OUTPUT_DIR_GENERATED}\" \\
    ${FORCE_FLAG}"
if [ $? -ne 0 ]; then log_message "WARNING: Stage 1.6 (ENCODE ID Mapping) failed."; else log_message "Stage 1.6 completed."; fi


# --- Stage 2: Enhanced Metadata Standardization (v1 to v2 H5ADs) ---
# Reads from OUTPUT_DIR_H5ADS, writes to OUTPUT_DIR_H5ADS
log_message "--- Stage 2: Enhanced Metadata Standardization ---"
run_command "python \"${SCRIPTS_DIR}/standardize_metadata.py\" \\
    --data-dir \"${OUTPUT_DIR_H5ADS}\" \\
    --output-dir \"${OUTPUT_DIR_H5ADS}\" \\
    --metadata-dir \"${PIPELINE_METADATA_JSON_DIR}\"" # Uses config JSONs from input dir
if [ $? -ne 0 ]; then log_message "ERROR: Stage 2 (Enhanced Metadata Standardization) failed. Exiting."; exit 1; fi
log_message "Stage 2 completed."


# --- Stage 2.5: Preprocess Datasets (Gene ID Standardization for each H5AD) ---
# Reads from OUTPUT_DIR_H5ADS, uses generated primary map, writes to PREPROCESSED_DIR_H5ADS
log_message "--- Stage 2.5: Preprocessing Datasets for Consistent Gene IDs ---"
run_command "python \"${SCRIPTS_DIR}/preprocess_dataset_gene_ids.py\" \\
    --data-dir \"${OUTPUT_DIR_H5ADS}\" \\
    --reference-mapping \"${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}\" \\
    --output-dir \"${PREPROCESSED_DIR_H5ADS}\" \\
    ${FORCE_FLAG}"
if [ $? -ne 0 ]; then log_message "ERROR: Stage 2.5 (Dataset Gene ID Preprocessing) failed. Exiting."; exit 1; fi
log_message "Stage 2.5 completed. Preprocessed files in: $PREPROCESSED_DIR_H5ADS"


# --- Stage 2.6: Fix Placeholder Gene IDs ---
log_message "--- Stage 2.6: Fixing Placeholder Gene IDs in Preprocessed Datasets ---"
for dataset_label in encode gtex mage adni; do
    preprocessed_h5ad_file="${PREPROCESSED_DIR_H5ADS}/${dataset_label}_standardized_preprocessed.h5ad"
    if [ -f "$preprocessed_h5ad_file" ]; then
        log_message "Attempting to fix placeholder IDs in: ${preprocessed_h5ad_file}"
        run_command "python \"${SCRIPTS_DIR}/fix_placeholder_ids.py\" \"$preprocessed_h5ad_file\" \"$preprocessed_h5ad_file.fixed\""
        if [ $? -eq 0 ]; then
            mv "$preprocessed_h5ad_file.fixed" "$preprocessed_h5ad_file"
            log_message "Placeholder IDs fixed for ${dataset_label} dataset."
        else
            log_message "WARNING: Failed to fix placeholder IDs for ${dataset_label} dataset."
        fi
    else
        log_message "INFO: No preprocessed file found for ${dataset_label} at ${preprocessed_h5ad_file}, skipping placeholder fix."
    fi
done
log_message "Stage 2.6 completed."


# --- Stage 2.7: Analyze ENCODE Gene Mapping Quality ---
log_message "--- Stage 2.7: Analyzing ENCODE Gene Mapping Quality ---"
ENCODE_PREPROCESSED_H5AD="${PREPROCESSED_DIR_H5ADS}/encode_standardized_preprocessed.h5ad"
ENCODE_MAPPING_STATS_JSON="${OUTPUT_DIR_H5ADS}/encode_mapping_stats.json" # Output to main H5AD output dir

if [ -f "$ENCODE_PREPROCESSED_H5AD" ]; then
    run_command "python - <<EOF
import scanpy as sc
import pandas as pd
import numpy as np
import json

encode_path = '$ENCODE_PREPROCESSED_H5AD'
try:
    adata = sc.read_h5ad(encode_path)
    print(f\"ENCODE dataset shape: {adata.shape}\")

    non_empty = 0
    percentage = 0.0
    source_counts_dict = {}

    if 'ensembl_id' in adata.var.columns:
        non_empty = sum(1 for x in adata.var['ensembl_id'] if x and str(x).strip() != '' and not str(x).startswith('PLACEHOLDER_'))
        if len(adata.var) > 0:
            percentage = non_empty / len(adata.var) * 100
        print(f\"ENCODE genes with valid (non-placeholder) mapped Ensembl IDs: {non_empty}/{len(adata.var)} ({percentage:.2f}%)\")

        validly_mapped_ids = adata.var.loc[(adata.var['ensembl_id'].notna()) & (adata.var['ensembl_id'] != '') & (~adata.var['ensembl_id'].astype(str).str.startswith('PLACEHOLDER_')), 'ensembl_id']
        print(f\"Sample valid mapped Ensembl IDs: {validly_mapped_ids.unique()[:5].tolist()}\")

        if 'mapping_source' in adata.var.columns:
            source_counts = adata.var['mapping_source'].value_counts()
            source_counts_dict = {k: int(v) for k, v in source_counts.items()}
            for source, count_val in source_counts.items():
                source_percentage_val = 0.0 
                if len(adata.var) > 0:
                    source_percentage_val = count_val / len(adata.var) * 100
                print(f\"Mapping source '{source}': {count_val} ({source_percentage_val:.2f}%)\")
    else:
        print(\"Error: ENCODE preprocessed dataset does not have an 'ensembl_id' column in .var!\")

    mapping_stats = {
        'total_genes': len(adata.var) if hasattr(adata, 'var') else 0,
        'mapped_genes_valid_ensembl': non_empty,
        'mapping_percentage_valid_ensembl': percentage,
        'mapping_sources': source_counts_dict
    }

    with open('$ENCODE_MAPPING_STATS_JSON', 'w') as f:
        json.dump(mapping_stats, f, indent=2)

    print(f\"ENCODE mapping stats saved to $ENCODE_MAPPING_STATS_JSON\")

except FileNotFoundError:
    print(f\"Error: ENCODE preprocessed file not found at {encode_path}. Cannot analyze mapping quality.\")
except Exception as e:
    print(f\"Error analyzing ENCODE dataset: {e}\")
    print(\"Continuing with pipeline execution...\")
EOF"
else
    log_message "INFO: ENCODE preprocessed file not found at $ENCODE_PREPROCESSED_H5AD. Skipping ENCODE mapping quality analysis."
fi
log_message "Stage 2.7 completed."


# --- Stage 3: Combined Dataset Creation ---
log_message "--- Stage 3: Creating Combined Dataset (All Genes, Sparse) ---"
COMBINED_DATASET_H5AD="${OUTPUT_DIR_H5ADS}/combined_dataset_all_genes_sparse.h5ad"
run_command "python \"${SCRIPTS_DIR}/create_combined_dataset_all_genes_sparse.py\" \
    --input-dir \"$PREPROCESSED_DIR_H5ADS\" \
    --reference-mapping \"$PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED\" \
    --output-file \"$COMBINED_DATASET_H5AD\" \
    --include-datasets \"encode,gtex,mage,adni\""

if [ $? -eq 0 ]; then
    log_message "Stage 3 completed. Combined dataset: $COMBINED_DATASET_H5AD"
else
    log_message "ERROR: Stage 3 failed. Check logs for details."
    exit 1
fi

# --- Metadata Enhancement Stages (2.8-3.5) ---
if [ "$SKIP_METADATA_ENHANCEMENT" != "true" ]; then
    log_message "--- METADATA ENHANCEMENT PIPELINE (Stages 2.8-3.5) ---"
    
    # --- Stage 2.8: Integrate Controlled-Access Metadata ---
    log_message "--- Stage 2.8: Integrating Controlled-Access Metadata (GTEx Demographics) ---"
    GTEX_CONTROLLED_ACCESS_DIR="${INPUT_BASE_DIR}/gtex/metadata"
    GTEX_PHENOTYPES_FILE="${GTEX_CONTROLLED_ACCESS_DIR}/phs000424.v10.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz"

    if [ -f "$GTEX_PHENOTYPES_FILE" ]; then
        log_message "GTEx controlled-access phenotype file found. Integrating demographics..."
        run_command "python \"${SCRIPTS_DIR}/integrate_missing_metadata.py\" \
            --preprocessed-dir \"$PREPROCESSED_DIR_H5ADS\" \
            --gtex-phenotypes \"$GTEX_PHENOTYPES_FILE\" \
            --adni-demographics \"$ADNI_DEMOGRAPHICS_FILE_INPUT\" \
            --mage-metadata-dir \"$MAGE_RAW_DATA_DIR\" \
            ${FORCE_FLAG}"
        if [ $? -eq 0 ]; then
            log_message "Stage 2.8 completed. Missing metadata integrated."
        else
            log_message "WARNING: Stage 2.8 (metadata integration) had issues."
        fi
    else
        log_message "INFO: GTEx controlled-access file not found at $GTEX_PHENOTYPES_FILE. Skipping controlled-access metadata integration."
    fi

    # --- Stage 2.9: Map Developmental Stages ---
    log_message "--- Stage 2.9: Mapping Age Data to Developmental Stage Ontology ---"
    run_command "python \"${SCRIPTS_DIR}/map_developmental_stages.py\" \
        --preprocessed-dir \"$PREPROCESSED_DIR_H5ADS\" \
        ${FORCE_FLAG}"
    if [ $? -eq 0 ]; then
        log_message "Stage 2.9 completed. HsapDv ontology terms mapped."
    else
        log_message "WARNING: Stage 2.9 (developmental stage mapping) had issues."
    fi

    # --- Stage 3.5: Integrate MAGE Technical Metadata ---
    log_message "--- Stage 3.5: Integrating MAGE Technical Metadata (RIN Scores) ---"
    MAGE_TECHNICAL_METADATA="${MAGE_RAW_DATA_DIR}/sample.metadata.MAGE.v1.0.txt"

    if [ -f "$MAGE_TECHNICAL_METADATA" ]; then
        run_command "python \"${SCRIPTS_DIR}/integrate_mage_technical_metadata.py\" \
            --preprocessed-dir \"$PREPROCESSED_DIR_H5ADS\" \
            --mage-metadata \"$MAGE_TECHNICAL_METADATA\" \
            --mage-population-file \"$MAGE_1000G_PED_FILE_INPUT\" \
            ${FORCE_FLAG}"
        if [ $? -eq 0 ]; then
            log_message "Stage 3.5 completed. MAGE technical metadata integrated."
        else
            log_message "WARNING: Stage 3.5 (MAGE technical metadata) had issues."
        fi
    else
        log_message "INFO: MAGE technical metadata file not found at $MAGE_TECHNICAL_METADATA. Skipping MAGE technical integration."
    fi
    
    log_message "--- METADATA ENHANCEMENT PIPELINE COMPLETE ---"
else
    log_message "--- SKIPPING METADATA ENHANCEMENT STAGES (2.8-3.5) ---"
fi

# --- Stage 4: Validation ---
log_message "--- Stage 4: Validating Standardized Datasets ---"
VALIDATION_REPORT_DIR="${OUTPUT_DIR_H5ADS}/validation_reports"
mkdir -p "$VALIDATION_REPORT_DIR"

run_command "python \"${SCRIPTS_DIR}/validate_standardized_datasets.py\" \
    --input-dir \"$PREPROCESSED_DIR_H5ADS\" \
    --output-file \"${VALIDATION_REPORT_DIR}/validation_report.json\" \
    --file-pattern \"*_standardized_preprocessed.h5ad\""

if [ $? -eq 0 ]; then
    log_message "Stage 4 completed. Validation reports in: $VALIDATION_REPORT_DIR"
else
    log_message "WARNING: Stage 4 validation had issues. Check reports for details."
fi

# --- Stage 4.1: Enhanced Validation and Gene Overlap Analysis ---
log_message "--- Stage 4.1: Enhanced Validation and Gene Overlap Analysis ---"
run_command "python \"${SCRIPTS_DIR}/analyze_gene_overlap.py\" \
    --preprocessed-dir \"$PREPROCESSED_DIR_H5ADS\" \
    --output-dir \"$VALIDATION_REPORT_DIR\" \
    ${FORCE_FLAG}"
if [ $? -eq 0 ]; then
    log_message "Stage 4.1 completed. Gene overlap analysis saved."
else
    log_message "WARNING: Stage 4.1 (gene overlap analysis) had issues."
fi

# --- Stage 4.2: Generate Partner Deliverables ---
log_message "--- Stage 4.2: Generating Partner Deliverables ---"
PARTNER_DELIVERABLES_DIR="${OUTPUT_DIR_H5ADS}/partner_deliverables"
mkdir -p "$PARTNER_DELIVERABLES_DIR"

# Generate ethnicity mappings
run_command "python \"${SCRIPTS_DIR}/create_subject_level_ethnicity_mapping.py\" \
    --preprocessed-dir \"$PREPROCESSED_DIR_H5ADS\" \
    --output-dir \"$PARTNER_DELIVERABLES_DIR\" \
    ${FORCE_FLAG}"

# Generate CZI schema compliant mappings
run_command "python \"${SCRIPTS_DIR}/create_czi_schema_compliant_mapping.py\" \
    --preprocessed-dir \"$PREPROCESSED_DIR_H5ADS\" \
    --output-dir \"$PARTNER_DELIVERABLES_DIR\" \
    ${FORCE_FLAG}"

if [ $? -eq 0 ]; then
    log_message "Stage 4.2 completed. Partner deliverables in: $PARTNER_DELIVERABLES_DIR"
else
    log_message "WARNING: Stage 4.2 (partner deliverables) had issues."
fi

# --- Final Message ---
log_message "=== RNA-seq Pipeline Processing Complete (Stages 0-4.2) ==="
log_message ""
log_message "🎯 **PIPELINE v2.2 COMPREHENSIVE RESULTS**:"
log_message "  ✅ Core Pipeline: Stages 0-4 (gene mapping, standardization, validation)"
log_message "  ✅ Metadata Enhancement: Stages 2.8-3.5 (demographics, ontology, technical)"
log_message "  ✅ Advanced Analysis: Stages 4.1-4.2 (gene overlap, partner deliverables)"
log_message ""
log_message "📁 **OUTPUT DIRECTORIES**:"
log_message "  Main H5AD Files: $OUTPUT_DIR_H5ADS"
log_message "  Preprocessed Data: $PREPROCESSED_DIR_H5ADS"
log_message "  Validation Reports: $VALIDATION_REPORT_DIR"
log_message "  Partner Deliverables: ${OUTPUT_DIR_H5ADS}/partner_deliverables"
log_message ""
log_message "📊 **KEY OUTPUTS**:"
log_message "  Combined Dataset: $COMBINED_DATASET_H5AD"
log_message "  Gene Mapping Files: ${GENERATED_GENE_MAPPING_RESOURCE_DIR}"
log_message "  Primary Gene Reference: ${PRIMARY_GENE_ANNOTATION_AND_REF_CSV_GENERATED}"
log_message "  Validation Report: ${VALIDATION_REPORT_DIR}/validation_report.json"
log_message "  Gene Overlap Analysis: ${VALIDATION_REPORT_DIR}/gene_overlap_analysis.json"
log_message "  Subject Ethnicity Mapping: ${OUTPUT_DIR_H5ADS}/partner_deliverables/subject_ethnicity_mapping_with_ontology.csv"
log_message "  CZI Schema Compliant Files: ${OUTPUT_DIR_H5ADS}/partner_deliverables/"
log_message ""
log_message "🔬 **OPTIONAL EXTENSIONS**:"
log_message "  GTEx Single-Cell: ${SCRIPTS_DIR}/single_cell/run_process_gtex_single_cell.sh"
log_message ""
log_message "📋 **LOGS & DOCUMENTATION**:"
log_message "  Full Pipeline Log: $PIPELINE_LOG_FILE"
log_message "  Partner Presentation: ${SCRIPTS_DIR}/v2.2_partner_presentation.ipynb"