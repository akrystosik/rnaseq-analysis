#!/bin/bash
# ENTEx Data Integration Pipeline
# This script orchestrates the downloading and processing of ENTEx data
# for integration into the gene expression analysis pipeline.

set -e
cd /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts
# Default values
INPUT_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/entex/download_list.txt"
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/entex"
METADATA_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/metadata"
THREADS=4
STANDARDIZATION_SCRIPT="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py"
STANDARDIZED_DATA_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data"

# Check for required Python packages
echo "Checking for required Python packages..."
python3 -c "
required_packages = ['pandas', 'matplotlib', 'seaborn', 'requests', 'anndata']
missing_packages = []
for package in required_packages:
    try:
        __import__(package)
        print(f'✓ {package} is installed')
    except ImportError:
        missing_packages.append(package)
        print(f'✗ {package} is NOT installed')
if missing_packages:
    print(f'\nERROR: Missing required packages: {missing_packages}')
    print('Please run setup_workspace.sh to install dependencies')
    exit(1)
else:
    print('\nAll required packages are installed!')
"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --input|-i)
      INPUT_FILE="$2"
      shift 2
      ;;
    --output-dir|-o)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --metadata-dir|-m)
      METADATA_DIR="$2"
      shift 2
      ;;
    --threads|-t)
      THREADS="$2"
      shift 2
      ;;
    --standardization-script|-s)
      STANDARDIZATION_SCRIPT="$2"
      shift 2
      ;;
    --standardized-data-dir|-d)
      STANDARDIZED_DATA_DIR="$2"
      shift 2
      ;;
    --help|-h)
      echo "Usage: $0 [OPTIONS]"
      echo "Options:"
      echo "  --input, -i FILE           Input file containing URLs (default: paste.txt)"
      echo "  --output-dir, -o DIR       Output directory for downloaded files"
      echo "  --metadata-dir, -m DIR     Directory for metadata files"
      echo "  --threads, -t NUM          Number of download threads (default: 4)"
      echo "  --standardization-script, -s FILE  Path to standardization script"
      echo "  --standardized-data-dir, -d DIR    Directory for standardized data"
      echo "  --help, -h                 Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$METADATA_DIR"
mkdir -p "$STANDARDIZED_DATA_DIR"

echo "=== ENTEx Data Integration Pipeline ==="
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Metadata directory: $METADATA_DIR"
echo "Threads: $THREADS"
echo "Standardization script: $STANDARDIZATION_SCRIPT"
echo "Standardized data directory: $STANDARDIZED_DATA_DIR"
echo ""

# Step 1: Download ENTEx files
echo "=== Step 1: Downloading ENTEx files ==="
RAW_METADATA="$METADATA_DIR/entex_raw_metadata.json"
python3 download_entex_files.py \
  --input "$INPUT_FILE" \
  --output-dir "$OUTPUT_DIR" \
  --threads "$THREADS" \
  --metadata-file "$RAW_METADATA"

# Step 2: Generate standardized metadata
echo "=== Step 2: Generating standardized metadata ==="
EXISTING_METADATA="$METADATA_DIR/entex_manual_metadata.json"
STANDARDIZED_METADATA="$METADATA_DIR/entex_metadata.json"

# Check if existing metadata exists
if [ -f "$EXISTING_METADATA" ]; then
  echo "Found existing metadata at $EXISTING_METADATA"
  python3 generate_entex_metadata.py \
    --input-dir "$OUTPUT_DIR" \
    --metadata-file "$RAW_METADATA" \
    --output-file "$STANDARDIZED_METADATA" \
    --existing-metadata "$EXISTING_METADATA"
else
  echo "No existing metadata found, creating new file"
  python3 generate_entex_metadata.py \
    --input-dir "$OUTPUT_DIR" \
    --metadata-file "$RAW_METADATA" \
    --output-file "$STANDARDIZED_METADATA"
fi

# Step 3: Run standardization pipeline with updated parameters
echo "=== Step 3: Running standardization pipeline ==="
if [ -f "$STANDARDIZATION_SCRIPT" ]; then
  # Check if the script has the new --metadata and --dataset-type parameters
  if grep -q "metadata" "$STANDARDIZATION_SCRIPT" && grep -q "dataset-type" "$STANDARDIZATION_SCRIPT"; then
    echo "Using updated standardization script with direct ENTEx support"
    python3 "$STANDARDIZATION_SCRIPT" \
      --metadata "$STANDARDIZED_METADATA" \
      --dataset-type "entex" \
      --output-dir "$STANDARDIZED_DATA_DIR"
  else
    echo "Using original standardization script with encode-entex-dir parameter"
    python3 "$STANDARDIZATION_SCRIPT" \
      --encode-entex-dir "$OUTPUT_DIR" \
      --output-dir "$STANDARDIZED_DATA_DIR"
  fi
  
  echo "Standardization complete! Data available at $STANDARDIZED_DATA_DIR"
else
  echo "ERROR: Standardization script not found at $STANDARDIZATION_SCRIPT"
  echo "Please check the path to your standardization script"
fi

# Step 4: Verify the integration
echo "=== Step 4: Verifying integration ==="
ENTEX_H5AD="$STANDARDIZED_DATA_DIR/entex_standardized.h5ad"
if [ -f "$ENTEX_H5AD" ]; then
  echo "Verification: ENTEx standardized file created successfully at $ENTEX_H5AD"
  
  # Get basic stats about the h5ad file
  python3 -c "
import sys
try:
    import anndata as ad
    import os
    
    file_path = '$ENTEX_H5AD'
    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
    
    print(f'File size: {file_size_mb:.2f} MB')
    
    adata = ad.read_h5ad(file_path)
    print(f'Samples: {adata.n_obs}')
    print(f'Genes: {adata.n_vars}')
    
    if 'tissue' in adata.obs.columns:
        tissue_counts = adata.obs['tissue'].value_counts()
        print(f'Tissues: {len(tissue_counts)}')
        print('Top tissues:')
        for tissue, count in tissue_counts.head(5).items():
            print(f'  - {tissue}: {count} samples')
    
    if 'subject_id' in adata.obs.columns:
        donor_counts = adata.obs['subject_id'].value_counts()
        print(f'Donors: {len(donor_counts)}')
    
except ImportError:
    print('Error: anndata package not available')
    sys.exit(1)
except Exception as e:
    print(f'Error verifying h5ad file: {e}')
    sys.exit(1)
"
else
  echo "Warning: ENTEx standardized file not found at $ENTEX_H5AD"
  
  # Check for any standardized files
  STANDARDIZED_FILES=$(find "$STANDARDIZED_DATA_DIR" -name "*.h5ad" -or -name "*.h5" | wc -l)
  if [ "$STANDARDIZED_FILES" -gt 0 ]; then
    echo "Found $STANDARDIZED_FILES standardized files in $STANDARDIZED_DATA_DIR"
    find "$STANDARDIZED_DATA_DIR" -name "*.h5ad" -or -name "*.h5" -exec ls -lh {} \;
  else
    echo "No standardized files found. Standardization may have failed."
  fi
fi

echo "=== Pipeline completed successfully ==="
echo "Raw data: $OUTPUT_DIR"
echo "Metadata: $STANDARDIZED_METADATA"
echo "Standardized data: $STANDARDIZED_DATA_DIR"