
#!/bin/bash
# RNA-seq Analysis Pipeline Execution Script
# This script automates the execution of the entire pipeline

set -e  # Exit on error

# Default values
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis"
INCLUDE_ENTEX=false
NUM_DATASETS=""
USE_LOG=false
CORRELATION_TYPE="pearson"
CLUSTER_THRESHOLD=0.9
SKIP_QC=false
DATASETS=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --base-dir)
      BASE_DIR="$2"
      shift 2
      ;;
    --include-entex)
      INCLUDE_ENTEX=true
      shift
      ;;
    --num-datasets)
      NUM_DATASETS="$2"
      shift 2
      ;;
    --datasets)
      DATASETS="$2"
      shift 2
      ;;
    --use-log)
      USE_LOG=true
      shift
      ;;
    --correlation-type)
      CORRELATION_TYPE="$2"
      shift 2
      ;;
    --cluster-threshold)
      CLUSTER_THRESHOLD="$2"
      shift 2
      ;;
    --skip-qc)
      SKIP_QC=true
      shift
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Print configuration
echo "===== RNA-seq Analysis Pipeline ====="
echo "Base Directory: $BASE_DIR"
echo "Include ENTEx: $INCLUDE_ENTEX"
echo "Number of Datasets Limit: $NUM_DATASETS"
echo "Use Log Transform: $USE_LOG"
echo "Correlation Type: $CORRELATION_TYPE"
echo "Cluster Threshold: $CLUSTER_THRESHOLD"
echo "Skip QC: $SKIP_QC"
echo "Specific Datasets: $DATASETS"
echo "======================================"

# Create timestamp for logs
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="$BASE_DIR/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/pipeline_${TIMESTAMP}.log"

# Function to log messages
log() {
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1" | tee -a "$LOG_FILE"
}

# Start logging
log "Starting RNA-seq analysis pipeline"
log "Logs will be saved to $LOG_FILE"

# Step 1: Setup
log "Step 1: Setting up pipeline"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$BASE_DIR/scripts/setup.py" --base-dir "$BASE_DIR" --script-dir "$SCRIPT_DIR" --skip-package-check | tee -a "$LOG_FILE"


# Step 2: Build gene mapping
log "Step 2: Building gene mapping database"
python "$BASE_DIR/scripts/build_gene_mapping.py" --use-biomart | tee -a "$LOG_FILE"

# Step 3: Run analysis
log "Step 3: Running analysis"

# Build command
ANALYSIS_CMD="python $BASE_DIR/scripts/analyze_cell_lines.py"

# Add options
if $INCLUDE_ENTEX; then
  ANALYSIS_CMD="$ANALYSIS_CMD --include-entex"
fi

if [ -n "$NUM_DATASETS" ]; then
  # Get limited number of datasets
  ALL_DATASETS=$(python -c "import json; print(' '.join(json.load(open('$BASE_DIR/config/settings.json'))['cell_lines'].keys()))")
  DATASETS=$(echo $ALL_DATASETS | tr ' ' '\n' | head -n $NUM_DATASETS | tr '\n' ' ')
  ANALYSIS_CMD="$ANALYSIS_CMD --datasets $DATASETS"
elif [ -n "$DATASETS" ]; then
  ANALYSIS_CMD="$ANALYSIS_CMD --datasets $DATASETS"
fi

if $USE_LOG; then
  ANALYSIS_CMD="$ANALYSIS_CMD --use-log"
fi

if [ "$CORRELATION_TYPE" != "pearson" ]; then
  ANALYSIS_CMD="$ANALYSIS_CMD --correlation-type $CORRELATION_TYPE"
fi

if [ "$CLUSTER_THRESHOLD" != "0.9" ]; then
  ANALYSIS_CMD="$ANALYSIS_CMD --cluster-threshold $CLUSTER_THRESHOLD"
fi

if $SKIP_QC; then
  ANALYSIS_CMD="$ANALYSIS_CMD --skip-qc"
fi

# Run analysis
log "Running command: $ANALYSIS_CMD"
eval "$ANALYSIS_CMD" | tee -a "$LOG_FILE"

# Check if analysis was successful
if [ $? -eq 0 ]; then
  log "Analysis completed successfully"
else
  log "Analysis failed"
  exit 1
fi

# Step 4: Print summary and next steps
log "Step 4: Pipeline execution complete"
log "Check analysis results in $BASE_DIR/analysis directory"

# Find the latest analysis directory
LATEST_ANALYSIS=$(ls -td "$BASE_DIR/analysis"/*/ 2>/dev/null | head -n 1)
if [ -n "$LATEST_ANALYSIS" ]; then
  log "Latest analysis results: $LATEST_ANALYSIS"
  
  # Check if visualizations were generated
  VIZ_DIR="$LATEST_ANALYSIS/visualizations"
  if [ -d "$VIZ_DIR" ]; then
    log "Visualizations available in: $VIZ_DIR"
    
    # List generated visualizations
    VIZS=$(find "$VIZ_DIR" -type f -name "*.png" | sort)
    if [ -n "$VIZS" ]; then
      log "Generated visualizations:"
      for VIZ in $VIZS; do
        log "  - $(basename "$VIZ")"
      done
    fi
  fi
fi

log "Pipeline execution completed at $(date)"
echo "===== Pipeline Execution Complete ====="
echo "Log file: $LOG_FILE"
