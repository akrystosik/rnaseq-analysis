#!/bin/bash
# Script to create symbolic links with correct filenames

# Set base directory
LATEST_RUN="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250428_234225"

# Create symbolic links
ln -sf "$LATEST_RUN/encode_standardized_v2.h5ad" "$LATEST_RUN/encode_standardized.h5ad"
ln -sf "$LATEST_RUN/gtex_standardized_v2.h5ad" "$LATEST_RUN/gtex_standardized.h5ad"
ln -sf "$LATEST_RUN/mage_standardized_v2.h5ad" "$LATEST_RUN/mage_standardized.h5ad"

echo "Created symbolic links with expected filenames"
echo "Next step: re-run the combined dataset generation script"
