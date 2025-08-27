#!/bin/bash
# Script to update run_rnaseq_pipeline.sh to use only cell_lines for ENCODE

# Create a backup of the original file
cp /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/run_rnaseq_pipeline.sh /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/run_rnaseq_pipeline.sh.bak

# Create the patch
cat > /tmp/pipeline_encode_fix.patch << 'EOL'
--- run_rnaseq_pipeline.sh	2025-04-30 07:15:00.000000000 +0000
+++ run_rnaseq_pipeline.sh.fixed	2025-04-30 07:15:00.000000000 +0000
@@ -108,7 +108,7 @@
 # Step 1: Initial Data Conversion (existing)
 log_message "=== Stage 1: Initial Data Conversion ==="
 run_command "python ${SCRIPTS_DIR}/standardize_datasets.py \\
-    --encode-dir \"${BASE_DIR}/encode/raw_data\" \\
+    --encode-dir \"${BASE_DIR}/encode/raw_data/cell_lines\" \\
     --encode-entex-dir \"${BASE_DIR}/encode/entex\" \\
     --entex-metadata-file \"${BASE_DIR}/encode/metadata/entex_metadata.json\" \\
     --gtex-file \"${BASE_DIR}/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz\" \\
EOL

# Apply the patch
patch /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/run_rnaseq_pipeline.sh /tmp/pipeline_encode_fix.patch

echo "Pipeline ENCODE path fix applied successfully."