#!/usr/bin/env python3
"""
MVP Check Script for Standardized RNA-seq Datasets (v8 - Added Head/Tissue Output)

Checks essential requirements for MVP:
1. adata.var_names (index) must have < 5% non-Ensembl IDs. (RELAXED)
2. adata.var['gene_id'] must exist and match the index.
3. adata.obs['tissue'] must exist (reports missing/empty values).
4. adata.obs['subject_id'] must exist (reports missing/empty values).
5. adata.obs['dataset'] must exist and contain ONLY the dataset name (case-insensitive).
"""

import scanpy as sc
import pandas as pd
import os
import sys
import glob
import logging
import re # Import regex
import numpy as np # Import numpy for isnan checks and array operations

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('mvp_checker')

# --- Helper to update deprecated pandas check ---
def is_categorical(dtype):
    """Check if dtype is categorical, handling deprecation."""
    if hasattr(pd.api.types, 'is_categorical_dtype'):
        # Suppress the warning locally if possible, or just use the check
        with pd.option_context('mode.chained_assignment', None): # Context might not suppress this specific warning
             # Use isinstance for forward compatibility
             return isinstance(dtype, pd.CategoricalDtype)
    else: # Fallback
        return isinstance(dtype, pd.CategoricalDtype)
# --- End helper ---


def check_adata_mvp(file_path, dataset_name):
    """Checks a single AnnData file for RELAXED MVP requirements."""
    logger.info(f"--- Checking {dataset_name} ({os.path.basename(file_path)}) ---")
    required_obs_cols = ['tissue', 'subject_id', 'dataset']
    passed = True # Assume passing initially
    messages = []
    diagnostic_outputs = {} # To store head outputs and unique tissues
    dataset_name_lower = dataset_name.lower() # Expected value (lowercase)
    non_ensembl_percentage = -1.0 # Initialize in case of error

    # Define the threshold for non-Ensembl IDs
    NON_ENSEMBL_THRESHOLD = 5.0 # Percentage

    try:
        adata = sc.read_h5ad(file_path)
        messages.append(f"  Shape: {adata.shape}")
        total_genes = adata.n_vars

        # --- Capture Diagnostic Outputs ---
        diagnostic_outputs['obs_head'] = adata.obs.head().to_string()
        diagnostic_outputs['var_head'] = adata.var.head().to_string()
        if 'tissue' in adata.obs.columns:
            unique_tissues = adata.obs['tissue'].unique()
            if is_categorical(unique_tissues.dtype): # Handle categorical output
                unique_tissues = unique_tissues.tolist()
            diagnostic_outputs['unique_tissues'] = sorted([str(t) for t in unique_tissues if pd.notna(t)]) # Sort and stringify non-NA
        else:
            diagnostic_outputs['unique_tissues'] = "Column 'tissue' not found."
        # --- End Capture ---


        # --- Requirement 1: Var Index must have < THRESHOLD % non-Ensembl IDs ---
        var_names_sample = adata.var_names[:10].tolist()
        messages.append(f"  Sample var_names (index): {var_names_sample}")

        ensembl_pattern = re.compile(r"ENSG\d+(\.\d+)?$")
        non_ensembl_count = 0
        non_ensembl_examples = []
        for x in adata.var_names:
            x_str = str(x)
            if not pd.isna(x) and not ensembl_pattern.fullmatch(x_str):
                 non_ensembl_count += 1
                 if len(non_ensembl_examples) < 5:
                     non_ensembl_examples.append(x_str)

        if total_genes > 0:
            non_ensembl_percentage = (non_ensembl_count / total_genes) * 100
        else:
            non_ensembl_percentage = 0.0

        if non_ensembl_percentage < NON_ENSEMBL_THRESHOLD:
            pass_msg = f"  [PASS] var_names index format: {non_ensembl_percentage:.2f}% non-Ensembl IDs (Threshold: < {NON_ENSEMBL_THRESHOLD}%)."
            if non_ensembl_count > 0:
                pass_msg += f" Examples: {non_ensembl_examples}"
            messages.append(pass_msg)
        else:
            messages.append(f"  [FAIL] var_names index format: {non_ensembl_percentage:.2f}% non-Ensembl IDs (Threshold: < {NON_ENSEMBL_THRESHOLD}%). Examples: {non_ensembl_examples}")
            passed = False

        # --- Requirement 2: adata.var['gene_id'] must exist and match index ---
        if 'gene_id' not in adata.var.columns:
            messages.append("  [FAIL] Missing required var column: 'gene_id'")
            passed = False
        elif not np.array_equal(adata.var['gene_id'].astype(str).values, adata.var_names.astype(str).values):
             mismatch_indices = np.where(adata.var['gene_id'].astype(str).values != adata.var_names.astype(str).values)[0]
             if len(mismatch_indices) > 0:
                  first_mismatch_loc = mismatch_indices[0]
                  mismatch_example_idx_val = adata.var_names[first_mismatch_loc]
                  mismatch_gene_id_val = adata.var['gene_id'].iloc[first_mismatch_loc]
                  messages.append(f"  [FAIL] adata.var['gene_id'] does not match adata.var_names index (e.g., at index position {first_mismatch_loc}, index='{mismatch_example_idx_val}', gene_id='{mismatch_gene_id_val}')")
             else:
                  messages.append(f"  [FAIL] adata.var['gene_id'] does not match adata.var_names index (mismatch detected but example could not be found).")
             passed = False
        else:
             messages.append("  [PASS] 'gene_id' column exists and matches index.")

        # --- Requirement 3, 4, 5: Required obs columns must exist ---
        missing_obs = [col for col in required_obs_cols if col not in adata.obs.columns]
        if missing_obs:
            messages.append(f"  [FAIL] Missing required obs columns: {missing_obs}")
            passed = False
        else:
            messages.append(f"  Required obs columns present: {required_obs_cols}")

            # --- Check dataset column content (Case-Insensitive) ---
            if 'dataset' in adata.obs.columns:
                 col_series = adata.obs['dataset']
                 if is_categorical(col_series.dtype):
                     unique_datasets_in_col = col_series.astype(str).str.lower().unique().tolist()
                 else:
                     unique_datasets_in_col = col_series.astype(str).str.lower().unique()

                 is_correct_dataset = (len(unique_datasets_in_col) == 1 and unique_datasets_in_col[0] == dataset_name_lower)

                 if not is_correct_dataset:
                      original_unique_vals = col_series.unique()
                      if is_categorical(original_unique_vals.dtype):
                          original_unique_vals = original_unique_vals.tolist()
                      messages.append(f"  [FAIL] 'dataset' column value mismatch. Expected only '{dataset_name_lower}', found: {original_unique_vals}")
                      passed = False
                 else:
                      messages.append(f"  [PASS] 'dataset' column contains correct value (case-insensitive).")

            # --- Check for NAs/empties in tissue and subject_id (Detailed Info) ---
            for col in ['tissue', 'subject_id']:
                 if col in adata.obs.columns:
                    col_series = adata.obs[col]
                    is_na_mask_np = col_series.isna().to_numpy(dtype=bool)
                    is_empty_series = col_series.apply(
                        lambda x: isinstance(x, str) and x.strip() == ''
                    )
                    is_empty_mask_np = is_empty_series.to_numpy(dtype=bool)
                    missing_mask_np = np.logical_or(is_na_mask_np, is_empty_mask_np)

                    try:
                        total_missing = np.sum(missing_mask_np)
                    except Exception as sum_err:
                        logger.warning(f"NumPy sum failed for combined mask in column '{col}': {sum_err}. Falling back to manual count.")
                        total_missing = sum(1 for x in missing_mask_np if x is True)

                    if total_missing > 0:
                         missing_values_unique = pd.unique(col_series.iloc[missing_mask_np])
                         missing_examples = [repr(v) for v in missing_values_unique[:5]]
                         messages.append(f"  [INFO] Obs column '{col}' has {total_missing}/{adata.n_obs} missing/empty values. Examples found: {', '.join(missing_examples)}")
                    else:
                         messages.append(f"  [INFO] Obs column '{col}' has no missing (NaN) or empty values.")
                 else:
                     messages.append(f"  [INFO] Required obs column '{col}' is missing (already marked as FAIL).")

    except Exception as e:
        error_msg = f"  [ERROR] Could not complete check for file due to error: {e}"
        messages.append(error_msg)
        logger.error(error_msg, exc_info=True) # Log full traceback
        passed = False

    final_status = "PASSED (MVP)" if passed else "FAILED (MVP)"
    logger.info(f"  Status: {final_status}")
    for msg in messages:
        if "[ERROR]" not in msg:
            logger.info(msg)

    # --- Console Output Section ---
    print(f"\n=== {dataset_name} MVP Check ===")
    print(f"File: {file_path}")
    print(f"Status: {final_status}")
    if non_ensembl_percentage >= 0:
         print(f"INFO: Non-Ensembl Gene IDs: {non_ensembl_percentage:.2f}% (Threshold: < {NON_ENSEMBL_THRESHOLD}%)")

    # Print diagnostic heads and tissues
    print("\n--- obs.head() ---")
    print(diagnostic_outputs.get('obs_head', 'N/A'))
    print("\n--- var.head() ---")
    print(diagnostic_outputs.get('var_head', 'N/A'))
    print("\n--- Unique Tissues ---")
    print(diagnostic_outputs.get('unique_tissues', 'N/A'))
    print("\n--- Check Messages ---")

    # Print check messages
    for msg in messages[1:]: # Skip shape message for console
        cleaned_msg = msg.strip().replace("[PASS]","").replace("[FAIL]","FAIL:").replace("[WARNING]","WARN:").replace("[INFO]","INFO:").replace("[ERROR]", "ERROR:")
        if cleaned_msg.startswith("Examples:") or cleaned_msg.startswith("Non-Ensembl examples:"):
             print(f"    {cleaned_msg}")
        elif "FAIL:" in cleaned_msg or "ERROR:" in cleaned_msg:
            print(f"  {cleaned_msg}")
        else:
             print(cleaned_msg)
    print("="*(len(dataset_name)+14))
    # --- End Console Output Section ---

    return passed


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_mvp.py <preprocessed_data_directory>")
        sys.exit(1)

    preprocessed_dir = sys.argv[1]

    if not os.path.isdir(preprocessed_dir):
        logger.error(f"Directory not found: {preprocessed_dir}")
        sys.exit(1)

    logger.info(f"Checking MVP requirements in directory: {preprocessed_dir}")

    # Find preprocessed files
    file_pattern = os.path.join(preprocessed_dir, "*_standardized_preprocessed.h5ad")
    preprocessed_files = glob.glob(file_pattern)

    if not preprocessed_files:
        logger.error(f"No preprocessed files found matching pattern: {file_pattern}")
        sys.exit(1)

    logger.info(f"Found {len(preprocessed_files)} preprocessed files to check.")

    all_passed = True
    found_any = False

    for file_path in preprocessed_files:
        dataset_name = os.path.basename(file_path).split('_')[0]
        if not dataset_name:
            logger.warning(f"Could not determine dataset name from file: {file_path}. Skipping.")
            continue

        found_any = True
        if not check_adata_mvp(file_path, dataset_name.upper()):
            all_passed = False

    logger.info("--- MVP Check Summary ---")
    if not found_any:
         logger.error("No datasets could be checked in the specified directory.")
         sys.exit(1)

    if all_passed:
        logger.info("All checked datasets PASSED the RELAXED MVP requirements.")
        print("\nOverall MVP Status (Relaxed): PASSED")
        sys.exit(0)
    else:
        logger.error("One or more datasets FAILED the RELAXED MVP requirements. See details above.")
        print("\nOverall MVP Status (Relaxed): FAILED")
        sys.exit(1)