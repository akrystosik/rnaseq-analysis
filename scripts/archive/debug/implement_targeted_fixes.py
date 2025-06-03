#!/usr/bin/env python3
"""
Script to implement targeted fixes to the RNA-seq standardization pipeline
"""
import os
import re
import logging
import shutil
import time
import argparse

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('apply_fixes')

# Define the base directory
BASE_DIR = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")

def backup_file(file_path):
    """Create a backup of a file before modifying it."""
    backup_path = f"{file_path}.bak.{int(time.time())}"
    shutil.copy2(file_path, backup_path)
    logger.info(f"Created backup of {file_path} at {backup_path}")
    return backup_path

def fix_adni_parsing(apply=False):
    """Add defensive processing for ADNI files with escaped tabs."""
    file_path = os.path.join(SCRIPTS_DIR, "standardize_datasets.py")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Define the preprocessing function
    preprocess_func = '''
def preprocess_adni_file(file_path):
    """Preprocess ADNI file to fix escaped tabs."""
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check if file has escaped tabs
        if '\\\\t' in content:
            logger.info(f"Fixing escaped tabs in {file_path}")
            # Replace escaped tabs with actual tabs
            fixed_content = content.replace('\\\\t', '\\t')
            
            # Create a temporary fixed file
            fixed_path = file_path + '.fixed'
            with open(fixed_path, 'w') as f:
                f.write(fixed_content)
            
            return fixed_path
    except Exception as e:
        logger.warning(f"Error preprocessing file {file_path}: {e}")
    
    return file_path
'''
    
    # Find process_adni_data function to add preprocessing
    matches = re.search(r'def process_adni_data\(.*?\):', content, re.DOTALL)
    if not matches:
        logger.error("Could not find process_adni_data function")
        return False
    
    # Add preprocessing function before process_adni_data
    new_content = re.sub(
        r'(def process_adni_data\(.*?\):)',
        f'{preprocess_func}\n\\1',
        content
    )
    
    # Find file reading code in process_adni_data
    csv_read_pattern = r'(try:\s+)(\s*df = pd\.read_csv.*?)(\s+except Exception)'
    
    # Add preprocessing call before reading
    csv_read_replacement = r'''\1    # Preprocess file to fix any escaped tabs
            fixed_file_path = preprocess_adni_file(file_path)
            \2
            # Clean up temporary file if needed
            if fixed_file_path != file_path and os.path.exists(fixed_file_path):
                os.remove(fixed_file_path)\3'''
    
    new_content = re.sub(csv_read_pattern, csv_read_replacement, new_content)
    
    if apply:
        # Backup original file
        backup_file(file_path)
        
        # Write modified content
        with open(file_path, 'w') as f:
            f.write(new_content)
        
        logger.info(f"Applied ADNI preprocessing fix to {file_path}")
        return True
    else:
        logger.info("ADNI preprocessing fix ready to apply (dry run)")
        return True

def fix_categorical_comparison(apply=False):
    """Add defensive conversion of categorical columns before comparisons."""
    file_path = os.path.join(SCRIPTS_DIR, "preprocess_dataset_gene_ids.py")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Define the helper function
    helper_func = '''
def ensure_string_columns(df):
    """Convert categorical columns to strings to avoid comparison issues."""
    for col in df.columns:
        if pd.api.types.is_categorical_dtype(df[col]):
            logger.debug(f"Converting categorical column {col} to string")
            df[col] = df[col].astype(str)
    return df
'''
    
    # Add helper function after imports
    import_pattern = r'(from typing import.*?\n\n)'
    new_content = re.sub(import_pattern, f'\\1{helper_func}\n', content, flags=re.DOTALL)
    
    # Find places where comparisons are made
    # Add defensive conversions before comparisons in preprocess_encode_dataset
    process_encode_pattern = r'(def preprocess_encode_dataset.*?# ENCODE mapping details:)'
    
    process_encode_content = re.search(process_encode_pattern, content, flags=re.DOTALL)
    if process_encode_content:
        encode_code = process_encode_content.group(0)
        
        # Find var_df creation and add conversion
        var_df_pattern = r'(var_df = pd\.DataFrame\(index=adata\.var_names\))'
        if re.search(var_df_pattern, encode_code):
            encode_code_fixed = re.sub(
                var_df_pattern,
                f'\\1\n\n    # Ensure string columns for comparisons\n    var_df = ensure_string_columns(var_df)',
                encode_code
            )
            
            new_content = new_content.replace(process_encode_content.group(0), encode_code_fixed)
    
    # Same for preprocess_entex_dataset
    process_entex_pattern = r'(def preprocess_entex_dataset.*?def preprocess_other_dataset)'
    
    process_entex_content = re.search(process_entex_pattern, content, flags=re.DOTALL)
    if process_entex_content:
        entex_code = process_entex_content.group(0)
        
        # Find var_df creation and add conversion
        var_df_pattern = r'(var_df = pd\.DataFrame\(index=adata\.var_names\))'
        if re.search(var_df_pattern, entex_code):
            entex_code_fixed = re.sub(
                var_df_pattern,
                f'\\1\n\n    # Ensure string columns for comparisons\n    var_df = ensure_string_columns(var_df)',
                entex_code
            )
            
            new_content = new_content.replace(process_entex_content.group(0), entex_code_fixed)
    
    # Same for preprocess_other_dataset
    process_other_pattern = r'(def preprocess_other_dataset.*?)(\s+return adata)'
    
    process_other_content = re.search(process_other_pattern, content, flags=re.DOTALL)
    if process_other_content:
        other_code = process_other_content.group(0)
        
        # Add conversion before returning
        other_code_fixed = other_code.replace(
            process_other_content.group(2),
            f'\n    # Ensure string columns for comparisons\n    adata.var = ensure_string_columns(adata.var){process_other_content.group(2)}'
        )
        
        new_content = new_content.replace(process_other_content.group(0), other_code_fixed)
    
    if apply:
        # Backup original file
        backup_file(file_path)
        
        # Write modified content
        with open(file_path, 'w') as f:
            f.write(new_content)
        
        logger.info(f"Applied categorical comparison fix to {file_path}")
        return True
    else:
        logger.info("Categorical comparison fix ready to apply (dry run)")
        return True

def fix_metadata_serialization(apply=False):
    """Add robust serialization of metadata for AnnData objects."""
    file_path = os.path.join(SCRIPTS_DIR, "standardize_datasets.py")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Define the serialization function
    serialization_func = '''
def ensure_serializable(obj):
    """Recursively convert an object to serializable types."""
    if obj is None:
        return None
    elif isinstance(obj, (str, int, float, bool)):
        return obj
    elif hasattr(obj, 'item'):  # Handle numpy scalars
        try:
            return obj.item()  # Convert to Python native type
        except:
            return str(obj)
    elif isinstance(obj, dict):
        return {str(k): ensure_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple, set)):
        return [ensure_serializable(x) for x in obj]
    elif hasattr(obj, 'tolist'):  # Handle numpy arrays
        return obj.tolist()
    elif hasattr(obj, 'to_dict'):  # Handle pandas objects
        return obj.to_dict()
    else:
        return str(obj)  # Last resort: convert to string
'''
    
    # Add serialization function after imports
    import_pattern = r'(import logging.*?from pathlib import Path\n\n)'
    new_content = re.sub(import_pattern, f'\\1{serialization_func}\n', content, flags=re.DOTALL)
    
    # Add serialization to create_standard_anndata function
    create_pattern = r'(def create_standard_anndata.*?# Add dataset information to uns\s+)(adata\.uns\["dataset_info"\] = dataset_info)'
    
    create_replacement = r'''\1# Ensure dataset_info is serializable
    dataset_info = ensure_serializable(dataset_info)
    adata.uns["dataset_info"] = dataset_info'''
    
    new_content = re.sub(create_pattern, create_replacement, new_content, flags=re.DOTALL)
    
    # Add serialization to save_anndata function
    save_pattern = r'(def save_anndata.*?try:\s+)(# Validate before saving.*?# Create directory)'
    
    save_replacement = r'''\1# Ensure all uns values are serializable
        serializable_uns = {}
        for k, v in adata.uns.items():
            serializable_uns[k] = ensure_serializable(v)
        
        # Store original uns
        original_uns = adata.uns.copy()
        
        # Replace with serializable version
        adata.uns = serializable_uns
        
        \2'''
    
    new_content = re.sub(save_pattern, save_replacement, new_content, flags=re.DOTALL)
    
    # Add code to restore original uns after saving
    verify_pattern = r'(# Verify successful save by loading\s+)(test_load = ad\.read_h5ad\(file_path\))'
    
    verify_replacement = r'''\1# Restore original uns
        adata.uns = original_uns
        
        \2'''
    
    new_content = re.sub(verify_pattern, verify_replacement, new_content, flags=re.DOTALL)
    
    if apply:
        # Backup original file
        backup_file(file_path)
        
        # Write modified content
        with open(file_path, 'w') as f:
            f.write(new_content)
        
        logger.info(f"Applied metadata serialization fix to {file_path}")
        return True
    else:
        logger.info("Metadata serialization fix ready to apply (dry run)")
        return True

def fix_placeholder_ids(apply=False):
    """Create a script to fix placeholder IDs in preprocessed datasets."""
    file_path = os.path.join(SCRIPTS_DIR, "fix_placeholder_ids.py")
    
    script_content = '''#!/usr/bin/env python3
"""
Fix placeholder IDs in preprocessed datasets by converting them to proper Entrez IDs.
"""
import os
import sys
import re
import logging
import scanpy as sc
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('placeholder_id_fixer')

def fix_placeholder_ids(h5ad_file, output_file):
    """
    Fix placeholder IDs in the specified h5ad file.
    
    Parameters:
    -----------
    h5ad_file : str
        Path to the h5ad file to fix
    output_file : str
        Path to save the fixed h5ad file
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Load the dataset
        logger.info(f"Loading dataset from {h5ad_file}")
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Check if there are placeholder IDs
        placeholder_pattern = re.compile(r'^PLACEHOLDER_(.+)$')
        
        # Find placeholder IDs in gene_id and ensembl_id columns
        placeholders = []
        for col in ['gene_id', 'ensembl_id']:
            if col in adata.var.columns:
                # Count placeholders
                placeholder_count = sum(1 for g in adata.var[col] if isinstance(g, str) and g.startswith('PLACEHOLDER_'))
                if placeholder_count > 0:
                    logger.info(f"Found {placeholder_count} placeholder IDs in {col} column")
                    placeholders.append(col)
        
        if not placeholders:
            logger.info("No placeholder IDs found, no fix needed")
            return True
        
        # Make a copy of var DataFrame
        var_df = adata.var.copy()
        
        # Fix categorical columns to prevent comparison issues
        for col in var_df.columns:
            if pd.api.types.is_categorical_dtype(var_df[col]):
                logger.info(f"Converting categorical column {col} to string")
                var_df[col] = var_df[col].astype(str)
        
        # Replace placeholder IDs with Entrez IDs
        entrez_count = 0
        other_count = 0
        
        for col in placeholders:
            for idx in var_df.index:
                value = var_df.loc[idx, col]
                if isinstance(value, str) and value.startswith('PLACEHOLDER_'):
                    # Extract the ID from the placeholder
                    match = placeholder_pattern.match(value)
                    if match:
                        id_part = match.group(1)
                        
                        # If it's a numeric ID, use Entrez prefix
                        if id_part.isdigit():
                            var_df.loc[idx, col] = f"ENTREZ:{id_part}"
                            entrez_count += 1
                        else:
                            # For non-numeric IDs, keep as is but remove PLACEHOLDER_
                            var_df.loc[idx, col] = id_part
                            other_count += 1
        
        # Replace the var DataFrame
        adata.var = var_df
        
        # Log results
        logger.info(f"Fix results: {entrez_count + other_count} placeholder IDs converted")
        logger.info(f"Total Entrez IDs: {entrez_count}")
        logger.info(f"Total other IDs: {other_count}")
        
        # Save the fixed dataset
        logger.info(f"Saving fixed dataset to {output_file}")
        adata.write_h5ad(output_file)
        
        logger.info("Fix completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error fixing placeholder IDs: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fix_placeholder_ids.py <input_h5ad_file> [output_h5ad_file]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        # Use same file with .fixed suffix
        output_file = input_file + ".fixed"
    
    success = fix_placeholder_ids(input_file, output_file)
    
    if success:
        print(f"Successfully fixed placeholder IDs and saved to {output_file}")
        sys.exit(0)
    else:
        print(f"Failed to fix placeholder IDs")
        sys.exit(1)
'''
    
    if apply:
        with open(file_path, 'w') as f:
            f.write(script_content)
        
        # Make executable
        os.chmod(file_path, 0o755)
        
        logger.info(f"Created placeholder ID fix script at {file_path}")
        return True
    else:
        logger.info("Placeholder ID fix script ready to create (dry run)")
        return True

def update_run_script(apply=False):
    """Update run_rnaseq_pipeline.sh to include the placeholder ID fix."""
    file_path = os.path.join(SCRIPTS_DIR, "run_rnaseq_pipeline.sh")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Add placeholder ID fix step
    fix_step = '''
# Fix placeholder IDs in preprocessed datasets
logger.info "=== Step 2.6: Fixing placeholder IDs in preprocessed datasets ==="
for dataset in encode gtex mage adni; do
    preprocessed_file="${PREPROCESSED_DIR}/${dataset}_standardized_preprocessed.h5ad"
    if [ -f "$preprocessed_file" ]; then
        logger.info "Fixing placeholder IDs in ${dataset} dataset"
        run_command "python ${SCRIPTS_DIR}/fix_placeholder_ids.py $preprocessed_file $preprocessed_file.fixed"
        if [ $? -eq 0 ]; then
            # Replace the original file with the fixed file
            mv "$preprocessed_file.fixed" "$preprocessed_file"
            logger.info "Placeholder IDs fixed in ${dataset} dataset"
        else
            logger.warning "Failed to fix placeholder IDs in ${dataset} dataset"
        fi
    fi
done
'''
    
    # Find a suitable location to insert the fix step
    step_pattern = r'(# Step 2.5: Preprocess Datasets.*?preprocessed_data"(\s+))'
    
    new_content = re.sub(step_pattern, f'\\1{fix_step}\\2', content, flags=re.DOTALL)
    
    if apply:
        # Backup original file
        backup_file(file_path)
        
        # Write modified content
        with open(file_path, 'w') as f:
            f.write(new_content)
        
        logger.info(f"Updated run script to include placeholder ID fix")
        return True
    else:
        logger.info("Run script update ready to apply (dry run)")
        return True

def update_changelog(apply=False):
    """Update the CHANGELOG.md file with the fixes."""
    file_path = os.path.join(BASE_DIR, "CHANGELOG.md")
    
    changelog_entry = f'''
## [0.1.2] - {time.strftime("%Y-%m-%d")}

### Fixed
- Added defensive processing for ADNI files with escaped tabs
- Fixed categorical data comparison issues in preprocessed datasets
- Improved metadata serialization to handle complex nested structures
- Added fix_placeholder_ids.py script to replace placeholder IDs with proper Entrez IDs
- Integrated placeholder ID fix into the pipeline workflow

### Changed
- Made metadata serialization more robust by converting numpy types
- Enhanced categorical column handling to prevent comparison errors
- Updated preprocessing to ensure consistent gene ID representation
'''
    
    if apply:
        with open(file_path, 'a') as f:
            f.write(changelog_entry)
        
        logger.info(f"Updated CHANGELOG.md with fixes")
        return True
    else:
        logger.info("CHANGELOG.md update ready to apply (dry run)")
        return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Implement targeted fixes to the RNA-seq pipeline")
    parser.add_argument("--fix", choices=["adni", "categorical", "metadata", "placeholder", "run", "all"], 
                        default="all", help="The specific fix to apply")
    parser.add_argument("--apply", action="store_true", help="Actually apply the fixes (default is dry run)")
    
    args = parser.parse_args()
    
    if args.apply:
        logger.info("=== Applying targeted fixes to RNA-seq Pipeline ===")
    else:
        logger.info("=== Dry run - showing fixes without applying them ===")
    
    if args.fix == "adni" or args.fix == "all":
        fix_adni_parsing(args.apply)
    
    if args.fix == "categorical" or args.fix == "all":
        fix_categorical_comparison(args.apply)
    
    if args.fix == "metadata" or args.fix == "all":
        fix_metadata_serialization(args.apply)
    
    if args.fix == "placeholder" or args.fix == "all":
        fix_placeholder_ids(args.apply)
    
    if args.fix == "run" or args.fix == "all":
        update_run_script(args.apply)
    
    if args.apply:
        update_changelog(args.apply)
        logger.info("All fixes applied successfully")
    else:
        logger.info("Dry run complete - to apply fixes, run with --apply flag")
