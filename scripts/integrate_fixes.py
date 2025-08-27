#!/usr/bin/env python3
"""
Script to integrate fixes into the RNA-seq standardization pipeline
"""
import os
import sys
import logging
import re
import shutil

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('integrate_fixes')

# Define the base directory
BASE_DIR = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")

def backup_file(file_path):
    """Create a backup of a file before modifying it."""
    backup_path = file_path + ".bak." + str(int(time.time()))
    shutil.copy2(file_path, backup_path)
    logger.info(f"Created backup of {file_path} at {backup_path}")
    return backup_path

def integrate_adni_parsing_fix():
    """Integrate the ADNI file parsing fix into standardize_datasets.py."""
    file_path = os.path.join(SCRIPTS_DIR, "standardize_datasets.py")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    logger.info(f"Adding ADNI file parsing fix to {file_path}")
    
    # Read the file content
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Look for the process_adni_file function to modify
    if "def process_adni_file(" not in content:
        logger.error("Could not find process_adni_file function in standardize_datasets.py")
        return False
    
    # Define the fix code
    fix_code = """
def preprocess_adni_file(file_path):
    \"\"\"Preprocess ADNI file to fix escaped tab characters.\"\"\"
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
    
    return file_path
    """
    
    # Insert the fix function before process_adni_file
    content = re.sub(r'(def process_adni_file\(.*?\):)', fix_code + r'\n\1', content, flags=re.DOTALL)
    
    # Modify the process_adni_file function to use the preprocessing
    process_adni_file_fix = """
            # Preprocess the file to fix escaped tabs
            fixed_file_path = preprocess_adni_file(file_path)
            try:
                # Try reading with tab delimiter first
                df = pd.read_csv(fixed_file_path, sep='\\t')
                # If successful, clean up temporary file if it was created
                if fixed_file_path != file_path and os.path.exists(fixed_file_path):
                    os.remove(fixed_file_path)
            except Exception as e:
                # If tab delimiter fails, try comma
                logger.warning(f"Failed to parse with tab delimiter: {e}")
                df = pd.read_csv(fixed_file_path)
                # Clean up temporary file
                if fixed_file_path != file_path and os.path.exists(fixed_file_path):
                    os.remove(fixed_file_path)
"""
    
    # Replace the old file reading code
    content = re.sub(r'(try:[\s\n]+)(\s*df = pd\.read_csv.*?)(\s*except Exception)', 
                    r'\1' + process_adni_file_fix + r'\3', 
                    content)
    
    # Create backup
    backup_file(file_path)
    
    # Write updated content
    with open(file_path, 'w') as f:
        f.write(content)
    
    logger.info("ADNI file parsing fix integrated")
    return True

def integrate_categorical_fix():
    """Integrate the categorical data comparison fix into preprocess_dataset_gene_ids.py."""
    file_path = os.path.join(SCRIPTS_DIR, "preprocess_dataset_gene_ids.py")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    logger.info(f"Adding categorical data comparison fix to {file_path}")
    
    # Read the file content
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Define the fix function
    fix_function = """
def fix_categorical_columns(df):
    \"\"\"
    Fix categorical columns by converting them to string type.
    This prevents the "Categoricals can only be compared if 'categories' are the same" error.
    \"\"\"
    for col in df.columns:
        if pd.api.types.is_categorical_dtype(df[col]):
            logger.debug(f"Converting categorical column {col} to string")
            df[col] = df[col].astype(str)
    return df
"""
    
    # Insert the fix function at the end of the imports
    content = re.sub(r'(import .*?\n\n)', r'\1\n' + fix_function, content, flags=re.DOTALL)
    
    # Find all places where var_df is used in comparisons
    # Add the fix before var operations
    var_operations = [
        r'var_df\["gene_id"\] == original_id',
        r'var_df\["gene_id"\]',
        r'var_df\.loc\[.*?\]'
    ]
    
    for pattern in var_operations:
        if re.search(pattern, content):
            # Add the fix before first occurrence
            fix_call = "    # Fix categorical columns to prevent comparison issues\n    var_df = fix_categorical_columns(var_df)\n\n"
            # Find the appropriate place to insert the fix (after var_df is defined)
            var_def_pattern = r'(var_df = .*?\n)'
            content = re.sub(var_def_pattern, r'\1\n' + fix_call, content, count=1)
            break
    
    # Create backup
    backup_file(file_path)
    
    # Write updated content
    with open(file_path, 'w') as f:
        f.write(content)
    
    logger.info("Categorical data comparison fix integrated")
    return True

def integrate_metadata_serialization_fix():
    """Integrate the metadata serialization fix into standardize_datasets.py."""
    file_path = os.path.join(SCRIPTS_DIR, "standardize_datasets.py")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    logger.info(f"Adding metadata serialization fix to {file_path}")
    
    # Read the file content
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Look for the create_standard_anndata and save_anndata functions
    if "def create_standard_anndata(" not in content or "def save_anndata(" not in content:
        logger.error("Could not find required functions in standardize_datasets.py")
        return False
    
    # Define the safe serialization function
    safe_serialization_function = """
def safe_str_serialization(uns_dict):
    \"\"\"
    Convert all values in a dictionary to serializable types for AnnData.uns
    \"\"\"
    if uns_dict is None:
        return {}
    
    result = {}
    for key, value in uns_dict.items():
        # Convert key to string just to be safe
        str_key = str(key)
        
        # Handle different value types
        if value is None:
            result[str_key] = ''
        elif isinstance(value, (str, int, float, bool)):
            result[str_key] = value
        elif hasattr(value, 'item'):  # numpy scalar types
            try:
                result[str_key] = value.item()  # Convert numpy scalar to Python type
            except:
                result[str_key] = str(value)
        elif isinstance(value, dict):
            result[str_key] = safe_str_serialization(value)  # Recursive call for nested dicts
        elif isinstance(value, (list, tuple, set)):
            # Handle lists, but make sure all elements are serializable
            result[str_key] = []
            for item in value:
                if isinstance(item, dict):
                    result[str_key].append(safe_str_serialization(item))
                elif hasattr(item, 'item'):
                    try:
                        result[str_key].append(item.item())
                    except:
                        result[str_key].append(str(item))
                elif isinstance(item, (str, int, float, bool)):
                    result[str_key].append(item)
                else:
                    result[str_key].append(str(item))  # Convert other types to string
        else:
            # Fallback for any other type
            result[str_key] = str(value)
    
    return result
"""
    
    # Insert the safe serialization function after the imports
    content = re.sub(r'(import .*?\n\n)', r'\1\n' + safe_serialization_function, content, flags=re.DOTALL)
    
    # Modify the create_standard_anndata function to use safe serialization
    create_standard_anndata_fix = """
    # Ensure all uns values are serializable
    dataset_info = safe_str_serialization(dataset_info)
"""
    
    # Insert the fix in create_standard_anndata before creating the AnnData object
    content = re.sub(r'(\s*# Create AnnData object)', 
                    create_standard_anndata_fix + r'\1', 
                    content)
    
    # Modify the save_anndata function to use safe serialization for all values in uns
    save_anndata_fix = """
        # Convert all values in uns to serializable types
        for key in list(adata.uns.keys()):
            if isinstance(adata.uns[key], dict):
                adata.uns[key] = safe_str_serialization(adata.uns[key])
"""
    
    # Insert the fix in save_anndata before the save operation
    content = re.sub(r'(\s*# Create directory if it doesn\'t exist)', 
                    save_anndata_fix + r'\1', 
                    content)
    
    # Create backup
    backup_file(file_path)
    
    # Write updated content
    with open(file_path, 'w') as f:
        f.write(content)
    
    logger.info("Metadata serialization fix integrated")
    return True

if __name__ == "__main__":
    import time
    
    logger.info("=== Integrating Fixes into RNA-seq Pipeline ===")
    
    # Integrate each fix
    adni_fix = integrate_adni_parsing_fix()
    categorical_fix = integrate_categorical_fix()
    metadata_fix = integrate_metadata_serialization_fix()
    
    # Report results
    logger.info("=== Integration Results ===")
    logger.info(f"ADNI parsing fix: {'SUCCESS' if adni_fix else 'FAILED'}")
    logger.info(f"Categorical data fix: {'SUCCESS' if categorical_fix else 'FAILED'}")
    logger.info(f"Metadata serialization fix: {'SUCCESS' if metadata_fix else 'FAILED'}")
    
    # Update changelog
    if adni_fix or categorical_fix or metadata_fix:
        with open(os.path.join(BASE_DIR, "CHANGELOG.md"), 'a') as f:
            f.write("\n## [0.1.2] - " + time.strftime("%Y-%m-%d") + "\n\n")
            f.write("### Fixed\n")
            if adni_fix:
                f.write("- Fixed parsing of tab-delimited files in ADNI dataset processing\n")
                f.write("- Added handling for escaped tabs in ADNI CSV files\n")
            if categorical_fix:
                f.write("- Fixed categorical data comparison issues in preprocessed datasets\n")
            if metadata_fix:
                f.write("- Fixed metadata serialization issues in standardize_datasets.py\n")
                f.write("- Improved handling of complex nested metadata structures\n")
        
        logger.info("Updated CHANGELOG.md with fixes")
