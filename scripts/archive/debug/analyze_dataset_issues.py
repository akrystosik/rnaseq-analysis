#!/usr/bin/env python3
"""
Script to analyze specific issues in RNA-seq datasets
"""
import os
import sys
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import re

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('dataset_analysis')

def analyze_adni_file(file_path):
    """Analyze an ADNI file for tab delimiter issues."""
    logger.info(f"Analyzing ADNI file: {file_path}")
    
    try:
        # Check for escaped tabs
        with open(file_path, 'r') as f:
            content = f.read()
        
        escaped_tabs = content.count('\\t')
        real_tabs = content.count('\t')
        
        logger.info(f"File content: Escaped tabs: {escaped_tabs}, Real tabs: {real_tabs}")
        
        # Try parsing with tab delimiter
        try:
            df = pd.read_csv(file_path, sep='\t')
            logger.info(f"Successfully parsed with tab delimiter: {len(df)} rows, {len(df.columns)} columns")
            logger.info(f"Columns: {', '.join(df.columns)}")
            return True
        except Exception as e:
            logger.warning(f"Failed to parse with tab delimiter: {e}")
        
        # Try parsing with comma delimiter
        try:
            df = pd.read_csv(file_path)
            logger.info(f"Successfully parsed with comma delimiter: {len(df)} rows, {len(df.columns)} columns")
            logger.info(f"Columns: {', '.join(df.columns)}")
            return True
        except Exception as e:
            logger.warning(f"Failed to parse with comma delimiter: {e}")
        
        logger.error("Failed to parse file with any delimiter")
        return False
    
    except Exception as e:
        logger.error(f"Error analyzing ADNI file: {e}")
        return False

def analyze_categorical_issues(h5ad_file):
    """Analyze categorical column issues in an h5ad file."""
    logger.info(f"Analyzing categorical issues in: {h5ad_file}")
    
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Check for categorical columns in var
        categorical_cols = []
        for col in adata.var.columns:
            if pd.api.types.is_categorical_dtype(adata.var[col]):
                categorical_cols.append(col)
        
        logger.info(f"Found {len(categorical_cols)} categorical columns in var: {', '.join(categorical_cols)}")
        
        # Try a comparison with each categorical column
        for col in categorical_cols:
            test_value = adata.var[col].iloc[0]
            try:
                result = adata.var[col] == test_value
                logger.info(f"Successfully compared {col} with value '{test_value}': {sum(result)} matches")
            except Exception as e:
                logger.error(f"Error comparing {col} with value '{test_value}': {e}")
                
                # Try fixing by converting to string
                logger.info(f"Attempting to fix by converting {col} to string")
                string_col = adata.var[col].astype(str)
                test_value_str = str(test_value)
                try:
                    result = string_col == test_value_str
                    logger.info(f"After conversion to string, comparison works: {sum(result)} matches")
                except Exception as e:
                    logger.error(f"Even after conversion to string, comparison fails: {e}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error analyzing categorical issues: {e}")
        return False

def analyze_metadata_serialization(h5ad_file):
    """Analyze metadata serialization issues in an h5ad file."""
    logger.info(f"Analyzing metadata serialization in: {h5ad_file}")
    
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Check uns keys
        uns_keys = list(adata.uns.keys())
        logger.info(f"Found {len(uns_keys)} keys in uns: {', '.join(uns_keys)}")
        
        # Try to serialize each key
        import json
        
        for key in uns_keys:
            value = adata.uns[key]
            try:
                # Check if value is directly serializable
                json.dumps(value)
                logger.info(f"Key '{key}' is directly serializable")
            except TypeError:
                logger.warning(f"Key '{key}' is not directly serializable")
                
                # Check if it's a dict
                if isinstance(value, dict):
                    logger.info(f"Key '{key}' is a dictionary with {len(value)} entries")
                    
                    # Sample some values
                    sample_values = list(value.items())[:3]
                    for k, v in sample_values:
                        logger.info(f"  Sample: {k} -> {type(v)}")
                elif isinstance(value, (list, tuple)):
                    logger.info(f"Key '{key}' is a {type(value).__name__} with {len(value)} elements")
                    
                    # Sample some values
                    sample_values = value[:3]
                    for i, v in enumerate(sample_values):
                        logger.info(f"  Sample[{i}]: {type(v)}")
                else:
                    logger.info(f"Key '{key}' is a {type(value).__name__}")
        
        # Check if dataset can be saved and loaded
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.h5ad') as tmp:
            try:
                logger.info(f"Attempting to save dataset to temporary file: {tmp.name}")
                adata.write_h5ad(tmp.name)
                logger.info("Successfully saved dataset")
                
                try:
                    logger.info("Attempting to reload saved dataset")
                    reloaded = sc.read_h5ad(tmp.name)
                    logger.info(f"Successfully reloaded dataset: {reloaded.n_obs} samples, {reloaded.n_vars} genes")
                    logger.info(f"Reloaded uns keys: {', '.join(reloaded.uns.keys())}")
                    return True
                except Exception as e:
                    logger.error(f"Error reloading dataset: {e}")
            except Exception as e:
                logger.error(f"Error saving dataset: {e}")
                
                # Try to identify which key is causing the issue
                for key in uns_keys:
                    test_adata = adata.copy()
                    del test_adata.uns[key]
                    try:
                        logger.info(f"Attempting to save dataset without key '{key}'")
                        test_adata.write_h5ad(tmp.name)
                        logger.info(f"Successfully saved dataset without key '{key}' - this key may be the problem")
                    except Exception as e:
                        logger.info(f"Error saving dataset without key '{key}': {e}")
        
        return False
    
    except Exception as e:
        logger.error(f"Error analyzing metadata serialization: {e}")
        return False

def analyze_placeholder_ids(h5ad_file):
    """Analyze placeholder ID issues in an h5ad file."""
    logger.info(f"Analyzing placeholder IDs in: {h5ad_file}")
    
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        
        # Check for gene_id column
        if 'gene_id' not in adata.var.columns:
            logger.warning("No gene_id column found in var")
            return False
        
        # Count placeholder IDs
        placeholder_pattern = re.compile(r'^PLACEHOLDER_')
        placeholder_ids = [g for g in adata.var['gene_id'] if placeholder_pattern.match(str(g))]
        
        logger.info(f"Found {len(placeholder_ids)} placeholder IDs out of {adata.n_vars} genes ({len(placeholder_ids)/adata.n_vars*100:.2f}%)")
        
        # Sample some placeholder IDs
        if placeholder_ids:
            logger.info("Sample placeholder IDs:")
            for i, pid in enumerate(placeholder_ids[:5]):
                logger.info(f"  {i+1}: {pid}")
            
            # Check if there's a pattern to these IDs
            numeric_placeholder = sum(1 for p in placeholder_ids if re.match(r'^PLACEHOLDER_\d+$', str(p)))
            logger.info(f"Numeric placeholders: {numeric_placeholder} ({numeric_placeholder/len(placeholder_ids)*100:.2f}% of placeholders)")
        
        return True
    
    except Exception as e:
        logger.error(f"Error analyzing placeholder IDs: {e}")
        return False

def find_dataset_file(dataset_type):
    """Find a dataset file of the specified type."""
    base_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data"
    
    # Look for the specified dataset type
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.startswith(f"{dataset_type}_") and file.endswith(".h5ad"):
                return os.path.join(root, file)
    
    return None

def find_adni_file():
    """Find an ADNI file for testing."""
    adni_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni_microarray"
    
    for subdir in os.listdir(adni_dir):
        if subdir.startswith("0") and "_S_" in subdir:
            subj_dir = os.path.join(adni_dir, subdir)
            if os.path.isdir(subj_dir):
                for file in os.listdir(subj_dir):
                    if file.endswith(".csv"):
                        return os.path.join(subj_dir, file)
    
    return None

if __name__ == "__main__":
    logger.info("=== Analyzing RNA-seq Dataset Issues ===")
    
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze specific issues in RNA-seq datasets")
    parser.add_argument("--issue", choices=["adni", "categorical", "metadata", "placeholder", "all"], required=True,
                       help="The specific issue to analyze")
    parser.add_argument("--file", type=str, help="Path to a specific file to analyze")
    
    args = parser.parse_args()
    
    if args.issue == "adni" or args.issue == "all":
        adni_file = args.file or find_adni_file()
        if adni_file:
            logger.info(f"Found ADNI file: {adni_file}")
            analyze_adni_file(adni_file)
        else:
            logger.error("No ADNI file found")
    
    if args.issue == "categorical" or args.issue == "all":
        encode_file = args.file or find_dataset_file("encode")
        if encode_file:
            logger.info(f"Found ENCODE file: {encode_file}")
            analyze_categorical_issues(encode_file)
        else:
            logger.error("No ENCODE file found")
        
        preprocessed_dir = os.path.join("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data", "preprocessed_data")
        if os.path.exists(preprocessed_dir):
            for file in os.listdir(preprocessed_dir):
                if file.endswith("_preprocessed.h5ad"):
                    logger.info(f"Found preprocessed file: {os.path.join(preprocessed_dir, file)}")
                    analyze_categorical_issues(os.path.join(preprocessed_dir, file))
                    break
    
    if args.issue == "metadata" or args.issue == "all":
        datasets = ["encode", "gtex", "mage"]
        for dataset in datasets:
            dataset_file = args.file or find_dataset_file(dataset)
            if dataset_file:
                logger.info(f"Found {dataset.upper()} file: {dataset_file}")
                analyze_metadata_serialization(dataset_file)
            else:
                logger.error(f"No {dataset.upper()} file found")
    
    if args.issue == "placeholder" or args.issue == "all":
        preprocessed_dir = os.path.join("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data", "preprocessed_data")
        if os.path.exists(preprocessed_dir):
            for file in os.listdir(preprocessed_dir):
                if file.endswith("_preprocessed.h5ad"):
                    logger.info(f"Found preprocessed file: {os.path.join(preprocessed_dir, file)}")
                    analyze_placeholder_ids(os.path.join(preprocessed_dir, file))
                    break
