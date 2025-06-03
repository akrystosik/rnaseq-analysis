#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/generate_tissue_mapping.py
#!/usr/bin/env python3
"""
Generate Complete Tissue to UBERON Mapping File

This script extracts unique tissue names from all datasets and creates
a CSV mapping file with comprehensive UBERON mappings.
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('tissue_mapping')

def load_tissue_mapping():
    """Load tissue to UBERON mapping from JSON file."""
    mapping_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json"
    
    try:
        with open(mapping_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error loading tissue mapping: {e}")
        # Return empty dictionary as fallback
        return {}

# Replace the hardcoded dictionary with the loaded one
TISSUE_TO_UBERON = load_tissue_mapping()

def extract_tissues_from_datasets(data_dir):
    """
    Extract unique tissue names from all datasets.
    
    Args:
        data_dir: Directory containing h5ad files
        
    Returns:
        List of unique tissue names
    """
    data_dir = Path(data_dir)
    all_tissues = set()
    
    # Find all h5ad files
    h5ad_files = list(data_dir.glob("*_standardized*.h5ad"))
    
    logger.info(f"Found {len(h5ad_files)} h5ad files")
    
    for h5ad_file in h5ad_files:
        try:
            logger.info(f"Processing {h5ad_file}")
            
            # Load the AnnData object
            adata = sc.read_h5ad(h5ad_file)
            
            # Extract tissues if available
            if 'tissue' in adata.obs.columns:
                if isinstance(adata.obs['tissue'], pd.Categorical):
                    tissues = set(adata.obs['tissue'].cat.categories)
                else:
                    tissues = set(adata.obs['tissue'].unique())
                    
                # Remove empty or NaN tissues
                tissues = {t for t in tissues if t and pd.notna(t) and str(t).strip() != ''}
                
                logger.info(f"Found {len(tissues)} unique tissues in {h5ad_file.name}")
                all_tissues.update(tissues)
        except Exception as e:
            logger.error(f"Error processing {h5ad_file}: {e}")
    
    # Convert to list and sort
    all_tissues = sorted(list(all_tissues))
    logger.info(f"Found {len(all_tissues)} unique tissues across all datasets")
    
    return all_tissues

def generate_mapping_file(tissues, output_file):
    """
    Generate a mapping file with comprehensive UBERON mappings.
    
    Args:
        tissues: List of tissue names
        output_file: Output CSV file
    """
    # Create DataFrame with tissues and empty mappings
    df = pd.DataFrame({
        'tissue_name': tissues,
        'ontology_id': '',
        'confidence': '',
        'source': ''
    })
    
    # Fill in mappings from our comprehensive dictionary
    for i, tissue in enumerate(tissues):
        if tissue in TISSUE_TO_UBERON:
            df.loc[i, 'ontology_id'] = TISSUE_TO_UBERON[tissue]
            df.loc[i, 'confidence'] = 'high'
            df.loc[i, 'source'] = 'built-in'
    
    # Count how many we have mappings for
    mapped_count = df['ontology_id'].notna().sum()
    logger.info(f"Generated mappings for {mapped_count}/{len(tissues)} tissues ({mapped_count/len(tissues)*100:.1f}%)")
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    logger.info(f"Saved mapping file to {output_file}")
    
    return df

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate tissue to UBERON mapping file")
    parser.add_argument("--data-dir", required=True, help="Directory containing h5ad files")
    parser.add_argument("--output-file", required=True, help="Output CSV file")
    
    args = parser.parse_args()
    
    # Extract tissues
    tissues = extract_tissues_from_datasets(args.data_dir)
    
    # Generate mapping file
    generate_mapping_file(tissues, args.output_file)

if __name__ == "__main__":
    main()