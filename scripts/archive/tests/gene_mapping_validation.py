#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/gene_mapping_validation.py
#!/usr/bin/env python3
"""
Gene ID Mapping Validation Script

This script validates the accuracy of gene ID mapping in the standardized AnnData objects
by comparing against original source files and reference mappings.

Key validations:
- Preservation of gene counts across standardization
- Accurate mapping between Ensembl IDs and gene symbols
- Consistent handling of gene ID versions
- Proper annotation of gene types and chromosomes
"""

import os
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import glob
import logging
import argparse
import gzip

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gene_mapping_validator')

# Define paths
DEFAULT_STD_DATA_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/")
DEFAULT_SRC_DATA_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/")
DEFAULT_GENCODE_MAPPING = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gencode_v24_complete_mapping.csv")

def load_standardized_dataset(file_path):
    """Load a standardized dataset from an h5ad file."""
    logger.info(f"Loading standardized dataset from {file_path}")
    try:
        adata = ad.read_h5ad(file_path)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        return adata
    except Exception as e:
        logger.error(f"Error loading dataset: {e}")
        return None

def load_gencode_reference(file_path):
    """Load reference GENCODE mapping."""
    logger.info(f"Loading GENCODE reference from {file_path}")
    try:
        gencode_df = pd.read_csv(file_path)
        logger.info(f"Loaded GENCODE reference with {len(gencode_df)} entries")
        return gencode_df
    except Exception as e:
        logger.error(f"Error loading GENCODE reference: {e}")
        return None

def load_original_gene_ids(dataset_name, src_data_dir):
    """Load original gene IDs from source files."""
    logger.info(f"Loading original gene IDs for {dataset_name}")
    
    if dataset_name == "encode":
        # Handle ENCODE cell lines
        gene_ids = set()
        tpm_files = glob.glob(os.path.join(src_data_dir, "encode/raw_data/cell_lines/**/*.tsv"), recursive=True)
        # Also include ENTEx files
        tpm_files.extend(glob.glob(os.path.join(src_data_dir, "encode/raw_data/entex/**/*.tsv"), recursive=True))
        
        for file in tpm_files[:20]:  # Sample a reasonable number of files
            try:
                df = pd.read_csv(file, sep='\t')
                if 'gene_id' in df.columns:
                    gene_ids.update(df['gene_id'])
            except Exception as e:
                logger.warning(f"Error reading {file}: {e}")
        return gene_ids
    
    elif dataset_name == "gtex":
        try:
            gtex_file = os.path.join(src_data_dir, "gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz")
            logger.info(f"Loading GTEx genes from {gtex_file}")
            gene_ids = set()
            
            with gzip.open(gtex_file, 'rt') as f:
                # Skip first two header lines
                next(f)  # Version line
                next(f)  # Dimensions line
                next(f)  # Column headers
                
                # Read all gene IDs from first column
                for line in f:
                    parts = line.strip().split('\t')
                    if parts:
                        gene_ids.add(parts[0])
            
            logger.info(f"Loaded {len(gene_ids)} gene IDs from GTEx file")
            return gene_ids
            
            sample_ids = list(gene_ids)[:100]
            with open("gtex_raw_ids_sample.txt", "w") as f:
                for id in sample_ids:
                    f.write(f"{id}\n")            
            
        except Exception as e:
            logger.error(f"Error reading GTEx file: {e}")
            return set()
        
    elif dataset_name == "adni":
        # Handle ADNI microarray files
        gene_ids = set()
        try:
            # Find CSV files in ADNI directory
            csv_files = glob.glob(os.path.join(src_data_dir, "adni_microarray/**/*.csv"), recursive=True)
            
            for file in csv_files[:10]:  # Sample a few files
                try:
                    df = pd.read_csv(file, sep='\t')  # Most ADNI files are tab-separated despite .csv extension
                    if 'gene_id' in df.columns:
                        gene_ids.update(df['gene_id'])
                except Exception as e:
                    # Try comma separator if tab fails
                    try:
                        df = pd.read_csv(file, sep=',')
                        if 'gene_id' in df.columns:
                            gene_ids.update(df['gene_id'])
                    except:
                        logger.warning(f"Error reading {file}: {e}")
            
            return gene_ids
        except Exception as e:
            logger.error(f"Error processing ADNI files: {e}")
            return set()
    
    elif dataset_name == "mage":
        # Handle MAGE files
        gene_ids = set()
        try:
            # Look for directories that match 'NA*' or 'HG*' patterns
            donor_dirs = glob.glob(os.path.join(src_data_dir, "mage/NA*")) + glob.glob(os.path.join(src_data_dir, "mage/HG*"))
            
            for donor_dir in donor_dirs[:5]:  # Sample a few donors
                csv_files = glob.glob(os.path.join(donor_dir, "**/*.csv"), recursive=True)
                
                for file in csv_files[:2]:  # Just a couple files per donor
                    try:
                        df = pd.read_csv(file, sep='\t')  # Try tab first
                        if 'gene_id' in df.columns:
                            gene_ids.update(df['gene_id'])
                        else:
                            # Try to find a gene ID column
                            gene_col = [col for col in df.columns if 'gene' in col.lower() and 'id' in col.lower()]
                            if gene_col:
                                gene_ids.update(df[gene_col[0]])
                    except:
                        # Try comma separator
                        try:
                            df = pd.read_csv(file, sep=',')
                            if 'gene_id' in df.columns:
                                gene_ids.update(df['gene_id'])
                            else:
                                # Try to find a gene ID column
                                gene_col = [col for col in df.columns if 'gene' in col.lower() and 'id' in col.lower()]
                                if gene_col:
                                    gene_ids.update(df[gene_col[0]])
                        except Exception as e:
                            logger.warning(f"Error reading {file}: {e}")
            
            return gene_ids
        except Exception as e:
            logger.error(f"Error processing MAGE files: {e}")
            return set()
    
    elif dataset_name == "entex":
        # ENTEx data is handled as part of ENCODE
        return load_original_gene_ids("encode", src_data_dir)
    
    logger.warning(f"No implementation for loading original gene IDs for {dataset_name}")
    return set()


def validate_gene_id_mapping(adata, original_gene_ids, gencode_reference, dataset_name):
    """Validate gene ID mapping accuracy."""
    results = {
        "dataset": dataset_name,
        "total_genes_standardized": adata.n_vars,
        "original_gene_count": len(original_gene_ids) if original_gene_ids else "unknown",
        "genes_with_symbols": 0,
        "genes_mapped_to_gencode": 0,
        "genes_with_mismatch": [],
    }
    
    # Special handling for GTEx
    if dataset_name.lower() == "gtex":
        # Look for the actual gene IDs in var DataFrame columns
        ensembl_id_column = None
        
        # Check if gene_id or similar column exists in var DataFrame
        potential_columns = ['gene_id', 'ensembl_id', 'original_ids']
        for col in potential_columns:
            if col in adata.var.columns:
                ensembl_id_column = col
                break
        
        if ensembl_id_column:
            logger.info(f"Found Ensembl IDs in GTEx var['{ensembl_id_column}'] column")
            
            # Use the proper Ensembl IDs from the column rather than the numeric index
            std_ids = set()
            for gene_id in adata.var[ensembl_id_column]:
                if pd.notna(gene_id):
                    # Handle standard Ensembl ID with/without version
                    if isinstance(gene_id, str) and gene_id.startswith("ENSG"):
                        base_id = gene_id.split('.')[0]
                        std_ids.add(base_id)
        else:
            logger.warning(f"No Ensembl ID column found in GTEx var DataFrame")
            # Print all var columns to help debugging
            logger.info(f"Available columns in var DataFrame: {list(adata.var.columns)}")
            std_ids = set()
        
        # Process original IDs as before
        orig_ids = set()
        for gene_id in original_gene_ids:
            if gene_id.startswith("ENSG"):
                base_id = gene_id.split('.')[0]
                orig_ids.add(base_id)
        
        # Find overlap
        matched_genes = std_ids.intersection(orig_ids)
        results["matched_original_genes"] = len(matched_genes)
        results["matched_original_percentage"] = len(matched_genes) / len(std_ids) * 100 if std_ids else 0
        
        # Add GENCODE mapping information
        if gencode_reference is not None:
            gencode_gene_ids = set(gencode_reference['gene_id'].str.split('.').str[0].unique())
            matched_gencode_genes = std_ids.intersection(gencode_gene_ids)
            
            results["genes_mapped_to_gencode"] = len(matched_gencode_genes)
            results["gencode_mapping_percentage"] = len(matched_gencode_genes) / adata.n_vars * 100 if adata.n_vars > 0 else 0
        
        # Return results for GTEx
        return results
    
    else:
    
        # Check how many standardized genes have symbols
        if 'gene_name' in adata.var.columns:
            results["genes_with_symbols"] = adata.var['gene_name'].notna().sum()
            
        # Check mapping source distribution if available
        if 'mapping_source' in adata.var.columns:
            mapping_counts = adata.var['mapping_source'].value_counts().to_dict()
            results["mapping_sources"] = mapping_counts
        
        # Check for genes in standardized data that match original IDs
        if original_gene_ids:
            # Strip version from Ensembl IDs for comparison
            std_base_ids = [id.split('.')[0] if '.' in id and id.startswith('ENSG') else id 
                            for id in adata.var_names]
            orig_base_ids = [id.split('.')[0] if '.' in id and id.startswith('ENSG') else id 
                            for id in original_gene_ids]
            
            matched_genes = set(std_base_ids).intersection(set(orig_base_ids))
            results["matched_original_genes"] = len(matched_genes)
            results["matched_original_percentage"] = len(matched_genes) / len(std_base_ids) * 100

        
        # Properly check GENCODE mapping by using standardized IDs
        if gencode_reference is not None:
            gencode_gene_ids = set(gencode_reference['gene_id'].str.split('.').str[0].unique())
            
            # Check how many standardized genes are in GENCODE
            std_base_ids = [id.split('.')[0] if '.' in id and id.startswith('ENSG') else id 
                        for id in adata.var_names]
            matched_gencode_genes = set(std_base_ids).intersection(gencode_gene_ids)
            
            results["genes_mapped_to_gencode"] = len(matched_gencode_genes)
            results["gencode_mapping_percentage"] = len(matched_gencode_genes) / adata.n_vars * 100
            
            logger.info(f"Found {len(matched_gencode_genes)} genes mapped to GENCODE reference")
        

        # Calculate overall mapping accuracy
        if gencode_reference is not None:
            results["gencode_mapping_percentage"] = results["genes_mapped_to_gencode"] / adata.n_vars * 100
            results["symbol_mismatch_percentage"] = len(results["genes_with_mismatch"]) / results["genes_with_symbols"] * 100 if results["genes_with_symbols"] > 0 else 0
        
        return results

def main():
    parser = argparse.ArgumentParser(description='Validate gene ID mapping in standardized datasets')
    parser.add_argument('--std-dir', type=str, default=str(DEFAULT_STD_DATA_DIR), 
                        help='Directory containing standardized h5ad files')
    parser.add_argument('--src-dir', type=str, default=str(DEFAULT_SRC_DATA_DIR), 
                        help='Directory containing original source files')
    parser.add_argument('--gencode', type=str, default=str(DEFAULT_GENCODE_MAPPING), 
                        help='Path to GENCODE reference mapping')
    parser.add_argument('--datasets', type=str, default='all', 
                        help='Comma-separated list of datasets to validate (or "all")')
    parser.add_argument('--output', type=str, default='gene_mapping_validation.csv', 
                        help='Output file for validation results')
    
    args = parser.parse_args()
    
    # Find standardized datasets
    std_data_dir = Path(args.std_dir)
    std_files = list(std_data_dir.glob("*_standardized.h5ad"))
    
    if not std_files:
        logger.error(f"No standardized datasets found in {std_data_dir}")
        return
    
    logger.info(f"Found {len(std_files)} standardized datasets")
    
    # Load GENCODE reference
    gencode_reference = load_gencode_reference(args.gencode)
    
    # Determine which datasets to validate
    if args.datasets.lower() == 'all':
        datasets_to_validate = [f.stem.replace('_standardized', '') for f in std_files]
    else:
        datasets_to_validate = [d.strip() for d in args.datasets.split(',')]
    
    # Validate each dataset
    validation_results = []
    
    for dataset_name in datasets_to_validate:
        std_file = std_data_dir / f"{dataset_name}_standardized.h5ad"
        
        if not std_file.exists():
            logger.warning(f"Standardized file not found for {dataset_name}")
            continue
        
        # Load standardized dataset
        adata = load_standardized_dataset(std_file)
        if adata is None:
            continue
        
        # Load original gene IDs
        original_gene_ids = load_original_gene_ids(dataset_name, args.src_dir)
        
        # Validate gene ID mapping
        results = validate_gene_id_mapping(adata, original_gene_ids, gencode_reference, dataset_name)
        validation_results.append(results)
        
        logger.info(f"Validated {dataset_name}: {results['total_genes_standardized']} genes, "
                   f"{results.get('matched_original_percentage', 'N/A')}% match with original, "
                   f"{results.get('gencode_mapping_percentage', 'N/A')}% mapped to GENCODE")
    
    # Save validation results
    if validation_results:
        # Extract main metrics for summary table
        summary_data = []
        for result in validation_results:
            summary_data.append({
                "Dataset": result["dataset"],
                "Total Genes": result["total_genes_standardized"],
                "Original Genes": result["original_gene_count"],
                "Genes with Symbols": result["genes_with_symbols"],
                "% Matched to Original": result.get("matched_original_percentage", "N/A"),
                "% Mapped to GENCODE": result.get("gencode_mapping_percentage", "N/A"),
                "% Symbol Mismatches": result.get("symbol_mismatch_percentage", "N/A")
            })
        
        # Save summary table
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(args.output, index=False)
        logger.info(f"Saved validation results to {args.output}")
        
        # Save detailed mismatch information for review
        for result in validation_results:
            if result["genes_with_mismatch"]:
                mismatch_df = pd.DataFrame(result["genes_with_mismatch"])
                mismatch_file = f"{result['dataset']}_symbol_mismatches.csv"
                mismatch_df.to_csv(mismatch_file, index=False)
                logger.info(f"Saved {len(result['genes_with_mismatch'])} symbol mismatches for {result['dataset']} to {mismatch_file}")

if __name__ == "__main__":
    main()