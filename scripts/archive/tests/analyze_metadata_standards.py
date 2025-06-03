#!/usr/bin/env python3
"""
Analyze Metadata Standards in RNA-seq Datasets

This script examines the metadata in standardized RNA-seq datasets to find
examples where ontology standards are met or not met.
"""

import os
import sys
import logging
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('metadata_analyzer')

# Define paths
DATA_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data")
METADATA_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json")

# Custom JSON encoder to handle NumPy types
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return super(NumpyEncoder, self).default(obj)

def load_dataset(file_path):
    """Load a dataset from an h5ad file."""
    try:
        logger.info(f"Loading dataset from {file_path}")
        adata = sc.read_h5ad(file_path)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        return adata
    except Exception as e:
        logger.error(f"Error loading dataset: {e}")
        return None

def load_mapping_files():
    """Load the ontology mapping files."""
    mappings = {}
    mapping_files = {
        'tissue_to_uberon': 'tissue_to_uberon.json',
        'assay_to_efo': 'assay_to_efo.json',
        'age_to_hsapdv': 'age_to_hsapdv.json',
        'species_to_taxon': 'species_to_taxon.json',
        'sex_standardization': 'sex_standardization.json'
    }
    
    for key, filename in mapping_files.items():
        file_path = METADATA_DIR / filename
        try:
            if file_path.exists():
                with open(file_path, 'r') as f:
                    mappings[key] = json.load(f)
                logger.info(f"Loaded {key} mapping with {len(mappings[key])} entries")
            else:
                logger.warning(f"Mapping file not found: {file_path}")
                mappings[key] = {}
        except Exception as e:
            logger.error(f"Error loading mapping file {file_path}: {e}")
            mappings[key] = {}
    
    return mappings

def analyze_tissue_ontology(adata, mappings):
    """Analyze the tissue ontology mappings."""
    results = {
        'total_samples': int(adata.n_obs),
        'has_tissue_field': 'tissue' in adata.obs.columns,
        'has_ontology_field': 'tissue_ontology' in adata.obs.columns,
        'examples': {
            'mapped': [],
            'unmapped': []
        }
    }
    
    if not results['has_tissue_field'] or not results['has_ontology_field']:
        return results
    
    # Count mapped and unmapped
    mapped_mask = (adata.obs['tissue_ontology'].notnull() & 
                  (adata.obs['tissue_ontology'] != '') & 
                  (adata.obs['tissue_ontology'] != 'nan'))
    
    results['mapped_count'] = int(mapped_mask.sum())
    results['unmapped_count'] = int(adata.n_obs - results['mapped_count'])
    results['mapped_percentage'] = float((results['mapped_count'] / adata.n_obs) * 100)
    
    # Find examples of mapped tissues
    if results['mapped_count'] > 0:
        mapped_samples = adata.obs[mapped_mask].iloc[:10]  # Get first 10 mapped samples
        for idx, row in mapped_samples.iterrows():
            results['examples']['mapped'].append({
                'sample_id': str(idx),
                'tissue': str(row['tissue']),
                'tissue_ontology': str(row['tissue_ontology'])
            })
    
    # Find examples of unmapped tissues
    if results['unmapped_count'] > 0:
        unmapped_mask = ~mapped_mask
        unmapped_samples = adata.obs[unmapped_mask].iloc[:10]  # Get first 10 unmapped samples
        for idx, row in unmapped_samples.iterrows():
            results['examples']['unmapped'].append({
                'sample_id': str(idx),
                'tissue': str(row['tissue']),
                'tissue_ontology': str(row['tissue_ontology'] if 'tissue_ontology' in row else '')
            })
    
    # Find unique tissues and their mapping status
    unique_tissues = adata.obs['tissue'].unique()
    tissues_status = {}
    
    for tissue in unique_tissues:
        if pd.isna(tissue) or tissue == '':
            continue
            
        tissue_samples = adata.obs[adata.obs['tissue'] == tissue]
        tissue_ontology_values = tissue_samples['tissue_ontology'].unique()
        mapped_values = [str(v) for v in tissue_ontology_values if v and v != '' and not pd.isna(v)]
        
        tissues_status[str(tissue)] = {
            'mapped': len(mapped_values) > 0,
            'ontology_values': mapped_values,
            'sample_count': int(len(tissue_samples))
        }
    
    results['unique_tissues'] = len(tissues_status)
    results['tissues_status'] = tissues_status
    
    return results

def analyze_species_ontology(adata):
    """Analyze the species ontology mappings."""
    results = {
        'total_samples': int(adata.n_obs),
        'has_species_field': 'species' in adata.obs.columns,
        'has_ontology_field': 'species_ontology' in adata.obs.columns,
        'examples': {
            'mapped': [],
            'unmapped': []
        }
    }
    
    if not results['has_species_field'] or not results['has_ontology_field']:
        return results
    
    # Count mapped and unmapped
    mapped_mask = (adata.obs['species_ontology'].notnull() & 
                  (adata.obs['species_ontology'] != '') & 
                  (adata.obs['species_ontology'] != 'nan'))
    
    results['mapped_count'] = int(mapped_mask.sum())
    results['unmapped_count'] = int(adata.n_obs - results['mapped_count'])
    results['mapped_percentage'] = float((results['mapped_count'] / adata.n_obs) * 100)
    
    # Find examples of mapped species
    if results['mapped_count'] > 0:
        mapped_samples = adata.obs[mapped_mask].iloc[:10]  # Get first 10 mapped samples
        for idx, row in mapped_samples.iterrows():
            results['examples']['mapped'].append({
                'sample_id': str(idx),
                'species': str(row['species']),
                'species_ontology': str(row['species_ontology'])
            })
    
    # Find examples of unmapped species
    if results['unmapped_count'] > 0:
        unmapped_mask = ~mapped_mask
        unmapped_samples = adata.obs[unmapped_mask].iloc[:10]  # Get first 10 unmapped samples
        for idx, row in unmapped_samples.iterrows():
            results['examples']['unmapped'].append({
                'sample_id': str(idx),
                'species': str(row['species']),
                'species_ontology': str(row['species_ontology'] if 'species_ontology' in row else '')
            })
    
    # Find unique species and their mapping status
    unique_species = adata.obs['species'].unique()
    species_status = {}
    
    for species in unique_species:
        if pd.isna(species) or species == '':
            continue
            
        species_samples = adata.obs[adata.obs['species'] == species]
        species_ontology_values = species_samples['species_ontology'].unique()
        mapped_values = [str(v) for v in species_ontology_values if v and v != '' and not pd.isna(v)]
        
        species_status[str(species)] = {
            'mapped': len(mapped_values) > 0,
            'ontology_values': mapped_values,
            'sample_count': int(len(species_samples))
        }
    
    results['unique_species'] = len(species_status)
    results['species_status'] = species_status
    
    return results

def analyze_assay_ontology(adata):
    """Analyze the assay ontology mappings."""
    results = {
        'total_samples': int(adata.n_obs),
        'has_data_type_field': 'data_type' in adata.obs.columns,
        'has_ontology_field': 'assay_ontology' in adata.obs.columns,
        'examples': {
            'mapped': [],
            'unmapped': []
        }
    }
    
    if not results['has_data_type_field'] or not results['has_ontology_field']:
        return results
    
    # Count mapped and unmapped
    mapped_mask = (adata.obs['assay_ontology'].notnull() & 
                  (adata.obs['assay_ontology'] != '') & 
                  (adata.obs['assay_ontology'] != 'nan'))
    
    results['mapped_count'] = int(mapped_mask.sum())
    results['unmapped_count'] = int(adata.n_obs - results['mapped_count'])
    results['mapped_percentage'] = float((results['mapped_count'] / adata.n_obs) * 100)
    
    # Find examples of mapped assays
    if results['mapped_count'] > 0:
        mapped_samples = adata.obs[mapped_mask].iloc[:10]  # Get first 10 mapped samples
        for idx, row in mapped_samples.iterrows():
            results['examples']['mapped'].append({
                'sample_id': str(idx),
                'data_type': str(row['data_type']),
                'assay_ontology': str(row['assay_ontology'])
            })
    
    # Find examples of unmapped assays
    if results['unmapped_count'] > 0:
        unmapped_mask = ~mapped_mask
        unmapped_samples = adata.obs[unmapped_mask].iloc[:10]  # Get first 10 unmapped samples
        for idx, row in unmapped_samples.iterrows():
            results['examples']['unmapped'].append({
                'sample_id': str(idx),
                'data_type': str(row['data_type']),
                'assay_ontology': str(row['assay_ontology'] if 'assay_ontology' in row else '')
            })
    
    # Find unique data types and their mapping status
    unique_data_types = adata.obs['data_type'].unique()
    data_type_status = {}
    
    for data_type in unique_data_types:
        if pd.isna(data_type) or data_type == '':
            continue
            
        data_type_samples = adata.obs[adata.obs['data_type'] == data_type]
        assay_ontology_values = data_type_samples['assay_ontology'].unique()
        mapped_values = [str(v) for v in assay_ontology_values if v and v != '' and not pd.isna(v)]
        
        data_type_status[str(data_type)] = {
            'mapped': len(mapped_values) > 0,
            'ontology_values': mapped_values,
            'sample_count': int(len(data_type_samples))
        }
    
    results['unique_data_types'] = len(data_type_status)
    results['data_type_status'] = data_type_status
    
    return results

def analyze_age_ontology(adata):
    """Analyze the developmental stage ontology mappings."""
    results = {
        'total_samples': int(adata.n_obs),
        'has_age_field': 'age' in adata.obs.columns,
        'has_ontology_field': 'developmental_stage_ontology' in adata.obs.columns,
        'examples': {
            'mapped': [],
            'unmapped': []
        }
    }
    
    if not results['has_age_field'] or not results['has_ontology_field']:
        return results
    
    # Count mapped and unmapped
    mapped_mask = (adata.obs['developmental_stage_ontology'].notnull() & 
                  (adata.obs['developmental_stage_ontology'] != '') & 
                  (adata.obs['developmental_stage_ontology'] != 'nan'))
    
    results['mapped_count'] = int(mapped_mask.sum())
    results['unmapped_count'] = int(adata.n_obs - results['mapped_count'])
    results['mapped_percentage'] = float((results['mapped_count'] / adata.n_obs) * 100)
    
    # Find examples of mapped ages
    if results['mapped_count'] > 0:
        mapped_samples = adata.obs[mapped_mask].iloc[:10]  # Get first 10 mapped samples
        for idx, row in mapped_samples.iterrows():
            results['examples']['mapped'].append({
                'sample_id': str(idx),
                'age': str(row['age']),
                'developmental_stage_ontology': str(row['developmental_stage_ontology'])
            })
    
    # Find examples of unmapped ages
    if results['unmapped_count'] > 0:
        unmapped_mask = ~mapped_mask
        unmapped_samples = adata.obs[unmapped_mask].iloc[:10]  # Get first 10 unmapped samples
        for idx, row in unmapped_samples.iterrows():
            results['examples']['unmapped'].append({
                'sample_id': str(idx),
                'age': str(row['age']),
                'developmental_stage_ontology': str(row['developmental_stage_ontology'] if 'developmental_stage_ontology' in row else '')
            })
    
    # Find unique ages and their mapping status
    unique_ages = adata.obs['age'].unique()
    age_status = {}
    
    for age in unique_ages:
        if pd.isna(age) or age == '':
            continue
            
        age_samples = adata.obs[adata.obs['age'] == age]
        age_ontology_values = age_samples['developmental_stage_ontology'].unique()
        mapped_values = [str(v) for v in age_ontology_values if v and v != '' and not pd.isna(v)]
        
        age_status[str(age)] = {
            'mapped': len(mapped_values) > 0,
            'ontology_values': mapped_values,
            'sample_count': int(len(age_samples))
        }
    
    results['unique_ages'] = len(age_status)
    results['age_status'] = age_status
    
    return results

def analyze_sex_standardization(adata):
    """Analyze the sex standardization."""
    results = {
        'total_samples': int(adata.n_obs),
        'has_sex_field': 'sex' in adata.obs.columns,
        'examples': {}
    }
    
    if not results['has_sex_field']:
        return results
    
    # Count each sex value
    sex_value_counts = adata.obs['sex'].value_counts().to_dict()
    # Convert to regular Python types
    results['sex_values'] = {str(k): int(v) for k, v in sex_value_counts.items()}
    
    # Check if values are standardized
    standard_values = ['male', 'female', 'unknown']
    non_standard_values = [str(val) for val in sex_value_counts.keys() 
                          if str(val) not in standard_values and val and not pd.isna(val)]
    
    results['standardized'] = len(non_standard_values) == 0
    results['non_standard_values'] = non_standard_values
    
    # Examples for each sex value
    for sex_value in sex_value_counts.keys():
        if pd.isna(sex_value) or sex_value == '':
            continue
            
        samples = adata.obs[adata.obs['sex'] == sex_value].iloc[:5]
        results['examples'][str(sex_value)] = []
        
        for idx, row in samples.iterrows():
            results['examples'][str(sex_value)].append({
                'sample_id': str(idx),
                'sex': str(sex_value)
            })
    
    return results

def analyze_gencode_version(adata):
    """Analyze the GENCODE version."""
    results = {
        'has_gencode_version': 'gencode_version' in adata.uns,
        'has_harmonized_gencode_version': 'harmonized_gencode_version' in adata.uns
    }
    
    if results['has_gencode_version']:
        results['gencode_version'] = str(adata.uns['gencode_version'])
    
    if results['has_harmonized_gencode_version']:
        results['harmonized_gencode_version'] = str(adata.uns['harmonized_gencode_version'])
    
    # Check if GENCODE v24 is used
    is_v24 = False
    if results['has_harmonized_gencode_version']:
        harmonized_version = str(adata.uns['harmonized_gencode_version']).replace('v', '')
        is_v24 = harmonized_version == '24'
    
    results['is_gencode_v24'] = is_v24
    
    return results

def analyze_reference_genome(adata):
    """Analyze the reference genome."""
    results = {
        'has_reference_genome': 'reference_genome' in adata.uns,
        'has_harmonized_reference_genome': 'harmonized_reference_genome' in adata.uns
    }
    
    if results['has_reference_genome']:
        results['reference_genome'] = str(adata.uns['reference_genome'])
    
    if results['has_harmonized_reference_genome']:
        results['harmonized_reference_genome'] = str(adata.uns['harmonized_reference_genome'])
    
    # Check if hg38 is used
    is_hg38 = False
    if results['has_harmonized_reference_genome']:
        harmonized_genome = str(adata.uns['harmonized_reference_genome'])
        is_hg38 = harmonized_genome in ['hg38', 'GRCh38']
    
    results['is_hg38'] = is_hg38
    
    return results

# Add a new function to analyze cell type ontology
def analyze_cell_type_ontology(adata):
    """Analyze the cell type ontology mappings."""
    results = {
        'total_samples': int(adata.n_obs),
        'has_cell_type_field': 'cell_type' in adata.obs.columns,
        'has_ontology_field': 'cell_type_ontology' in adata.obs.columns,
        'examples': {
            'mapped': [],
            'unmapped': []
        }
    }
    
    if not results['has_cell_type_field'] or not results['has_ontology_field']:
        return results
    
    # Count mapped and unmapped
    mapped_mask = (adata.obs['cell_type_ontology'].notnull() & 
                  (adata.obs['cell_type_ontology'] != '') & 
                  (adata.obs['cell_type_ontology'] != 'nan'))
    
    results['mapped_count'] = int(mapped_mask.sum())
    results['unmapped_count'] = int(adata.n_obs - results['mapped_count'])
    results['mapped_percentage'] = float((results['mapped_count'] / adata.n_obs) * 100)
    
    # Add example samples
    # [Code to add examples similar to other ontology functions]
    
    return results


def analyze_dataset(dataset_name, file_path, mappings):
    """Analyze a dataset for metadata standards."""
    adata = load_dataset(file_path)
    if adata is None:
        return None
    
    results = {
        'dataset': dataset_name,
        'samples': int(adata.n_obs),
        'genes': int(adata.n_vars),
        'tissue_ontology': analyze_tissue_ontology(adata, mappings),
        'species_ontology': analyze_species_ontology(adata),
        'assay_ontology': analyze_assay_ontology(adata),
        'age_ontology': analyze_age_ontology(adata),
        'sex_standardization': analyze_sex_standardization(adata),
        'gencode_version': analyze_gencode_version(adata),
        'reference_genome': analyze_reference_genome(adata), 
        'cell_type_ontology': analyze_cell_type_ontology(adata), 
        
    }
    
    return results


def main():
    # Load mapping files
    mappings = load_mapping_files()
    
    # Define datasets to analyze
    dataset_files = {
        'adni': DATA_DIR / 'adni_standardized_v2.h5ad',
        'encode': DATA_DIR / 'encode_standardized_v2.h5ad',
        'entex': DATA_DIR / 'entex_standardized_v2.h5ad',
        'mage': DATA_DIR / 'mage_standardized_v2.h5ad',
        'gtex': DATA_DIR / 'gtex_standardized_v2.h5ad',
        'combined_all_genes': DATA_DIR / 'combined_all_genes_standardized.h5ad'
    }
    
    # Analyze each dataset
    results = {}
    
    for name, file_path in dataset_files.items():
        if file_path.exists():
            logger.info(f"Analyzing {name} dataset")
            dataset_results = analyze_dataset(name, file_path, mappings)
            results[name] = dataset_results
            logger.info(f"Completed analysis of {name}")
        else:
            logger.warning(f"Dataset file not found: {file_path}")
    
    # Save results to JSON file
    output_file = DATA_DIR / 'metadata_analysis_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)  # Use custom encoder
    
    logger.info(f"Analysis results saved to {output_file}")
    
    # Print summary
    logger.info("=== Metadata Standards Summary ===")
    for name, result in results.items():
        if result is None:
            continue
            
        logger.info(f"\nDataset: {name}")
        
        # GENCODE version
        gencode_info = result['gencode_version']
        if gencode_info['is_gencode_v24']:
            logger.info(f"  GENCODE v24: ✓ (Harmonized version: {gencode_info.get('harmonized_gencode_version', 'N/A')})")
        else:
            logger.info(f"  GENCODE v24: ✗ (Harmonized version: {gencode_info.get('harmonized_gencode_version', 'N/A')})")
        
        # Reference genome
        genome_info = result['reference_genome']
        if genome_info['is_hg38']:
            logger.info(f"  hg38: ✓ (Harmonized genome: {genome_info.get('harmonized_reference_genome', 'N/A')})")
        else:
            logger.info(f"  hg38: ✗ (Harmonized genome: {genome_info.get('harmonized_reference_genome', 'N/A')})")
        
        # Tissue ontology
        tissue_info = result['tissue_ontology']
        if tissue_info['has_tissue_field'] and tissue_info['has_ontology_field']:
            mapped_pct = tissue_info.get('mapped_percentage', 0)
            logger.info(f"  Tissue (UBERON): {mapped_pct:.1f}% mapped")
        else:
            logger.info("  Tissue (UBERON): Missing fields")


        # Cell Type  ontology
        logger.info(f"  Cell Type (CL): {cell_type_info.get('mapped_percentage', 0):.1f}% mapped")
 
        # Species ontology
        species_info = result['species_ontology']
        if species_info['has_species_field'] and species_info['has_ontology_field']:
            mapped_pct = species_info.get('mapped_percentage', 0)
            logger.info(f"  Species (NCBI Taxon): {mapped_pct:.1f}% mapped")
        else:
            logger.info("  Species (NCBI Taxon): Missing fields")
        
        # Assay ontology
        assay_info = result['assay_ontology']
        if assay_info['has_data_type_field'] and assay_info['has_ontology_field']:
            mapped_pct = assay_info.get('mapped_percentage', 0)
            logger.info(f"  Assay (EFO): {mapped_pct:.1f}% mapped")
        else:
            logger.info("  Assay (EFO): Missing fields")
        
        # Age ontology
        age_info = result['age_ontology']
        if age_info['has_age_field'] and age_info['has_ontology_field']:
            mapped_pct = age_info.get('mapped_percentage', 0)
            logger.info(f"  Age (HsapDv): {mapped_pct:.1f}% mapped")
        else:
            logger.info("  Age (HsapDv): Missing fields")
        
        # Sex standardization
        sex_info = result['sex_standardization']
        if sex_info['has_sex_field']:
            if sex_info['standardized']:
                logger.info("  Sex standardization: ✓")
            else:
                non_std = ', '.join(sex_info['non_standard_values'])
                logger.info(f"  Sex standardization: ✗ (Non-standard values: {non_std})")
        else:
            logger.info("  Sex standardization: Missing field")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())