#!/usr/bin/env python3
"""
RNA-seq Metadata Standardization Script

This script standardizes metadata for RNA-seq datasets, adding consistent
annotations, ontology mappings, and inferred metadata where possible.
It now uses JSON configuration files for ontology mappings and dataset-specific metadata.
"""

import os
import re
import logging
import argparse
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import numpy as np

# Import shared utilities
from rnaseq_utils import (
    standardize_metadata, load_dataset_specific_metadata, apply_dataset_specific_metadata,
    load_mappings
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('metadata_standardization')

# Define paths and constants
BASE_DIR = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
DEFAULT_METADATA_DIR = BASE_DIR / "metadata/json"
DEFAULT_OUTPUT_DIR = BASE_DIR / "standardized_data"


class TissueOntologyMapper:
    """
    A class to handle mapping of tissue terms to standardized ontology terms.
    """    
    def _load_json_mapping(self, mapping_file):
        """Load mappings from a JSON file"""
        try:
            import json
            with open(mapping_file, 'r') as f:
                mapping_data = json.load(f)
                for tissue, ontology_id in mapping_data.items():
                    # Handle either simple format or detailed format
                    if isinstance(ontology_id, dict) and 'id' in ontology_id:
                        self.tissue_mappings[tissue] = ontology_id['id']
                        self.mapping_confidence[tissue] = ontology_id.get('confidence', 'medium')
                    else:
                        self.tissue_mappings[tissue] = ontology_id
                        self.mapping_confidence[tissue] = 'medium'
            logger.info(f"Loaded {len(self.tissue_mappings)} tissue mappings from {mapping_file}")
        except Exception as e:
            logger.error(f"Error loading tissue mappings from JSON: {e}")
            
    def __init__(self, mapping_file=None, cache_file=None):
        """
        Initialize the tissue mapper with optional mapping files.
        
        Args:
            mapping_file: Path to mapping file (JSON or CSV)
            cache_file: Path to save/load cached mapping results
        """
        self.tissue_mappings = {}
        self.mapping_confidence = {}
        self.cache_file = cache_file
        
        # First try to load the default JSON mapping file
        default_mapping_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json"
        if os.path.exists(default_mapping_file):
            self._load_json_mapping(default_mapping_file)
            logger.info(f"Loaded default tissue mapping from {default_mapping_file}")
        
        # Then load custom mapping if provided (allows overrides)
        if mapping_file and os.path.exists(mapping_file):
            if mapping_file.endswith('.json'):
                self._load_json_mapping(mapping_file)
            else:
                self._load_csv_mapping(mapping_file)
            
        # Finally load cached mappings if available
        if cache_file and os.path.exists(cache_file):
            self._load_cache()
            
    
    def _lookup_uberon_term(self, tissue_name):
        """
        Lookup a tissue name in the loaded mappings.
        Returns (ontology_id, confidence) tuple.
        """
        # Check if we have this in our mappings
        if tissue_name in self.tissue_mappings:
            return self.tissue_mappings[tissue_name], self.mapping_confidence.get(tissue_name, 'medium')
        
        # Try with lowercase
        tissue_lower = tissue_name.lower()
        if tissue_lower in self.tissue_mappings:
            return self.tissue_mappings[tissue_lower], self.mapping_confidence.get(tissue_lower, 'medium')
        
        # Now try with a more flexible approach by checking for substrings
        for known_tissue, ontology_id in self.tissue_mappings.items():
            if known_tissue.lower() in tissue_lower or tissue_lower in known_tissue.lower():
                # Match with lower confidence since it's not an exact match
                return ontology_id, 'low'
        
        # For now, return None if not found
        return None, 'none'
    
    def _load_existing_mappings(self, mapping_file):
        """Load mappings from a CSV file"""
        try:
            df = pd.read_csv(mapping_file)
            for _, row in df.iterrows():
                tissue = row['tissue_name']
                ontology_id = row['ontology_id']
                confidence = row.get('confidence', 'medium')
                
                self.tissue_mappings[tissue] = ontology_id
                self.mapping_confidence[tissue] = confidence
                
            logger.info(f"Loaded {len(self.tissue_mappings)} tissue mappings from {mapping_file}")
        except Exception as e:
            logger.error(f"Error loading tissue mappings: {e}")
    
    def _load_cache(self):
        """Load cached mapping results"""
        try:
            import json
            with open(self.cache_file, 'r') as f:
                cache_data = json.load(f)
                self.tissue_mappings.update(cache_data.get('mappings', {}))
                self.mapping_confidence.update(cache_data.get('confidence', {}))
            logger.info(f"Loaded {len(self.tissue_mappings)} tissue mappings from cache")
        except Exception as e:
            logger.error(f"Error loading mapping cache: {e}")
    
    def _save_cache(self):
        """Save current mappings to cache file"""
        if not self.cache_file:
            return
            
        try:
            import json
            cache_data = {
                'mappings': self.tissue_mappings,
                'confidence': self.mapping_confidence
            }
            with open(self.cache_file, 'w') as f:
                json.dump(cache_data, f, indent=2)
            logger.info(f"Saved {len(self.tissue_mappings)} tissue mappings to cache")
        except Exception as e:
            logger.error(f"Error saving mapping cache: {e}")
    
    def map_tissue(self, tissue_name):
        """
        Map a tissue name to an ontology ID.
        Returns (ontology_id, confidence) tuple.
        """
        # Handle edge cases
        if pd.isna(tissue_name) or not tissue_name or tissue_name.strip() == '':
            return '', 'none'
        
        # Normalize tissue name
        tissue_name = str(tissue_name).strip()
        
        # Return cached result if available
        if tissue_name in self.tissue_mappings:
            return (self.tissue_mappings[tissue_name], 
                    self.mapping_confidence.get(tissue_name, 'medium'))
                    
        # Try with lowercase version
        tissue_lower = tissue_name.lower()
        if tissue_lower in self.tissue_mappings:
            return (self.tissue_mappings[tissue_lower],
                    self.mapping_confidence.get(tissue_lower, 'medium'))
        
        # Lookup in UBERON
        ontology_id, confidence = self._lookup_uberon_term(tissue_name)
        
        # Cache the result
        if ontology_id:
            self.tissue_mappings[tissue_name] = ontology_id
            self.mapping_confidence[tissue_name] = confidence
            if self.cache_file:
                self._save_cache()
        
        return ontology_id or '', confidence
    
    def map_tissues_in_adata(self, adata, tissue_field='tissue'):
        """Add ontology mappings to an AnnData object."""
        if tissue_field not in adata.obs.columns:
            logger.warning(f"Tissue field '{tissue_field}' not found in adata.obs")
            return adata
        
        # Initialize tissue ontology columns - use string dtype to avoid implicit conversion issues
        adata.obs['tissue_ontology'] = pd.Series(index=adata.obs.index, dtype='str')
        adata.obs['tissue_ontology_confidence'] = pd.Series(index=adata.obs.index, dtype='str')
        
        # Track mapping statistics
        total_samples = adata.n_obs
        mapped_samples = 0
        confidence_counts = {'high': 0, 'medium': 0, 'low': 0, 'none': 0}
        
        try:
            # Filter out empty or NaN tissue values
            if isinstance(adata.obs[tissue_field], pd.Categorical):
                # Get unique non-empty categories
                unique_tissues = [t for t in adata.obs[tissue_field].cat.categories if t and pd.notna(t) and t != '']
            else:
                unique_tissues = adata.obs[tissue_field].dropna().unique()
                unique_tissues = [t for t in unique_tissues if t and str(t).strip() != '']
            
            logger.info(f"Mapping {len(unique_tissues)} unique tissues to ontology terms")
            
            for tissue in unique_tissues:
                try:
                    ontology_id, confidence = self.map_tissue(str(tissue))
                    
                    # Create mask for this tissue, ensuring tissue is treated as a string
                    mask = adata.obs[tissue_field].astype(str) == str(tissue)
                    
                    if ontology_id:
                        # Avoid implicit conversion by converting everything to strings
                        adata.obs.loc[mask, 'tissue_ontology'] = str(ontology_id)
                        adata.obs.loc[mask, 'tissue_ontology_confidence'] = str(confidence)
                        mapped_samples += mask.sum()
                        confidence_counts[confidence] += mask.sum()
                    else:
                        # Mark as unmapped but don't fail
                        adata.obs.loc[mask, 'tissue_ontology'] = ''
                        adata.obs.loc[mask, 'tissue_ontology_confidence'] = 'none'
                        confidence_counts['none'] += mask.sum()
                except Exception as e:
                    logger.warning(f"Error mapping tissue '{tissue}': {e}")
                    # Continue with next tissue
            
            # Handle empty tissue values
            empty_mask = (adata.obs[tissue_field].isna()) | (adata.obs[tissue_field].astype(str) == '')
            if empty_mask.any():
                logger.warning(f"Found {empty_mask.sum()} samples with empty tissue values")
                adata.obs.loc[empty_mask, 'tissue_ontology'] = ''
                adata.obs.loc[empty_mask, 'tissue_ontology_confidence'] = 'none'
                confidence_counts['none'] += empty_mask.sum()
        
        except Exception as e:
            logger.error(f"Error during tissue mapping: {e}")
            # Don't raise the exception, just log it and continue
        
        # Log mapping statistics
        mapping_percentage = (mapped_samples / total_samples) * 100 if total_samples > 0 else 0
        logger.info(f"Tissue ontology mapping complete: {mapped_samples}/{total_samples} ({mapping_percentage:.1f}%)")
        logger.info(f"Mapping confidence: high={confidence_counts['high']}, "
                f"medium={confidence_counts['medium']}, "
                f"low={confidence_counts['low']}, "
                f"none={confidence_counts['none']}")
        
        return adata
        
def infer_rnaseq_protocol(adata):
    """
    Infer RNA-seq protocol type from metadata.
    
    Args:
        adata: AnnData object
        
    Returns:
        Dictionary with protocol information
    """
    protocol_info = {
        'protocol_type': 'RNA-seq',
        'protocol_confidence': 'low',
        'selection_method': 'unknown',
        'sequencing_platform': 'unknown',
        'read_configuration': 'unknown',
        'strand_specificity': 'unknown'
    }
    
    # Extract metadata fields that might contain protocol information
    metadata_fields = adata.uns.keys()
    obs_fields = adata.obs.columns
    
    # Look for GTEx-specific protocol information
    if 'dataset' in adata.obs.columns and len(adata.obs) > 0:
        if 'gtex' in str(adata.obs['dataset'].iloc[0]).lower():
            protocol_info['protocol_type'] = 'RNA-seq polyA+'
            protocol_info['protocol_confidence'] = 'high'
            protocol_info['selection_method'] = 'polyA+ selection using Illumina TruSeq protocol'
            protocol_info['sequencing_platform'] = 'Illumina HiSeq 2000/2500'
            protocol_info['read_configuration'] = '76-bp paired-end reads'
            protocol_info['strand_specificity'] = 'Non-strand specific'
            return protocol_info
    
    # Check for protocol keywords in metadata fields
    protocol_keywords = {
        'polyA': {'protocol_type': 'RNA-seq polyA+', 'selection_method': 'polyA selection'},
        'ribo': {'protocol_type': 'RNA-seq ribo-depleted', 'selection_method': 'ribosomal RNA depletion'},
        'truseq': {'protocol_type': 'RNA-seq', 'selection_method': 'Illumina TruSeq protocol'},
        'total RNA': {'protocol_type': 'RNA-seq total', 'selection_method': 'total RNA'},
        'strand': {'strand_specificity': 'Strand-specific'},
        'unstranded': {'strand_specificity': 'Non-strand specific'},
        'illumina': {'sequencing_platform': 'Illumina'},
        'paired': {'read_configuration': 'paired-end'},
        'single': {'read_configuration': 'single-end'}
    }
    
    # Search all potential metadata locations
    for field in metadata_fields:
        field_value = str(adata.uns.get(field, '')).lower()
        for keyword, info in protocol_keywords.items():
            if keyword.lower() in field_value:
                for key, value in info.items():
                    protocol_info[key] = value
                protocol_info['protocol_confidence'] = 'medium'
    
    # Also check obs fields for protocol info
    for field in obs_fields:
        if len(adata.obs) > 0:
            sample_value = str(adata.obs[field].iloc[0]).lower()
            for keyword, info in protocol_keywords.items():
                if keyword.lower() in sample_value:
                    for key, value in info.items():
                        protocol_info[key] = value
                    protocol_info['protocol_confidence'] = 'medium'
    
    return protocol_info

def infer_gene_ids_format(adata):
    """
    Infer GENCODE version from gene ID format.
    
    Args:
        adata: AnnData object
        
    Returns:
        Dictionary with gene ID information
    """
    gene_info = {
        'id_format': 'unknown',
        'gencode_version': None,
        'confidence': 'low',
        'notes': ''
    }
    
    # Get a sample of gene IDs (first 5)
    if adata.n_vars > 0:
        sample_ids = adata.var_names[:5].tolist()
        logger.info(f"Sample gene IDs: {', '.join(map(str, sample_ids))}")
        
        # Check for Ensembl gene ID format (ENSG)
        if all(re.match(r'ENSG\d+', str(gene_id)) for gene_id in sample_ids):
            gene_info['id_format'] = 'ensembl'
            gene_info['confidence'] = 'high'
            gene_info['notes'] = 'Gene IDs follow ENSEMBL format (ENSG)'
            
            # Try to infer GENCODE version based on metadata
            if 'gencode_version' in adata.uns:
                gene_info['gencode_version'] = str(adata.uns['gencode_version'])
                gene_info['confidence'] = 'high'
                gene_info['notes'] += f", explicit GENCODE v{adata.uns['gencode_version']} annotation"
    
    return gene_info

def infer_species(adata):
    """
    Infer species from metadata or gene IDs.
    
    Args:
        adata: AnnData object
        
    Returns:
        Dictionary with species information
    """
    species_info = {
        'species': 'unknown',
        'ontology': '',
        'confidence': 'low'
    }
    
    # Check for existing species information
    if 'species' in adata.obs.columns and len(adata.obs) > 0:
        first_species = str(adata.obs['species'].iloc[0])
        if 'human' in first_species.lower() or 'homo sapiens' in first_species.lower():
            species_info['species'] = 'Homo sapiens'
            species_info['ontology'] = 'NCBITaxon:9606'
            species_info['confidence'] = 'high'
            return species_info
    
    # Infer from gene IDs
    gene_ids = adata.var_names[:10].tolist()
    if all(re.match(r'ENSG\d+', str(gene_id)) for gene_id in gene_ids):
        # Human ENSEMBL gene IDs start with ENSG
        species_info['species'] = 'Homo sapiens'
        species_info['ontology'] = 'NCBITaxon:9606'
        species_info['confidence'] = 'high'
    
    return species_info

def infer_genome_version(adata):
    """
    Infer genome version from metadata.
    
    Args:
        adata: AnnData object
        
    Returns:
        Dictionary with genome version information
    """
    genome_info = {
        'genome_version': 'unknown',
        'confidence': 'low',
        'notes': ''
    }
    
    # Check if chromosome naming convention exists
    chr_naming = None
    if 'var' in dir(adata) and adata.n_vars > 0:
        if 'chromosome' in adata.var.columns:
            sample_chr = adata.var['chromosome'][:5].tolist()
            if any(str(c).startswith('chr') for c in sample_chr):
                chr_naming = 'chr-prefixed'
                genome_info['confidence'] = 'medium'
                genome_info['notes'] = 'Chromosome names are prefixed with "chr"'
            else:
                chr_naming = 'ensembl-style'
                genome_info['confidence'] = 'medium'
                genome_info['notes'] = 'Chromosome names use Ensembl style (without "chr" prefix)'
    
    # Map chromosome naming to likely genome version
    if chr_naming == 'chr-prefixed':
        genome_info['genome_version'] = 'hg38'
        genome_info['confidence'] = 'medium'
        genome_info['notes'] += ', which is common in hg38'
    elif chr_naming == 'ensembl-style':
        genome_info['genome_version'] = 'GRCh37'
        genome_info['confidence'] = 'medium'
        genome_info['notes'] += ', which is common in GRCh37'
    
    # Check for explicit genome version in metadata
    if 'reference_genome' in adata.uns:
        genome_info['genome_version'] = str(adata.uns['reference_genome'])
        genome_info['confidence'] = 'high'
        genome_info['notes'] = 'Genome version explicitly stated in metadata'
    
    return genome_info

def enhance_standardized_metadata(dataset_name, adata, output_file=None, mapping_dir=None, metadata_dir=None):
    """
    Enhance metadata standardization beyond the basic standardization,
    adding inferred fields and detailed ontology mappings.
    
    Args:
        dataset_name: Name of the dataset
        adata: AnnData object with dataset
        output_file: Optional file to save standardized data
        mapping_dir: Directory containing mapping files
        metadata_dir: Directory containing dataset-specific metadata JSON files
        
    Returns:
        adata with enhanced standardized metadata
    """
    logger.info(f"Enhancing metadata for {dataset_name} dataset")
    
    # Initialize container for standardized metadata
    if 'dataset' not in adata.uns:
        adata.uns['dataset'] = dataset_name
    
    # Initialize gene_info if not present
    if 'gene_info' not in adata.uns:
        adata.uns['gene_info'] = {
            'id_format': 'unknown',
            'gencode_version': None,
            'confidence': 'low',
            'notes': 'Added during standardization'
        }
            
    # Check gene ID format to infer GENCODE version
    logger.info("Checking gene ID format to infer GENCODE version")
    gene_info = infer_gene_ids_format(adata)
    
    # Add dataset name to each cell if not present
    if 'dataset' not in adata.obs.columns:
        adata.obs['dataset'] = dataset_name
    
    # Add sample IDs if not present
    if 'sample_id' not in adata.obs.columns and adata.obs.index.name != 'sample_id':
        adata.obs['sample_id'] = adata.obs.index
    
    # Initialize tissue ontology mapper
    tissue_mapper = None
    if mapping_dir:
        uberon_mapping_file = os.path.join(mapping_dir, 'tissue_to_uberon.csv')
        cache_file = os.path.join(mapping_dir, 'tissue_mapping_cache.json')
        tissue_mapper = TissueOntologyMapper(
            uberon_mapping_file=uberon_mapping_file,
            cache_file=cache_file
        )
    else:
        # Use default mappings without files
        tissue_mapper = TissueOntologyMapper()
    
    # Apply dataset-specific metadata from JSON file if available
    # We do this early so that our metadata can influence the standardization process
    if metadata_dir:
        dataset_specific_metadata = load_dataset_specific_metadata(metadata_dir, dataset_name)
        if dataset_specific_metadata:
            adata = apply_dataset_specific_metadata(adata, dataset_specific_metadata)
 
 
    # Apply cell type information if available
    if 'cell_type_info' in dataset_specific_metadata:
        cell_info = dataset_specific_metadata['cell_type_info']
        
        # Add cell type columns if not present
        if 'cell_type' not in adata.obs.columns:
            adata.obs['cell_type'] = ''
        
        if 'cell_type_ontology' not in adata.obs.columns:
            adata.obs['cell_type_ontology'] = ''
        
        if 'anatomical_entity' not in adata.obs.columns:
            adata.obs['anatomical_entity'] = ''
        
        if 'anatomical_entity_ontology' not in adata.obs.columns:
            adata.obs['anatomical_entity_ontology'] = ''
        
        # Apply to all samples in the dataset
        if 'cell_type_id' in cell_info:
            adata.obs['cell_type_ontology'] = cell_info['cell_type_id']
        
        if 'anatomical_entity_id' in cell_info:
            adata.obs['anatomical_entity_ontology'] = cell_info['anatomical_entity_id']
        
        # Store in uns metadata
        adata.uns['cell_type_info'] = cell_info
            
    
    # Infer RNA-seq protocol
    logger.info(f"Inferring RNA-seq protocol for {dataset_name}")
    rnaseq_protocol = infer_rnaseq_protocol(adata)
    adata.uns['rna_seq_protocol'] = rnaseq_protocol
    if rnaseq_protocol.get('protocol_confidence', 'low') == 'low':
        logger.warning(f"Could not determine RNA-seq protocol for {dataset_name} from metadata.")
    
    # Infer species
    logger.info(f"Inferring species for {dataset_name}")
    species_info = infer_species(adata)
    
    # Add species information as string to avoid conversion issues
    if 'species' not in adata.obs.columns:
        adata.obs['species'] = species_info['species']
    
    if 'species_ontology' not in adata.obs.columns:
        adata.obs['species_ontology'] = species_info['ontology']
    
    logger.info(f"Dataset {dataset_name} gene IDs suggest {species_info['species']} species")
    
    # Infer genome version
    logger.info(f"Inferring genome version for {dataset_name}")
    genome_info = infer_genome_version(adata)
    if genome_info['genome_version'] != 'unknown':
        logger.info(f"Dataset {dataset_name} chromosome names suggest {genome_info['genome_version']} ({genome_info.get('notes', '')})")
        
        # Infer GENCODE version from genome version if not already known
        if not gene_info['gencode_version'] and genome_info['genome_version'] == 'hg38':
            logger.info(f"Dataset {dataset_name} uses hg38, which commonly uses GENCODE v24")
            gene_info['gencode_version'] = '24'
            gene_info['confidence'] = 'medium'
            gene_info['notes'] += ", inferred from genome version"
            
        if not gene_info['gencode_version'] and genome_info['genome_version'] == 'GRCh37':
            logger.info(f"Dataset {dataset_name} uses GRCh37, which commonly uses GENCODE v19")
            gene_info['gencode_version'] = '19'
            gene_info['confidence'] = 'medium'
            gene_info['notes'] += ", inferred from genome version"
    
    # Store inferred gene and genome information as strings
    adata.uns['gene_info'] = {k: str(v) if v is not None else None for k, v in gene_info.items()}
    adata.uns['genome_info'] = {k: str(v) if v is not None else None for k, v in genome_info.items()}
    
    # Add tissue ontology mapping
    if tissue_mapper and 'tissue' in adata.obs.columns:
        logger.info(f"Adding tissue ontology mappings for {dataset_name}")
        adata = tissue_mapper.map_tissues_in_adata(adata, tissue_field='tissue')
    
    # Add assay ontology if not present - ensure string dtype
    if 'assay_ontology' not in adata.obs.columns:
        adata.obs['assay_ontology'] = 'EFO:0009922'  # RNA-seq
    
    # Add developmental stage ontology if not present - ensure string dtype
    if 'developmental_stage_ontology' not in adata.obs.columns:
        # Set to human adult by default
        adata.obs['developmental_stage_ontology'] = 'HsapDv:0000087'
    
    # Handle GENCODE version harmonization
    if gene_info['gencode_version'] or 'gencode_version' in adata.uns:
        # Get original GENCODE version
        if 'gencode_version' in adata.uns:
            original_gencode = str(adata.uns['gencode_version']).replace('v', '')
        else:
            original_gencode = str(gene_info['gencode_version']).replace('v', '')
            
        adata.uns['original_gencode_version'] = original_gencode
        
        # Set harmonized GENCODE version
        if 'harmonized_gencode_version' in adata.uns:
            harmonized_gencode = str(adata.uns['harmonized_gencode_version']).replace('v', '')
        else:
            # Default harmonized version is v24
            harmonized_gencode = '24'
            
        adata.uns['harmonized_gencode_version'] = harmonized_gencode
        
        # Add mapping notes if versions differ
        if original_gencode != harmonized_gencode:
            adata.uns['gencode_mapping_notes'] = (
                f"Original {dataset_name} data was aligned using GENCODE "
                f"v{original_gencode} annotation. For harmonization "
                f"across datasets, gene IDs were mapped to GENCODE v{harmonized_gencode} "
                f"to ensure consistent comparison with other RNA-seq datasets."
            )
            
    
    # Handle genome version harmonization
    if 'reference_genome' in adata.uns:
        # Original version from metadata
        original_genome = str(adata.uns['reference_genome'])
    else:
        # Inferred version
        original_genome = str(genome_info.get('genome_version', 'unknown'))

    adata.uns['original_reference_genome'] = original_genome

    # Set harmonized genome version
    if 'harmonized_reference_genome' in adata.uns:
        harmonized_genome = str(adata.uns['harmonized_reference_genome'])
    else:
        # Default harmonized version
        harmonized_genome = 'hg38'

    adata.uns['harmonized_reference_genome'] = harmonized_genome

    # Add mapping notes if versions differ
    if original_genome != harmonized_genome and original_genome != 'unknown':
        adata.uns['genome_mapping_notes'] = (
            f"Original {dataset_name} data was aligned to {original_genome}. "
            f"For harmonization across datasets, coordinates were lifted over to {harmonized_genome} "
            f"to ensure consistent comparison with other datasets."
        )
    
    # Log standardization completion
    logger.info(f"Metadata standardization complete for {dataset_name}")
    logger.info(f"  Samples: {adata.n_obs}")
    logger.info(f"  Genes: {adata.n_vars}")
    logger.info(f"  Reference Genome: {original_genome} (harmonized to {harmonized_genome})")
    logger.info(f"  GENCODE Version: {original_gencode if 'original_gencode' in locals() else 'unknown'} (harmonized to v{harmonized_gencode if 'harmonized_gencode' in locals() else 'unknown'})")
    logger.info(f"  RNA-seq Type: {rnaseq_protocol.get('protocol_type', 'unknown')} (confidence: {rnaseq_protocol.get('protocol_confidence', 'low')})")
    
    # Calculate quality control metrics
    logger.info("Quality control checks:")
    qc_fields = {
        'subject_id': 'subject_id',
        'sex': 'sex',
        'age': 'age',
        'tissue': 'tissue',
        'species': 'species',
        'tissue_ontology': 'tissue_ontology',
        'assay_ontology': 'assay_ontology',
        'developmental_stage_ontology': 'developmental_stage_ontology'
    }
    
    for field_name, obs_name in qc_fields.items():
        if obs_name in adata.obs.columns:
            non_na_count = adata.obs[obs_name].notna().sum()
            percentage = (non_na_count / adata.n_obs) * 100
            logger.info(f"  {field_name}: {non_na_count}/{adata.n_obs} ({percentage:.1f}%)")
        else:
            logger.info(f"  {field_name}: 0/{adata.n_obs} (0.0%)")
    
    # Save standardized data if output file is specified
    if output_file:
        logger.info(f"Saving standardized {dataset_name} dataset to {output_file}")
        
        # Ensure categorical columns are properly handled - to avoid conversion errors
        # Convert columns to strings before saving when needed
        string_columns = ['tissue_ontology', 'tissue_ontology_confidence', 
                          'species', 'species_ontology', 'assay_ontology', 
                          'developmental_stage_ontology']
        
        for col in string_columns:
            if col in adata.obs.columns:
                # Make sure the column is actually string type
                adata.obs[col] = adata.obs[col].astype(str)
        
        try:            
            # Convert all uns values to strings for safety
            safe_uns = {}
            for k, v in adata.uns.items():
                if isinstance(v, dict):
                    safe_uns[k] = {sk: str(sv) if sv is not None else None 
                                for sk, sv in v.items()}
                elif isinstance(v, (list, np.ndarray)):
                    safe_uns[k] = [str(item) if item is not None else None for item in v]
                else:
                    safe_uns[k] = str(v) if v is not None else None
                    
            # Store the original uns
            original_uns = adata.uns.copy()
            # Replace with string versions for saving
            adata.uns = safe_uns
            # Save the file
            adata.write(output_file)
            # Restore original uns
            adata.uns = original_uns
            
            
        except Exception as e:
            logger.error(f"Error saving {dataset_name}: {e}")
            # Try to save with additional type safety
            try:
                # Convert all uns values to strings for safety
                safe_uns = {}
                for k, v in adata.uns.items():
                    if isinstance(v, dict):
                        safe_uns[k] = {sk: str(sv) if sv is not None else None 
                                     for sk, sv in v.items()}
                    else:
                        safe_uns[k] = str(v) if v is not None else None
                
                adata.uns = safe_uns
                adata.write(output_file)
            except Exception as e2:
                logger.error(f"Second attempt failed for {dataset_name}: {e2}")
    
    return adata

def process_all_datasets(data_dir, output_dir=None, mapping_dir=None, metadata_dir=None):
    """
    Process all datasets in a directory.
    
    Args:
        data_dir: Directory containing datasets
        output_dir: Optional directory to save standardized datasets
        mapping_dir: Optional directory containing mapping files
        metadata_dir: Optional directory containing dataset-specific metadata JSON files
    """
    data_dir = Path(data_dir)
    
    # Use same directory for output if not specified
    if not output_dir:
        output_dir = data_dir
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
    
    # Find all h5ad files
    h5ad_files = list(data_dir.glob("*_standardized_v1.h5ad"))
    
    # Also check for other h5ad files if no standardized ones found
    if not h5ad_files:
        h5ad_files = list(data_dir.glob("*.h5ad"))
    
    if not h5ad_files:
        logger.error(f"No h5ad files found in {data_dir}")
        return
    
    logger.info(f"Found {len(h5ad_files)} h5ad files to process")
    
    # Process each dataset
    for h5ad_file in h5ad_files:
        # Extract dataset name from filename
        dataset_name = h5ad_file.stem.split('_')[0]
        
        try:
            # Load the dataset
            logger.info(f"Loading {dataset_name} dataset from {h5ad_file}")
            adata = sc.read_h5ad(h5ad_file)
            
            # Determine output file path
            if "_standardized_v2.h5ad" in str(h5ad_file):
                output_file = output_dir / h5ad_file.name
            else:
                output_file = output_dir / f"{dataset_name.lower()}_standardized_v2.h5ad"
            
            # Enhance metadata standardization
            adata = enhance_standardized_metadata(dataset_name, adata, output_file, mapping_dir, metadata_dir)
            
        except Exception as e:
            logger.error(f"Error processing {dataset_name}: {e}")
            import traceback
            logger.error(traceback.format_exc())
    
    logger.info("All datasets processed successfully")

def main():
    parser = argparse.ArgumentParser(description="Standardize metadata for RNA-seq datasets")
    parser.add_argument("--data-dir", required=True, help="Directory containing h5ad files")
    parser.add_argument("--output-dir", help="Directory to save standardized datasets")
    parser.add_argument("--mapping-dir", help="Directory containing mapping files")
    parser.add_argument("--metadata-dir", default=str(DEFAULT_METADATA_DIR), 
                         help="Directory containing dataset-specific metadata JSON files")
    parser.add_argument("--dataset", help="Process only this specific dataset")
    
    args = parser.parse_args()
    
    # Check if dataset-specific processing is requested
    if args.dataset:
        # Find the dataset file
        data_dir = Path(args.data_dir)
        h5ad_files = list(data_dir.glob(f"{args.dataset.lower()}*_standardized_v1.h5ad"))
        
        # Fall back to any h5ad file matching the dataset name
        if not h5ad_files:
            h5ad_files = list(data_dir.glob(f"{args.dataset.lower()}*.h5ad"))
        
        if not h5ad_files:
            logger.error(f"No h5ad files found for dataset {args.dataset}")
            return
        
        # Use the first matching file
        h5ad_file = h5ad_files[0]
        
        # Determine output directory and file
        output_dir = Path(args.output_dir) if args.output_dir else data_dir
        output_file = output_dir / f"{args.dataset.lower()}_standardized_v2.h5ad"
        
        # Load the dataset
        logger.info(f"Loading {args.dataset} dataset from {h5ad_file}")
        adata = sc.read_h5ad(h5ad_file)
        
        # Enhance metadata standardization
        enhance_standardized_metadata(args.dataset, adata, output_file, args.mapping_dir, args.metadata_dir)
    else:
        # Process all datasets
        process_all_datasets(args.data_dir, args.output_dir, args.mapping_dir, args.metadata_dir)

if __name__ == "__main__":
    main()