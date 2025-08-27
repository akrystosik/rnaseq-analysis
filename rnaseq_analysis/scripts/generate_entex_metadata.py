#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/generate_entex_metadata.py
#!/usr/bin/env python3
"""
Script to generate comprehensive metadata for ENTEx files
and integrate with existing metadata for standardization.
"""

import os
import json
import argparse
import pandas as pd
import re
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Generate comprehensive metadata for ENTEx files')
    parser.add_argument('--input-dir', '-i', required=True, help='Directory containing downloaded ENTEx files')
    parser.add_argument('--metadata-file', '-m', required=True, help='Path to raw metadata JSON file')
    parser.add_argument('--output-file', '-o', required=True, help='Path to output standardized metadata JSON file')
    parser.add_argument('--existing-metadata', '-e', help='Path to existing metadata file to update/merge with')
    return parser.parse_args()

def read_tsv_metadata(file_path):
    """Extract metadata from TSV gene quantification file."""
    try:
        # Read the first few lines to extract metadata
        metadata = {}
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                if i >= 10:  # Check more lines for metadata
                    break
                if line.startswith('#'):
                    parts = line.strip('#').strip().split(':', 1)
                    if len(parts) == 2:
                        key, value = parts
                        metadata[key.strip()] = value.strip()
        
        # Try to extract gene ID type from file content
        df = pd.read_csv(file_path, sep='\t', comment='#', nrows=5)
        
        # Check for transcript vs gene quantification
        if 'transcript_id' in df.columns:
            metadata['quantification_level'] = 'transcript'
        elif 'gene_id' in df.columns:
            metadata['quantification_level'] = 'gene'
        else:
            # Try to infer from column contents
            for col in df.columns:
                if 'gene' in col.lower():
                    metadata['quantification_level'] = 'gene'
                    break
                elif 'transcript' in col.lower() or 'tx' in col.lower():
                    metadata['quantification_level'] = 'transcript'
                    break
            
            if 'quantification_level' not in metadata:
                # Try to infer from filename
                if 'gene' in file_path.lower() and not 'transcript' in file_path.lower():
                    metadata['quantification_level'] = 'gene'
                elif 'transcript' in file_path.lower():
                    metadata['quantification_level'] = 'transcript'
                else:
                    metadata['quantification_level'] = 'unknown'
        
        # Check gene ID type if we have gene_id column
        if 'gene_id' in df.columns:
            # Check first few gene_ids to determine format
            sample_ids = df['gene_id'].iloc[:5].tolist()
            if any('ENSG' in str(gene_id) for gene_id in sample_ids):
                metadata['gene_id_type'] = 'ensembl'
            elif any('ENTREZ' in str(gene_id) for gene_id in sample_ids):
                metadata['gene_id_type'] = 'entrez'
            else:
                metadata['gene_id_type'] = 'unknown'
        
        return metadata
    except Exception as e:
        print(f"Error reading TSV metadata from {file_path}: {e}")
        return {}

def create_standardized_metadata(raw_metadata, file_dir):
    """Create standardized metadata from raw metadata and file information."""
    standardized_metadata = {}
    
    for file_id, meta in raw_metadata.items():
        # Find the actual file path
        file_path = None
        for filename in os.listdir(file_dir):
            if file_id in filename:
                file_path = os.path.join(file_dir, filename)
                break
        
        if not file_path or not os.path.exists(file_path):
            print(f"File not found for {file_id}, skipping")
            continue
        
        # Extract additional metadata from TSV
        tsv_metadata = read_tsv_metadata(file_path)
        
        # Standardize tissue name
        tissue = standardize_tissue_name(meta.get('tissue', 'unknown_tissue'))
        
        # Check if this is a gene quantification (not transcript)
        quant_level = tsv_metadata.get('quantification_level', 'unknown')
        if quant_level != 'gene':
            print(f"Skipping non-gene quantification file: {file_id} (type: {quant_level})")
            continue
        
        # Create standardized entry
        entry = {
            'file_id': file_id,
            'file_path': file_path,
            'experiment_id': meta.get('experiment_id', ''),
            'tissue': tissue,
            'donor': meta.get('donor', 'unknown_donor'),
            'assay_type': meta.get('assay_type', 'RNA-seq'),
            'dataset_type': 'entex',
            'genome_assembly': meta.get('genome_assembly', 'GRCh38'),
            'quantification_method': meta.get('output_type', 'gene quantifications'),
            'quantification_level': quant_level,  # Add explicit field for gene vs transcript
            'biological_replicate': meta.get('biological_replicates', [0])[0] if meta.get('biological_replicates') else 0,
            'technical_replicate': meta.get('technical_replicates', [0])[0] if meta.get('technical_replicates') else 0,
            'gene_id_type': tsv_metadata.get('gene_id_type', 'ensembl'),
            'expression_unit': 'TPM'  # Assuming TPM for ENCODE data
        }
        
        standardized_metadata[file_id] = entry
    
    return standardized_metadata

def standardize_tissue_name(tissue):
    """Standardize tissue names for consistent mapping."""
    tissue = tissue.lower()
    
    # Map of raw tissue names to standardized names
    tissue_map = {
        'gastrocnemius medialis': 'muscle',
        'thyroid gland': 'thyroid',
        'adrenal gland': 'adrenal',
        'stomach': 'stomach',
        'upper lobe of left lung': 'lung',
        'transverse colon': 'colon',
        # Add more mappings as needed
    }
    
    return tissue_map.get(tissue, tissue)

def create_standardized_metadata(raw_metadata, file_dir):
    """Create standardized metadata from raw metadata and file information."""
    standardized_metadata = {}
    
    for file_id, meta in raw_metadata.items():
        # Find the actual file path
        file_path = None
        for filename in os.listdir(file_dir):
            if file_id in filename:
                file_path = os.path.join(file_dir, filename)
                break
        
        if not file_path or not os.path.exists(file_path):
            print(f"File not found for {file_id}, skipping")
            continue
        
        # Extract additional metadata from TSV
        tsv_metadata = read_tsv_metadata(file_path)
        
        # Standardize tissue name
        tissue = standardize_tissue_name(meta.get('tissue', 'unknown_tissue'))
        
        # Create standardized entry
        entry = {
            'file_id': file_id,
            'file_path': file_path,
            'experiment_id': meta.get('experiment_id', ''),
            'tissue': tissue,
            'donor': meta.get('donor', 'unknown_donor'),
            'assay_type': meta.get('assay_type', 'RNA-seq'),
            'dataset_type': 'entex',
            'genome_assembly': meta.get('genome_assembly', 'GRCh38'),
            'quantification_method': meta.get('output_type', 'gene quantifications'),
            'biological_replicate': meta.get('biological_replicates', [0])[0] if meta.get('biological_replicates') else 0,
            'technical_replicate': meta.get('technical_replicates', [0])[0] if meta.get('technical_replicates') else 0,
            'gene_id_type': tsv_metadata.get('gene_id_type', 'ensembl'),
            'expression_unit': 'TPM'  # Assuming TPM for ENCODE data
        }
        
        standardized_metadata[file_id] = entry
    
    return standardized_metadata

def update_existing_metadata(existing_metadata, new_metadata):
    """Update existing metadata with new metadata, preserving existing entries."""
    if not existing_metadata:
        return new_metadata
    
    # Merge dictionaries, new_metadata takes precedence for duplicates
    merged = existing_metadata.copy()
    for key, value in new_metadata.items():
        merged[key] = value
    
    return merged

def main():
    args = parse_args()
    
    # Read raw metadata
    with open(args.metadata_file, 'r') as f:
        raw_metadata = json.load(f)
    
    # Create standardized metadata
    standardized_metadata = create_standardized_metadata(raw_metadata, args.input_dir)
    
    # Update with existing metadata if provided
    if args.existing_metadata and os.path.exists(args.existing_metadata):
        with open(args.existing_metadata, 'r') as f:
            existing_metadata = json.load(f)
        standardized_metadata = update_existing_metadata(existing_metadata, standardized_metadata)
    
    # Save standardized metadata
    with open(args.output_file, 'w') as f:
        json.dump(standardized_metadata, f, indent=2)
    
    print(f"Generated metadata for {len(standardized_metadata)} files")
    print(f"Standardized metadata saved to {args.output_file}")

if __name__ == "__main__":
    main()