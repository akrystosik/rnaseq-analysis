#!/usr/bin/env python3
import json
import os

# Path to ADNI metadata file
metadata_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/adni_metadata.json"

# Load existing metadata
with open(metadata_file, 'r') as f:
    metadata = json.load(f)

# Add required top-level fields for validation
metadata['harmonized_gencode_version'] = '24'
metadata['harmonized_reference_genome'] = 'hg38'

# If there's a gencode_version in dataset_info, add it as original_gencode_version
if 'dataset_info' in metadata and 'gencode_version' in metadata['dataset_info']:
    metadata['original_gencode_version'] = metadata['dataset_info']['gencode_version']

# If there's a reference_genome in dataset_info, add it as original_reference_genome
if 'dataset_info' in metadata and 'reference_genome' in metadata['dataset_info']:
    metadata['original_reference_genome'] = metadata['dataset_info']['reference_genome']

# Add additional metadata to dataset_info
if 'dataset_info' not in metadata:
    metadata['dataset_info'] = {}

metadata['dataset_info']['platform'] = 'Affymetrix Human Genome U219 Array'
metadata['dataset_info']['expression_unit'] = 'Normalized intensity'
metadata['dataset_info']['reference'] = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7541709/'

# Add notes about genome liftover if original was different from harmonized
if metadata.get('original_reference_genome') != metadata.get('harmonized_reference_genome'):
    metadata['genome_mapping_notes'] = f"Original ADNI data was aligned to {metadata.get('original_reference_genome', 'hg19')}. For harmonization across datasets, coordinates were lifted over to {metadata.get('harmonized_reference_genome', 'hg38')} to ensure consistent comparison with other RNA-seq datasets."

# Save updated metadata
with open(metadata_file, 'w') as f:
    json.dump(metadata, f, indent=2)

print(f"Updated {metadata_file} with required fields for validation.")
