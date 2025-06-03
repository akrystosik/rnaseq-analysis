#!/usr/bin/env python3
"""
Create mapping from ENCODE ENCFF IDs to cell line names
"""

import anndata as ad
import pandas as pd
import json

def create_encff_to_cellline_mapping():
    """Create mapping from ENCFF IDs to cell line names"""
    
    # Load the ENCODE standardized H5AD file
    encode_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/run_20250528_235853/encode_standardized_v2.h5ad"
    adata = ad.read_h5ad(encode_file)
    
    # Create mapping from H5AD data
    h5ad_mapping = {}
    for idx, row in adata.obs.iterrows():
        encff_id = str(idx)  # Index contains ENCFF IDs
        cell_line = row['cell_line']
        h5ad_mapping[encff_id] = cell_line
    
    print("=== ENCFF to Cell Line Mapping from H5AD ===")
    for encff_id, cell_line in h5ad_mapping.items():
        print(f"{encff_id} -> {cell_line}")
    
    # Load metadata JSON to cross-reference
    metadata_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/encode_metadata.json"
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
    
    # Extract file mappings from metadata
    print("\n=== ENCFF to Cell Line Mapping from Metadata JSON ===")
    metadata_mapping = {}
    if 'dataset_info' in metadata and 'cell_lines' in metadata['dataset_info']:
        for cell_line, info in metadata['dataset_info']['cell_lines'].items():
            if 'file' in info and 'accession' in info['file']:
                encff_id = info['file']['accession']
                metadata_mapping[encff_id] = cell_line
                print(f"{encff_id} -> {cell_line}")
    
    # Combine both mappings (H5AD takes precedence as it's what was actually processed)
    final_mapping = {**metadata_mapping, **h5ad_mapping}
    
    print("\n=== Final Combined ENCFF to Cell Line Mapping ===")
    for encff_id in sorted(final_mapping.keys()):
        cell_line = final_mapping[encff_id]
        print(f"{encff_id} -> {cell_line}")
    
    # Save mapping to JSON file
    output_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/encff_to_cellline_mapping.json"
    with open(output_file, 'w') as f:
        json.dump(final_mapping, f, indent=2)
    
    print(f"\nMapping saved to: {output_file}")
    
    # Create reverse mapping (cell line to ENCFF IDs)
    cellline_to_encff = {}
    for encff_id, cell_line in final_mapping.items():
        if cell_line not in cellline_to_encff:
            cellline_to_encff[cell_line] = []
        cellline_to_encff[cell_line].append(encff_id)
    
    print("\n=== Cell Line to ENCFF IDs (Reverse Mapping) ===")
    for cell_line in sorted(cellline_to_encff.keys()):
        encff_ids = cellline_to_encff[cell_line]
        print(f"{cell_line} -> {', '.join(encff_ids)}")
    
    # Save reverse mapping
    reverse_output_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2/cellline_to_encff_mapping.json"
    with open(reverse_output_file, 'w') as f:
        json.dump(cellline_to_encff, f, indent=2)
    
    print(f"\nReverse mapping saved to: {reverse_output_file}")
    
    return final_mapping, cellline_to_encff

if __name__ == "__main__":
    create_encff_to_cellline_mapping()