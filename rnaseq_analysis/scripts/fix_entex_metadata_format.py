#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/fix_entex_metadata_format.py
#!/usr/bin/env python3
"""
Script to fix the format of ENTEx metadata JSON file by consolidating
the different sections into a unified format with proper file paths.
"""

import os
import json
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Fix ENTEx metadata format')
    parser.add_argument('--metadata', required=True, 
                        help='Path to ENTEx metadata JSON file')
    parser.add_argument('--data-dir', required=True,
                        help='Directory containing ENTEx data files')
    parser.add_argument('--output', 
                        help='Path to save updated metadata (defaults to overwriting input)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load the metadata
    print(f"Loading metadata from {args.metadata}")
    with open(args.metadata, 'r') as f:
        metadata = json.load(f)
    
    # Check if file entries exist in the data directory
    data_dir = Path(args.data_dir)
    
    # Create the new combined metadata
    combined_metadata = {}
    
    # Track the number of entries processed
    entex_entries_processed = 0
    file_entries_processed = 0
    
    # Process existing file entries (like ENCFF977PCO)
    for key, value in metadata.items():
        # Skip the donor_map and entex_metadata keys
        if key in ['donor_map', 'entex_metadata', 'sample_lookup']:
            continue
            
        if isinstance(value, dict) and 'file_id' in value:
            # This is a file entry
            file_id = value['file_id']
            file_path = value.get('file_path')
            
            # Check if the file exists
            if file_path and os.path.exists(file_path):
                # Keep the entry with existing path
                combined_metadata[file_id] = value
                file_entries_processed += 1
            else:
                # Try to find the file in the data directory
                possible_path = data_dir / f"{file_id}.tsv"
                if possible_path.exists():
                    value['file_path'] = str(possible_path)
                    combined_metadata[file_id] = value
                    file_entries_processed += 1
                else:
                    print(f"Warning: File not found for {file_id}")
    
    # Process entries from entex_metadata array
    if 'entex_metadata' in metadata:
        for entry in metadata['entex_metadata']:
            sample_id = entry.get('sample_id')
            
            if not sample_id:
                continue
                
            # Check if this entry is already in combined_metadata
            if sample_id in combined_metadata:
                continue
                
            # Try to find file
            possible_path = data_dir / f"{sample_id}.tsv"
            if possible_path.exists():
                # Create a new entry with proper file path
                new_entry = entry.copy()
                new_entry['file_id'] = sample_id
                new_entry['file_path'] = str(possible_path)
                
                # Add donor information if available
                donor_id = entry.get('donor_id')
                if donor_id and 'donor_map' in metadata and donor_id in metadata['donor_map']:
                    donor_info = metadata['donor_map'][donor_id]
                    if 'sex' in donor_info:
                        new_entry['sex'] = donor_info['sex']
                    if 'age' in donor_info:
                        new_entry['age'] = donor_info['age']
                
                # Add dataset info
                new_entry['dataset_type'] = 'entex'
                new_entry['data_type'] = 'RNA-seq'
                new_entry['expression_unit'] = 'TPM'
                
                combined_metadata[sample_id] = new_entry
                entex_entries_processed += 1
            else:
                print(f"Warning: File not found for entex_metadata entry {sample_id}")
    
    # Keep donor_map for reference
    if 'donor_map' in metadata:
        combined_metadata['donor_map'] = metadata['donor_map']
    
    print(f"Processed {file_entries_processed} file entries and {entex_entries_processed} entries from entex_metadata")
    print(f"Total entries in new metadata: {len(combined_metadata) - ('donor_map' in combined_metadata)}")
    
    # Save updated metadata
    output_path = args.output if args.output else args.metadata
    print(f"Saving updated metadata to {output_path}")
    with open(output_path, 'w') as f:
        json.dump(combined_metadata, f, indent=2)
    
    print("Metadata format has been fixed. You can now run the standardization pipeline.")

if __name__ == "__main__":
    main()