#!/usr/bin/env python3
"""
Validate Croissant JSON-LD metadata files.

This script validates the generated Croissant metadata files against the 
official schema using the mlcroissant library.
"""

import mlcroissant as mlc
from pathlib import Path
import json
import sys

def validate_croissant_file(file_path: Path) -> bool:
    """Validate a single Croissant metadata file."""
    
    print(f"\nValidating: {file_path.name}")
    print(f"{'='*50}")
    
    try:
        # Load and validate the metadata
        dataset = mlc.Dataset(str(file_path))
        
        # If we get here, basic validation passed
        print("✓ Basic JSON-LD structure is valid")
        
        # Check metadata fields
        metadata = dataset.metadata
        print(f"✓ Dataset name: {metadata.name}")
        print(f"✓ Dataset description length: {len(metadata.description)} characters")
        
        if hasattr(metadata, 'citation') and metadata.citation:
            print(f"✓ Citations: {len(metadata.citation)} found")
        
        if hasattr(metadata, 'license'):
            print(f"✓ License: {metadata.license}")
        
        # Check distribution files
        if hasattr(metadata, 'distribution') and metadata.distribution:
            print(f"✓ Distribution files: {len(metadata.distribution)} found")
        
        # Check record sets
        if hasattr(metadata, 'record_set') and metadata.record_set:
            print(f"✓ Record sets: {len(metadata.record_set)} found")
        
        print("✓ Validation successful!")
        return True
        
    except Exception as e:
        print(f"✗ Validation failed: {str(e)}")
        
        # Try to load as raw JSON to check for basic JSON validity
        try:
            with open(file_path, 'r') as f:
                json.load(f)
            print("✓ File is valid JSON")
        except json.JSONDecodeError as json_e:
            print(f"✗ JSON parsing error: {str(json_e)}")
        
        return False

def main():
    """Main function to validate all Croissant files."""
    
    # Find all Croissant files in the croissant_metadata directory
    metadata_dir = Path('../croissant_metadata')
    if not metadata_dir.exists():
        metadata_dir = Path('./croissant_metadata')  # Fallback for current dir
    if not metadata_dir.exists():
        metadata_dir = Path('.')  # Final fallback
    
    croissant_files = list(metadata_dir.glob('*_croissant.jsonld'))
    
    if not croissant_files:
        print("No Croissant metadata files found")
        return
    
    print(f"Found {len(croissant_files)} Croissant metadata files to validate")
    
    valid_files = 0
    total_files = len(croissant_files)
    
    for file_path in sorted(croissant_files):
        if validate_croissant_file(file_path):
            valid_files += 1
    
    print(f"\n{'='*60}")
    print(f"Validation Summary: {valid_files}/{total_files} files passed validation")
    
    if valid_files == total_files:
        print("🎉 All Croissant metadata files are valid!")
        sys.exit(0)
    else:
        print(f"⚠️  {total_files - valid_files} files failed validation")
        sys.exit(1)

if __name__ == "__main__":
    main()