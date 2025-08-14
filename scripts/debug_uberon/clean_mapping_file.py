#!/usr/bin/env python3
"""
Clean the tissue_to_uberon.json file:
1. Keep only the first 66 validated entries (lines 1-66)  
2. Convert all keys to lowercase for consistent lookups
3. Remove the redundant/inconsistent lines 67-133
"""

import json

def main():
    # Load current mapping file
    input_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json'
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    print("Original mapping file analysis:")
    print(f"Total entries: {len(data)}")
    
    # Get all keys and take only first 66 (the validated ones)
    all_keys = list(data.keys())
    first_66_keys = all_keys[:66]
    
    print(f"Keeping first 66 validated entries")
    print(f"Removing {len(all_keys) - 66} redundant/inconsistent entries")
    
    # Create new mapping with lowercase keys from first 66 entries
    cleaned_mapping = {}
    for key in first_66_keys:
        lowercase_key = key.lower()
        cleaned_mapping[lowercase_key] = data[key]
    
    print(f"\nCleaned mapping:")
    print(f"Total entries: {len(cleaned_mapping)}")
    
    # Show first few entries as example
    print(f"\nFirst 5 entries (now lowercase):")
    for i, (key, value) in enumerate(list(cleaned_mapping.items())[:5], 1):
        print(f"  {i}. \"{key}\": \"{value}\"")
    
    # Save cleaned mapping
    output_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon_cleaned.json'
    with open(output_file, 'w') as f:
        json.dump(cleaned_mapping, f, indent=2, sort_keys=True)
    
    print(f"\n✅ Cleaned mapping saved to: tissue_to_uberon_cleaned.json")
    
    # Show summary of our validated corrections that are now included
    corrections_included = [
        ('brain - caudate (basal ganglia)', 'UBERON:0002420'),
        ('brain - hippocampus', 'UBERON:0002421'),
        ('esophagus - muscularis', 'UBERON:0003832'),
        ('heart - atrial appendage', 'UBERON:0006618'),
        ('lymphoblast', 'CL:0017005'),
        ('pbmc', 'CL:2000001')
    ]
    
    print(f"\n✅ Validated corrections included:")
    for tissue, mapping in corrections_included:
        if tissue in cleaned_mapping and cleaned_mapping[tissue] == mapping:
            print(f"  ✓ {tissue}: {mapping}")
        else:
            print(f"  ⚠ {tissue}: Expected {mapping}, got {cleaned_mapping.get(tissue, 'NOT_FOUND')}")

if __name__ == "__main__":
    main()