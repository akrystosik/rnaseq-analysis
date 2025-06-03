import json
import os

# Test ENCODE metadata
encode_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/encode_metadata.json'
print(f"Testing ENCODE metadata at {encode_path}")
print(f"File exists: {os.path.exists(encode_path)}")
if os.path.exists(encode_path):
    try:
        with open(encode_path, 'r') as f:
            data = json.load(f)
        print(f"Successfully loaded encode_metadata.json with {len(data)} top-level keys")
    except Exception as e:
        print(f"Error loading encode_metadata.json: {e}")

# Test ENTEx metadata
entex_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/metadata/entex_metadata.json'
print(f"\nTesting ENTEx metadata at {entex_path}")
print(f"File exists: {os.path.exists(entex_path)}")
if os.path.exists(entex_path):
    try:
        with open(entex_path, 'r') as f:
            data = json.load(f)
        print(f"Successfully loaded entex_metadata.json with {len(data)} top-level keys")
        if 'sample_lookup' in data:
            print(f"Found 'sample_lookup' with {len(data['sample_lookup'])} samples")
        if 'donor_map' in data:
            print(f"Found 'donor_map' with {len(data['donor_map'])} donors")
    except Exception as e:
        print(f"Error loading entex_metadata.json: {e}")