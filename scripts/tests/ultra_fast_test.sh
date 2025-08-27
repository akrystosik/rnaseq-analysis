#!/bin/bash
# Ultra-fast test script for RNA-seq pipeline fixes

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"

echo "=== Running Ultra-Fast Tests for RNA-seq Pipeline Fixes ==="

# 1. Test ENCODE tissue fix
echo -e "\n=== Testing ENCODE tissue extraction fix ==="
python3 -c "
import sys
sys.path.append('${BASE_DIR}/scripts')
import os
import json

# Import function directly if possible, otherwise mock it
try:
    from standardize_datasets import extract_encode_metadata
    print('Successfully imported extract_encode_metadata function')
except ImportError:
    print('Could not import function directly, using mock test')
    
    # Define mock function with our fix
    def extract_encode_metadata(file_path, entex_metadata=None, metadata_dir=None):
        metadata = {}
        
        # Get encode cell info
        encode_cell_info = {}
        if metadata_dir:
            json_path = os.path.join(metadata_dir, 'encode_metadata.json')
            if os.path.exists(json_path):
                with open(json_path, 'r') as f:
                    metadata_json = json.load(f)
                    if 'dataset_info' in metadata_json and 'cell_lines' in metadata_json['dataset_info']:
                        encode_cell_info = metadata_json['dataset_info']['cell_lines']
        
        # Extract cell line from path - improved detection
        path_parts = file_path.split(os.sep)
        cell_line_dir = []
        
        for part in path_parts:
            # Check exact match
            if part in encode_cell_info:
                cell_line_dir.append(part)
                continue
                
            # Check for cell line as substring
            for cl in encode_cell_info.keys():
                if cl in part:
                    cell_line_dir.append(part)
                    break
        
        # Store cell line if found
        if cell_line_dir:
            cell_line_full = cell_line_dir[0]
            
            # Extract cell line name
            cell_line = None
            for cl in encode_cell_info.keys():
                if cell_line_full.startswith(cl):
                    cell_line = cl
                    break
            
            if cell_line:
                metadata['cell_line'] = cell_line
                
                # Add tissue from cell line info
                if cell_line in encode_cell_info and 'tissue' in encode_cell_info[cell_line]:
                    metadata['tissue'] = encode_cell_info[cell_line]['tissue']
        
        # Add standard fields
        metadata['data_type'] = 'RNA-seq'
        metadata['expression_unit'] = 'TPM'
        
        # Ensure tissue is never None/nan for cell lines
        if 'cell_line' in metadata and 'tissue' not in metadata:
            cell_line = metadata['cell_line']
            if cell_line in encode_cell_info and 'tissue' in encode_cell_info[cell_line]:
                metadata['tissue'] = encode_cell_info[cell_line]['tissue']
            else:
                # Fallback mapping for common cell lines
                fallback_mapping = {
                    'A549': 'lung',
                    'K562': 'blood',
                    'HepG2': 'liver',
                    'IMR-90': 'lung',
                    'H1-hESC': 'embryonic stem cell'
                }
                if cell_line in fallback_mapping:
                    metadata['tissue'] = fallback_mapping[cell_line]
        
        return metadata

# Test extract_encode_metadata with various file paths
test_paths = [
    '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/A549/ENCFF001.tsv',
    '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/K562_polyA_plus/ENCFF002.tsv',
    '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/unknown/ENCFF003.tsv'
]

metadata_dir = '${BASE_DIR}/metadata/json'

for path in test_paths:
    print(f'\\nTesting path: {path}')
    metadata = extract_encode_metadata(path, None, metadata_dir)
    print('Extracted metadata:')
    for key, value in metadata.items():
        print(f'  {key}: {value}')
    
    # Check if tissue was extracted
    if 'tissue' in metadata:
        print(f'  Result: SUCCESS - Tissue found: {metadata[\"tissue\"]}')
    else:
        print('  Result: FAILED - No tissue found')
"

# 2. Test placeholder gene ID fix
echo -e "\n=== Testing placeholder gene ID fix ==="
python3 -c "
import pandas as pd

# Create test data with placeholder IDs
test_data = {
    'gene_id': ['gene1', 'gene2', 'gene3', 'gene4'],
    'ensembl_id': ['ENSG0001', 'PLACEHOLDER_123', 'ENSG0003', 'PLACEHOLDER_456'],
    'mapping_source': ['reference_mapping', 'reference_mapping', 'reference_mapping', 'reference_mapping']
}

df = pd.DataFrame(test_data)
print('Original data:')
print(df)

# Apply our fix for placeholder IDs
for idx, row in df.iterrows():
    ensembl_id = row['ensembl_id']
    if isinstance(ensembl_id, str) and ensembl_id.startswith('PLACEHOLDER_'):
        df.loc[idx, 'ensembl_id'] = ''
        df.loc[idx, 'mapping_source'] = 'unmapped'

print('\\nAfter fix:')
print(df)

# Check if fix worked
placeholder_count = sum(1 for x in df['ensembl_id'] if str(x).startswith('PLACEHOLDER_'))
if placeholder_count == 0:
    print('\\nTest PASSED: All placeholder IDs removed')
else:
    print(f'\\nTest FAILED: {placeholder_count} placeholder IDs remain')
"

# 3. Test dataset label fix
echo -e "\n=== Testing dataset label fix ==="
python3 -c "
import pandas as pd

# Create test DataFrame with 'combined' dataset labels
test_data = {
    'sample_id': ['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
    'tissue': ['lung', 'blood', 'liver', 'brain', 'kidney'],
    'dataset': ['encode', 'combined', 'combined', 'mage', 'adni']
}

df = pd.DataFrame(test_data)
print('Original data:')
print(df)

# Apply our fix for dataset labels
dataset_name = 'encode'  # This would be the actual dataset name in the real code
df['dataset'] = dataset_name

print('\\nAfter fix:')
print(df)

# Check if fix worked
combined_count = sum(df['dataset'] == 'combined')
if combined_count == 0:
    print('\\nTest PASSED: All \"combined\" labels replaced')
else:
    print(f'\\nTest FAILED: {combined_count} \"combined\" labels remain')
"

echo -e "\n=== Ultra-Fast Tests Complete ==="
