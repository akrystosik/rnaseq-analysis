#!/usr/bin/env python3

import pandas as pd
import anndata as ad

# Load the original ENCODE data
print("Loading original ENCODE data...")
adata_orig = ad.read_h5ad('../../standardized_data/encode_standardized_v2.h5ad')

# Load reference mapping
print("Loading reference mapping...")
mapping_df = pd.read_csv('../../metadata/gene_mapping/gene_id_reference_mapping.csv')
numeric_to_ensembl = {}
ensembl_to_info = {}

for _, row in mapping_df.iterrows():
    gene_id = row['gene_id']
    numeric_id = row['numeric_id']
    
    if not pd.isna(numeric_id):
        numeric_to_ensembl[str(numeric_id)] = gene_id
    
    ensembl_to_info[gene_id] = {
        'gene_name': row['gene_name'],
        'gene_type': row['gene_type'],
        'chromosome': row['chromosome'],
        'mapping_confidence': row['mapping_confidence']
    }

# Load ENCODE mapping
print("Loading ENCODE mapping...")
encode_df = pd.read_csv('../../metadata/gene_mapping/encode_specific_mappings/encode_id_to_ensembl_mapping.csv')
encode_orig_to_target = dict(zip(encode_df['original_id'].astype(str), encode_df['ensembl_id'].astype(str)))

print(f"Reference numeric mapping: {len(numeric_to_ensembl)} entries")
print(f"ENCODE mapping: {len(encode_orig_to_target)} entries")

# Test specific problematic IDs
test_ids = ['1', '2', '9', '10', '0', '3', '4', '5']

for test_idx in test_ids:
    print(f"\n=== Processing ID {test_idx} ===")
    
    # Simulate the original data structure
    input_var_name = test_idx
    original_id_source = ''  # ENCODE doesn't have original_ids column
    
    print(f"  input_var_name: '{input_var_name}'")
    print(f"  original_id_source: '{original_id_source}'")
    
    target_id = ''
    mapping_source = 'unmapped'
    
    # Step 1: Try ENCODE mapping
    if original_id_source in encode_orig_to_target:
        target_id = encode_orig_to_target[original_id_source]
        mapping_source = 'encode_mapping'
        print(f"  Step 1 - ENCODE mapping: '{original_id_source}' -> '{target_id}'")
    else:
        print(f"  Step 1 - ENCODE mapping: '{original_id_source}' not found")
    
    # Step 2: If ENCODE gave ENTREZ, try reference
    if target_id.startswith('ENTREZ:'):
        entrez_numeric = target_id[7:]
        if entrez_numeric in numeric_to_ensembl:
            potential_ensembl = numeric_to_ensembl[entrez_numeric]
            if potential_ensembl.startswith('ENSG'):
                target_id = potential_ensembl
                mapping_source = 'encode_mapping -> ref_numeric'
                print(f"  Step 2 - ENTREZ upgrade: ENTREZ:{entrez_numeric} -> {potential_ensembl}")
    
    # Step 3: Reference numeric fallback
    if not target_id.startswith('ENSG') and input_var_name.isdigit():
        print(f"  Step 3 condition met: target_id='{target_id}', isdigit={input_var_name.isdigit()}")
        if input_var_name in numeric_to_ensembl:
            potential_ensembl = numeric_to_ensembl[input_var_name]
            print(f"  Step 3 - Found in reference: {input_var_name} -> {potential_ensembl}")
            if potential_ensembl.startswith('ENSG'):
                target_id = potential_ensembl
                mapping_source = 'ref_numeric_fallback'
                print(f"  Step 3 - SUCCESS: Using {potential_ensembl}")
            else:
                target_id = f"ENTREZ:{input_var_name}"
                mapping_source = 'entrez_id_fallback'
                print(f"  Step 3 - Non-ENSG reference: Using ENTREZ:{input_var_name}")
        else:
            target_id = f"ENTREZ:{input_var_name}"
            mapping_source = 'entrez_id_fallback'
            print(f"  Step 3 - Not in reference: Using ENTREZ:{input_var_name}")
    else:
        print(f"  Step 3 condition NOT met: target_id='{target_id}', isdigit={input_var_name.isdigit()}")
    
    # Final fallback
    if not target_id.startswith(('ENSG', 'ENTREZ:', 'gSpikein')):
        target_id = input_var_name
        mapping_source = 'unmapped'
        print(f"  Final fallback: Using {input_var_name}")
    
    print(f"  FINAL RESULT: '{target_id}' (source: {mapping_source})")