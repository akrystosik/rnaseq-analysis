#!/usr/bin/env python3
"""
Comprehensive sanity check of tissue ontology corrections before partner delivery
Validate all corrections, check for data integrity, and verify biological accuracy
"""

import pandas as pd
import anndata as ad
import json
from pathlib import Path
import requests
import time

def validate_ontology_ids_via_api(ontology_ids):
    """Validate ontology IDs using OLS4 API"""
    print("🔍 Validating ontology IDs via OLS4 API...")
    validation_results = {}
    
    for ont_id in ontology_ids:
        if not ont_id or ont_id == "":
            continue
            
        try:
            # Parse ontology type and ID
            if ont_id.startswith('UBERON:'):
                ont_type = 'uberon'
                term_id = ont_id.replace('UBERON:', '')
            elif ont_id.startswith('CL:'):
                ont_type = 'cl'
                term_id = ont_id.replace('CL:', '')
            else:
                validation_results[ont_id] = {'valid': False, 'error': 'Unknown ontology prefix'}
                continue
            
            # Query OLS4 API
            url = f"https://www.ebi.ac.uk/ols4/api/ontologies/{ont_type}/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{ont_id.replace(':', '_')}"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                validation_results[ont_id] = {
                    'valid': True,
                    'label': data.get('label', 'Unknown'),
                    'description': data.get('description', ['No description'])[0] if data.get('description') else 'No description'
                }
            else:
                validation_results[ont_id] = {'valid': False, 'error': f'API returned {response.status_code}'}
                
            time.sleep(0.1)  # Rate limiting
            
        except Exception as e:
            validation_results[ont_id] = {'valid': False, 'error': str(e)}
    
    return validation_results

def check_file_integrity(original_dir, corrected_dir):
    """Check data file integrity between original and corrected versions"""
    print("📊 Checking file integrity...")
    integrity_results = {}
    
    files_to_check = [
        'gtex_standardized_preprocessed.h5ad',
        'adni_standardized_preprocessed.h5ad', 
        'encode_standardized_preprocessed.h5ad',
        'mage_standardized_preprocessed.h5ad'
    ]
    
    for filename in files_to_check:
        orig_path = Path(original_dir) / filename
        corr_path = Path(corrected_dir) / filename
        
        if not orig_path.exists() or not corr_path.exists():
            integrity_results[filename] = {'status': 'missing_file'}
            continue
            
        try:
            # Load both files
            orig_adata = ad.read_h5ad(orig_path)
            corr_adata = ad.read_h5ad(corr_path)
            
            # Check basic dimensions
            sample_match = orig_adata.n_obs == corr_adata.n_obs
            gene_match = orig_adata.n_vars == corr_adata.n_vars
            
            # Check if gene expression data is identical (should be unchanged)
            import numpy as np
            if hasattr(orig_adata.X, 'toarray'):
                data_identical = np.allclose(orig_adata.X.toarray(), corr_adata.X.toarray())
            else:
                data_identical = np.allclose(orig_adata.X, corr_adata.X)
            
            # Check metadata columns exist
            required_cols = ['tissue', 'tissue_ontology', 'tissue_ontology_confidence']
            cols_present = all(col in corr_adata.obs.columns for col in required_cols)
            
            integrity_results[filename] = {
                'status': 'ok',
                'original_samples': int(orig_adata.n_obs),
                'corrected_samples': int(corr_adata.n_obs),
                'sample_count_match': sample_match,
                'gene_count_match': gene_match,
                'expression_data_unchanged': data_identical,
                'required_columns_present': cols_present
            }
            
        except Exception as e:
            integrity_results[filename] = {'status': 'error', 'error': str(e)}
    
    return integrity_results

def validate_specific_corrections(corrected_dir):
    """Validate that specific critical corrections were properly applied"""
    print("🎯 Validating specific corrections...")
    correction_validations = {}
    
    # Critical corrections to validate
    critical_corrections = [
        {
            'file': 'gtex_standardized_preprocessed.h5ad',
            'tissue': 'adipose - visceral (omentum)',
            'expected_ontology': 'UBERON:0014454',
            'expected_samples': 587,
            'description': "Partner's visceral adipose issue"
        },
        {
            'file': 'gtex_standardized_preprocessed.h5ad', 
            'tissue': 'brain - caudate (basal ganglia)',
            'expected_ontology': 'UBERON:0002420',
            'expected_samples': 300,
            'description': "Brain caudate correction"
        },
        {
            'file': 'mage_standardized_preprocessed.h5ad',
            'tissue': 'lymphoblast', 
            'expected_ontology': 'CL:0017005',
            'expected_samples': 731,
            'description': "MAGE lymphoblast developmental specificity"
        }
    ]
    
    for correction in critical_corrections:
        filepath = Path(corrected_dir) / correction['file']
        
        try:
            adata = ad.read_h5ad(filepath)
            
            # Find matching tissue samples
            tissue_match = adata.obs['tissue'].str.lower().str.strip() == correction['tissue']
            matching_samples = tissue_match.sum()
            
            # Check ontology assignment
            if matching_samples > 0:
                ontology_values = adata.obs[tissue_match]['tissue_ontology'].unique()
                correct_ontology = len(ontology_values) == 1 and ontology_values[0] == correction['expected_ontology']
                
                correction_validations[correction['description']] = {
                    'status': 'validated' if correct_ontology and matching_samples == correction['expected_samples'] else 'failed',
                    'expected_samples': correction['expected_samples'],
                    'found_samples': int(matching_samples),
                    'expected_ontology': correction['expected_ontology'],
                    'found_ontologies': list(ontology_values),
                    'sample_count_match': matching_samples == correction['expected_samples'],
                    'ontology_correct': correct_ontology
                }
            else:
                correction_validations[correction['description']] = {
                    'status': 'failed',
                    'error': 'No matching tissue samples found'
                }
                
        except Exception as e:
            correction_validations[correction['description']] = {
                'status': 'error',
                'error': str(e)
            }
    
    return correction_validations

def check_mapping_consistency():
    """Check that the corrected mapping file is consistent"""
    print("📝 Checking mapping file consistency...")
    
    mapping_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json"
    
    try:
        with open(mapping_file, 'r') as f:
            mappings = json.load(f)
        
        consistency_results = {
            'total_mappings': len(mappings),
            'duplicate_keys': [],
            'duplicate_values': [],
            'empty_values': [],
            'invalid_formats': []
        }
        
        # Check for duplicates and issues
        seen_keys = set()
        seen_values = {}
        
        for tissue, ontology_id in mappings.items():
            # Check duplicate keys (shouldn't happen in dict, but check anyway)
            if tissue in seen_keys:
                consistency_results['duplicate_keys'].append(tissue)
            seen_keys.add(tissue)
            
            # Check for duplicate ontology values (could be legitimate)
            if ontology_id in seen_values:
                seen_values[ontology_id].append(tissue)
            else:
                seen_values[ontology_id] = [tissue]
            
            # Check for empty values
            if not ontology_id or ontology_id == "":
                consistency_results['empty_values'].append(tissue)
            
            # Check format validity
            if ontology_id and not (ontology_id.startswith('UBERON:') or ontology_id.startswith('CL:')):
                consistency_results['invalid_formats'].append(f"{tissue}: {ontology_id}")
        
        # Identify duplicate values
        for ont_id, tissues in seen_values.items():
            if len(tissues) > 1:
                consistency_results['duplicate_values'].append(f"{ont_id}: {tissues}")
        
        consistency_results['status'] = 'ok' if not any([
            consistency_results['duplicate_keys'],
            consistency_results['empty_values'], 
            consistency_results['invalid_formats']
        ]) else 'issues_found'
        
        return consistency_results
        
    except Exception as e:
        return {'status': 'error', 'error': str(e)}

def main():
    print("🔬 COMPREHENSIVE SANITY CHECK")
    print("=" * 70)
    print("Validating corrections before partner delivery...")
    print("=" * 70)
    
    original_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/run_20250502_211754"
    corrected_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/run_20250502_211754_uberon_corrected"
    
    sanity_check_results = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'directories': {
            'original': original_dir,
            'corrected': corrected_dir
        }
    }
    
    # 1. Check file integrity
    integrity_results = check_file_integrity(original_dir, corrected_dir)
    sanity_check_results['file_integrity'] = integrity_results
    
    # 2. Validate specific critical corrections
    correction_validations = validate_specific_corrections(corrected_dir)
    sanity_check_results['critical_corrections'] = correction_validations
    
    # 3. Check mapping file consistency
    mapping_consistency = check_mapping_consistency()
    sanity_check_results['mapping_consistency'] = mapping_consistency
    
    # 4. Validate ontology IDs via API (sample of corrected IDs)
    critical_ontology_ids = ['UBERON:0014454', 'UBERON:0002420', 'CL:0017005']
    api_validation = validate_ontology_ids_via_api(critical_ontology_ids)
    sanity_check_results['ontology_validation'] = api_validation
    
    # Print summary
    print("\n📋 SANITY CHECK SUMMARY:")
    print("=" * 40)
    
    # File integrity summary
    print("\n📊 FILE INTEGRITY:")
    all_files_ok = True
    for filename, result in integrity_results.items():
        if result.get('status') == 'ok':
            intact = all([
                result.get('sample_count_match', False),
                result.get('gene_count_match', False), 
                result.get('expression_data_unchanged', False),
                result.get('required_columns_present', False)
            ])
            status = "✅" if intact else "⚠️"
            print(f"  {status} {filename}: {result.get('original_samples', 'N/A')} samples")
            if not intact:
                all_files_ok = False
        else:
            print(f"  ❌ {filename}: {result.get('error', 'Unknown error')}")
            all_files_ok = False
    
    # Critical corrections summary
    print("\n🎯 CRITICAL CORRECTIONS:")
    all_corrections_ok = True
    for description, result in correction_validations.items():
        if result.get('status') == 'validated':
            print(f"  ✅ {description}: {result['found_samples']} samples → {result['found_ontologies'][0]}")
        else:
            print(f"  ❌ {description}: {result.get('error', 'Validation failed')}")
            all_corrections_ok = False
    
    # Mapping consistency summary  
    print("\n📝 MAPPING CONSISTENCY:")
    mapping_ok = mapping_consistency.get('status') == 'ok'
    status = "✅" if mapping_ok else "❌"
    print(f"  {status} {mapping_consistency.get('total_mappings', 0)} mappings")
    if not mapping_ok:
        for issue_type in ['empty_values', 'invalid_formats', 'duplicate_keys']:
            if mapping_consistency.get(issue_type):
                print(f"    - {issue_type}: {len(mapping_consistency[issue_type])}")
    
    # API validation summary
    print("\n🔍 ONTOLOGY VALIDATION:")
    all_valid = all(result.get('valid', False) for result in api_validation.values())
    status = "✅" if all_valid else "❌"
    print(f"  {status} {len([r for r in api_validation.values() if r.get('valid')])} / {len(api_validation)} IDs valid")
    for ont_id, result in api_validation.items():
        if result.get('valid'):
            print(f"    ✅ {ont_id}: {result.get('label', 'Unknown')}")
        else:
            print(f"    ❌ {ont_id}: {result.get('error', 'Unknown error')}")
    
    # Overall assessment
    overall_pass = all_files_ok and all_corrections_ok and mapping_ok and all_valid
    
    print("\n" + "=" * 70)
    if overall_pass:
        print("🎉 SANITY CHECK PASSED - READY FOR PARTNER DELIVERY")
        print("✅ All corrections validated")
        print("✅ Data integrity preserved") 
        print("✅ Ontology IDs verified")
        print("✅ Partner's critical issue resolved")
    else:
        print("⚠️  SANITY CHECK FOUND ISSUES - REVIEW REQUIRED")
        if not all_files_ok:
            print("❌ File integrity issues detected")
        if not all_corrections_ok:
            print("❌ Critical correction validation failed")
        if not mapping_ok:
            print("❌ Mapping consistency issues found")
        if not all_valid:
            print("❌ Ontology ID validation failed")
    
    print("=" * 70)
    
    # Save detailed results
    output_file = Path(corrected_dir) / "sanity_check_results.json"
    with open(output_file, 'w') as f:
        json.dump(sanity_check_results, f, indent=2, default=str)
    
    print(f"\n📁 Detailed results saved: {output_file}")
    
    return overall_pass

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)