#!/usr/bin/env python3
"""
Systematic validation of tissue-to-UBERON mappings using OLS4 API
"""

import json
import requests
import time
from urllib.parse import quote

def validate_uberon_term(uberon_id):
    """Validate a UBERON term via OLS4 API"""
    if uberon_id.startswith('CL:'):
        # Cell Ontology terms - skip for now
        return {'valid': True, 'label': 'Cell Ontology term', 'description': 'Cell type'}
    
    if not uberon_id.startswith('UBERON:'):
        return {'valid': False, 'error': f'Invalid UBERON format: {uberon_id}'}
    
    # Convert UBERON:0000955 to IRI format
    iri = f"http://purl.obolibrary.org/obo/{uberon_id.replace(':', '_')}"
    double_encoded_iri = quote(quote(iri, safe=''), safe='')
    
    url = f"https://www.ebi.ac.uk/ols4/api/ontologies/uberon/terms/{double_encoded_iri}"
    
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return {
                'valid': True,
                'label': data.get('label', 'Unknown'),
                'description': data.get('description', ['No description'])[0] if data.get('description') else 'No description'
            }
        elif response.status_code == 404:
            return {'valid': False, 'error': f'UBERON term not found: {uberon_id}'}
        else:
            return {'valid': False, 'error': f'API error {response.status_code} for {uberon_id}'}
    except Exception as e:
        return {'valid': False, 'error': f'Request failed for {uberon_id}: {str(e)}'}

def check_semantic_accuracy(tissue_name, uberon_details):
    """Check if tissue name semantically matches UBERON term"""
    if not uberon_details.get('valid'):
        return True, []
    
    tissue_lower = tissue_name.lower()
    label_lower = uberon_details['label'].lower()
    
    issues = []
    
    # Check specific problematic cases we've identified
    if 'terminal ileum' in tissue_lower and 'peyer' in label_lower:
        issues.append(f"Terminal ileum tissue '{tissue_name}' mapped to Peyer's patch '{uberon_details['label']}'")
    elif 'skin' in tissue_lower and 'placenta' in label_lower:
        issues.append(f"Skin tissue '{tissue_name}' mapped to placenta '{uberon_details['label']}'")
    elif 'suprapubic' in tissue_lower and 'suprapubic' not in label_lower:
        issues.append(f"Suprapubic tissue '{tissue_name}' not mapped to suprapubic term '{uberon_details['label']}'")
    elif 'lower leg' in tissue_lower and 'lower leg' not in label_lower and 'leg' not in label_lower:
        issues.append(f"Lower leg tissue '{tissue_name}' not mapped to leg-related term '{uberon_details['label']}'")
    elif 'caudate (basal ganglia)' in tissue_lower and 'caudate nucleus' in label_lower:
        issues.append(f"Caudate (basal ganglia) region '{tissue_name}' mapped to specific nucleus '{uberon_details['label']}' instead of broader striatal region")
    
    return len(issues) == 0, issues

def main():
    # Load the tissue mapping file
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json', 'r') as f:
        mappings = json.load(f)
    
    print("Comprehensive Tissue-to-UBERON Mapping Validation")
    print("=" * 60)
    print(f"Total mappings: {len(mappings)}")
    print()
    
    errors = []
    semantic_issues = []
    validated = 0
    processed_ids = set()
    
    # Check each mapping for validity AND semantic accuracy
    for tissue_name, uberon_id in mappings.items():
        if uberon_id in processed_ids:
            continue
        processed_ids.add(uberon_id)
        
        result = validate_uberon_term(uberon_id)
        
        if result['valid']:
            print(f"✓ {uberon_id}: {result['label']}")
            validated += 1
            
            # Check semantic accuracy
            is_accurate, issues = check_semantic_accuracy(tissue_name, result)
            if not is_accurate:
                print(f"  ⚠️  SEMANTIC ISSUES:")
                for issue in issues:
                    print(f"    {issue}")
                    semantic_issues.append(issue)
                
                # Show all affected tissues
                affected_tissues = [t for t, id in mappings.items() if id == uberon_id]
                print(f"    Affects tissues: {', '.join(affected_tissues)}")
        else:
            print(f"✗ {uberon_id}: {result['error']}")
            errors.append((uberon_id, result['error']))
        
        print()
        time.sleep(0.1)  # Be nice to the API
    
    print("=" * 60)
    print("VALIDATION SUMMARY")
    print(f"- Total unique UBERON IDs: {len(processed_ids)}")
    print(f"- Valid terms: {validated}")
    print(f"- Invalid terms: {len(errors)}")
    print(f"- Semantic issues: {len(semantic_issues)}")
    
    if errors:
        print(f"\nINVALID TERMS:")
        for uberon_id, error in errors:
            print(f"  {uberon_id}: {error}")
    
    if semantic_issues:
        print(f"\nSEMANTIC ISSUES REQUIRING CORRECTION:")
        for i, issue in enumerate(semantic_issues, 1):
            print(f"  {i}. {issue}")
    
    # Save comprehensive validation results
    results = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'total_terms': len(processed_ids),
        'valid_terms': validated,
        'invalid_terms': len(errors),
        'semantic_issues': len(semantic_issues),
        'errors': [{'uberon_id': uid, 'error': err} for uid, err in errors],
        'semantic_problems': semantic_issues
    }
    
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/tissue_uberon_validation_results.json', 'w') as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    main()