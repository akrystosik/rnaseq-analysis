#!/usr/bin/env python3
"""
Semantic validation of tissue-to-UBERON mappings
Checks if the UBERON term actually matches the tissue description
"""

import json
import requests
import time
from urllib.parse import quote

def get_uberon_details(uberon_id):
    """Get UBERON term details via OLS4 API"""
    if uberon_id.startswith('CL:'):
        return {'label': 'Cell Ontology term', 'description': 'Cell type', 'valid': True}
    
    if not uberon_id.startswith('UBERON:'):
        return {'valid': False, 'error': f'Invalid format: {uberon_id}'}
    
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
        else:
            return {'valid': False, 'error': f'API error {response.status_code}'}
    except Exception as e:
        return {'valid': False, 'error': f'Request failed: {str(e)}'}

def check_semantic_match(tissue_name, uberon_details):
    """Check if tissue name semantically matches UBERON term"""
    if not uberon_details.get('valid'):
        return False, uberon_details.get('error', 'Invalid term')
    
    tissue_lower = tissue_name.lower()
    label_lower = uberon_details['label'].lower()
    desc_lower = uberon_details['description'].lower()
    
    # Known problematic cases
    semantic_issues = []
    
    # Skin tissue checks
    if 'skin' in tissue_lower:
        if 'placenta' in label_lower:
            semantic_issues.append(f"Skin tissue '{tissue_name}' mapped to placenta term '{uberon_details['label']}'")
        elif 'suprapubic' in tissue_lower and 'suprapubic' not in label_lower:
            semantic_issues.append(f"Suprapubic skin '{tissue_name}' not mapped to suprapubic term '{uberon_details['label']}'")
        elif 'lower leg' in tissue_lower and 'lower leg' not in label_lower:
            semantic_issues.append(f"Lower leg skin '{tissue_name}' not mapped to lower leg term '{uberon_details['label']}'")
    
    # Adipose tissue checks
    if 'adipose' in tissue_lower:
        if 'adipose' not in label_lower and 'fat' not in label_lower:
            semantic_issues.append(f"Adipose tissue '{tissue_name}' not mapped to adipose/fat term '{uberon_details['label']}'")
        if 'visceral' in tissue_lower and 'visceral' not in label_lower and 'abdominal' not in label_lower:
            semantic_issues.append(f"Visceral adipose '{tissue_name}' not mapped to visceral/abdominal term '{uberon_details['label']}'")
    
    # Brain tissue checks
    if 'brain' in tissue_lower:
        if 'brain' not in label_lower and 'cerebral' not in label_lower and 'cortex' not in label_lower and 'hippocampus' not in label_lower and 'amygdala' not in label_lower and 'caudate' not in label_lower and 'putamen' not in label_lower and 'hypothalamus' not in label_lower and 'cerebellum' not in label_lower and 'substantia nigra' not in label_lower and 'nucleus accumbens' not in label_lower and 'spinal cord' not in label_lower:
            semantic_issues.append(f"Brain tissue '{tissue_name}' not mapped to brain-related term '{uberon_details['label']}'")
    
    # General organ checks
    organs = ['heart', 'lung', 'liver', 'kidney', 'pancreas', 'stomach', 'colon', 'prostate', 'ovary', 'testis', 'thyroid', 'spleen']
    for organ in organs:
        if organ in tissue_lower and organ not in label_lower:
            semantic_issues.append(f"{organ.title()} tissue '{tissue_name}' not mapped to {organ} term '{uberon_details['label']}'")
    
    return len(semantic_issues) == 0, semantic_issues

def main():
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json', 'r') as f:
        mappings = json.load(f)
    
    print("Semantic Validation of Tissue-to-UBERON Mappings")
    print("=" * 60)
    
    all_issues = []
    processed = set()  # Track processed UBERON IDs to avoid duplicates
    
    for tissue_name, uberon_id in mappings.items():
        if uberon_id in processed:
            continue
        processed.add(uberon_id)
        
        print(f"\\nValidating {uberon_id} for tissue: {tissue_name}")
        uberon_details = get_uberon_details(uberon_id)
        
        if uberon_details.get('valid'):
            print(f"  → {uberon_details['label']}: {uberon_details['description'][:100]}{'...' if len(uberon_details['description']) > 100 else ''}")
            
            is_match, issues = check_semantic_match(tissue_name, uberon_details)
            if not is_match:
                all_issues.extend(issues)
                for issue in issues:
                    print(f"  ⚠️  SEMANTIC ISSUE: {issue}")
        else:
            print(f"  ❌ {uberon_details.get('error', 'Unknown error')}")
        
        time.sleep(0.1)
    
    print(f"\\n{'=' * 60}")
    print(f"SUMMARY")
    print(f"Total unique UBERON IDs checked: {len(processed)}")
    print(f"Semantic issues found: {len(all_issues)}")
    
    if all_issues:
        print(f"\\n🔍 ALL SEMANTIC ISSUES:")
        for i, issue in enumerate(all_issues, 1):
            print(f"{i}. {issue}")
            
        # Find affected tissues for each problematic UBERON ID
        print(f"\\n📋 TISSUES AFFECTED BY EACH ISSUE:")
        for issue in all_issues:
            print(f"\\n{issue}")
            # Find all tissues that use problematic mappings mentioned in the issue
            for tissue_name, uberon_id in mappings.items():
                uberon_details = get_uberon_details(uberon_id)
                if uberon_details.get('valid'):
                    if (uberon_details['label'].lower() in issue.lower() or 
                        uberon_id in issue):
                        print(f"  → {tissue_name}")

if __name__ == "__main__":
    main()