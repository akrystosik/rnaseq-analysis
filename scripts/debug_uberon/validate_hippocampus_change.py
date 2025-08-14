#!/usr/bin/env python3
"""
Validate the hippocampus mapping change from UBERON:0002310 to UBERON:0002421
"""

import requests
from urllib.parse import quote

def lookup_uberon_term(uberon_id):
    """Look up a specific UBERON term"""
    try:
        iri = f"http://purl.obolibrary.org/obo/{uberon_id.replace(':', '_')}"
        double_encoded_iri = quote(quote(iri, safe=''), safe='')
        url = f"https://www.ebi.ac.uk/ols4/api/ontologies/uberon/terms/{double_encoded_iri}"
        
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return {
                'valid': True,
                'uberon_id': uberon_id,
                'label': data.get('label', 'Unknown'),
                'description': data.get('description', [''])[0] if data.get('description') else '',
                'synonyms': data.get('synonyms', [])
            }
        else:
            return {'valid': False, 'uberon_id': uberon_id, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'valid': False, 'uberon_id': uberon_id, 'error': str(e)}

def main():
    print("Validating Hippocampus Mapping Change")
    print("=" * 40)
    
    # Check the old mapping
    print("❌ Previous mapping UBERON:0002310:")
    old_term = lookup_uberon_term("UBERON:0002310")
    if old_term['valid']:
        print(f"  Label: {old_term['label']}")
        print(f"  Description: {old_term['description']}")
    else:
        print(f"  Error: {old_term['error']}")
    
    print()
    
    # Check the new mapping
    print("✅ New mapping UBERON:0002421:")
    new_term = lookup_uberon_term("UBERON:0002421")
    if new_term['valid']:
        print(f"  Label: {new_term['label']}")
        print(f"  Description: {new_term['description']}")
        if new_term['synonyms']:
            print(f"  Synonyms: {', '.join(new_term['synonyms'])}")
    else:
        print(f"  Error: {new_term['error']}")
    
    print()
    
    # Analysis
    if old_term['valid'] and new_term['valid']:
        print("📊 ANALYSIS:")
        print(f"Old: {old_term['label']} - {old_term['description'][:100]}...")
        print(f"New: {new_term['label']} - {new_term['description'][:100]}...")
        
        # Check which is more appropriate for general "hippocampus" tissue
        old_is_substructure = any(term in old_term['label'].lower() for term in ['fimbria', 'alveus', 'field'])
        new_is_general = 'hippocampus' in new_term['label'].lower() and not any(term in new_term['label'].lower() for term in ['fimbria', 'alveus', 'field'])
        
        if old_is_substructure and new_is_general:
            print("✅ IMPROVEMENT: Changed from specific substructure to general hippocampus")
        elif new_is_general:
            print("✅ APPROPRIATE: New mapping is general hippocampus")
        else:
            print("❓ NEEDS REVIEW: Both terms may be too specific or inappropriate")

if __name__ == "__main__":
    main()