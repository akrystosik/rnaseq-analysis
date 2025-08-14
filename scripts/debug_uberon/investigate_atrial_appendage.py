#!/usr/bin/env python3
"""
Investigate Heart - Atrial Appendage mapping
Check if UBERON:0006631 is appropriately specific or too narrow
"""

import requests
import json
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

def search_atrial_terms():
    """Search for atrial appendage and related cardiac terms"""
    search_terms = [
        'atrial appendage',
        'atrial auricle', 
        'right atrial appendage',
        'left atrial appendage',
        'atrium appendage'
    ]
    
    results = []
    for term in search_terms:
        try:
            url = "https://www.ebi.ac.uk/ols4/api/search"
            params = {
                'q': term,
                'ontology': 'uberon',
                'rows': 5,
                'format': 'json'
            }
            
            print(f"Searching: '{term}'")
            response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                docs = data.get('response', {}).get('docs', [])
                
                for doc in docs:
                    if doc.get('ontology_name') == 'uberon':
                        iri = doc.get('iri', '')
                        if 'UBERON_' in iri:
                            uberon_id = iri.split('UBERON_')[-1]
                            uberon_id = f"UBERON:{uberon_id}"
                            
                            result = {
                                'uberon_id': uberon_id,
                                'label': doc.get('label', ''),
                                'description': doc.get('description', [''])[0] if doc.get('description') else '',
                                'search_term': term,
                                'score': doc.get('score', 0)
                            }
                            results.append(result)
        except Exception as e:
            print(f"  Error: {e}")
    
    return results

def main():
    print("Investigating Heart - Atrial Appendage Mapping")
    print("=" * 50)
    
    # Check current mapping
    current_uberon = "UBERON:0006631"
    print(f"Current mapping: {current_uberon}")
    
    current_term = lookup_uberon_term(current_uberon)
    if current_term['valid']:
        print(f"  Label: {current_term['label']}")
        print(f"  Description: {current_term['description']}")
        if current_term['synonyms']:
            print(f"  Synonyms: {', '.join(current_term['synonyms'])}")
    else:
        print(f"  Error: {current_term['error']}")
    
    print()
    
    # Analyze specificity
    if current_term['valid']:
        label_lower = current_term['label'].lower()
        desc_lower = current_term['description'].lower()
        
        # Check if it specifies left vs right
        is_lateralized = any(side in label_lower or side in desc_lower 
                           for side in ['left', 'right', 'dextra', 'sinistra'])
        
        # Check if it's general atrial appendage
        is_general_appendage = 'appendage' in label_lower and not is_lateralized
        
        print("🔍 SPECIFICITY ANALYSIS:")
        print(f"  Label: '{current_term['label']}'")
        print(f"  Specifies left/right: {is_lateralized}")
        print(f"  General atrial appendage: {is_general_appendage}")
        
        if is_lateralized:
            print("  ⚠️  Current term is lateralized (left/right specific)")
        elif is_general_appendage:
            print("  ✅ Current term is appropriately general")
        else:
            print("  ❓ Need to assess appropriateness")
    
    print()
    
    # Search for alternatives
    print("🔎 SEARCHING FOR ALTERNATIVES:")
    search_results = search_atrial_terms()
    
    # Remove duplicates and current term
    unique_results = {}
    for result in search_results:
        uberon_id = result['uberon_id']
        if uberon_id not in unique_results:
            unique_results[uberon_id] = result
    
    # Sort by relevance
    sorted_results = sorted(unique_results.values(), key=lambda x: x['score'], reverse=True)
    
    print(f"Found {len(sorted_results)} terms:")
    for i, result in enumerate(sorted_results[:6], 1):
        current_indicator = "👉 CURRENT" if result['uberon_id'] == current_uberon else ""
        print(f"\n{i}. {result['uberon_id']}: {result['label']} {current_indicator}")
        print(f"   Description: {result['description'][:100]}...")
        print(f"   Search term: '{result['search_term']}' (score: {result['score']:.2f})")
        
        # Analyze each result
        label_lower = result['label'].lower()
        desc_lower = result['description'].lower()
        is_lateralized = any(side in label_lower or side in desc_lower 
                           for side in ['left', 'right'])
        is_general = 'appendage' in label_lower and not is_lateralized
        
        if is_lateralized:
            print(f"   ⚠️  Lateralized (left/right specific)")
        elif is_general:
            print(f"   ✅ General atrial appendage")
    
    # Check related cardiac terms
    print(f"\n🔍 RELATED CARDIAC TERMS:")
    related_terms = [
        "UBERON:0002078",  # right atrium
        "UBERON:0002079",  # left atrium  
        "UBERON:0000948",  # heart
        "UBERON:0006631",  # current mapping
    ]
    
    for uberon_id in related_terms:
        term_info = lookup_uberon_term(uberon_id)
        if term_info['valid']:
            print(f"\n{uberon_id}: {term_info['label']}")
            print(f"  Description: {term_info['description'][:100]}...")
    
    print(f"\n📊 ASSESSMENT:")
    print("For tissue samples labeled 'Heart - Atrial Appendage':")
    print("• Could represent either left OR right atrial appendage tissue")
    print("• Or could represent general atrial appendage tissue without lateralization")
    print("• Need to determine if current mapping is appropriately general")
    
    if current_term['valid']:
        current_label = current_term['label'].lower()
        if 'left' in current_label or 'right' in current_label:
            print("❌ Current mapping is lateralized - may be too specific")
        else:
            print("✅ Current mapping appears appropriately general")

if __name__ == "__main__":
    main()