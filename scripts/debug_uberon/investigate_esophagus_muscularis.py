#!/usr/bin/env python3
"""
Investigate esophagus muscularis mapping - currently too specific
Need to find appropriate general esophageal muscle layer term
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

def search_esophageal_terms():
    """Search for appropriate esophageal muscle terms"""
    search_terms = [
        'esophageal muscle',
        'esophagus muscle', 
        'muscle layer of esophagus',
        'esophageal musculature'
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
    print("Investigating Esophagus Muscularis Mapping")
    print("=" * 45)
    
    # Check current mapping
    current_uberon = "UBERON:0004648"
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
    
    # Check if it's too specific
    if current_term['valid']:
        label_lower = current_term['label'].lower()
        desc_lower = current_term['description'].lower()
        
        is_specific = any(term in label_lower or term in desc_lower 
                         for term in ['mucosa', 'propria', 'externa', 'inner', 'outer', 'circular', 'longitudinal'])
        
        print("🔍 SPECIFICITY ANALYSIS:")
        if is_specific:
            print("  ⚠️  Current term appears to be layer-specific")
        else:
            print("  ✅ Current term appears appropriately general")
        
        print(f"  Label analysis: '{label_lower}'")
        print(f"  Description contains specific layer terms: {is_specific}")
    
    print()
    
    # Search for better alternatives
    print("🔎 SEARCHING FOR ALTERNATIVES:")
    search_results = search_esophageal_terms()
    
    # Remove duplicates and current term
    unique_results = {}
    for result in search_results:
        uberon_id = result['uberon_id']
        if uberon_id != current_uberon and uberon_id not in unique_results:
            unique_results[uberon_id] = result
    
    # Sort by relevance
    sorted_results = sorted(unique_results.values(), key=lambda x: x['score'], reverse=True)
    
    print(f"Found {len(sorted_results)} alternative terms:")
    for i, result in enumerate(sorted_results[:5], 1):
        print(f"\n{i}. {result['uberon_id']}: {result['label']}")
        print(f"   Description: {result['description'][:100]}...")
        print(f"   Search term: '{result['search_term']}' (score: {result['score']:.2f})")
        
        # Check if this is more general
        label_lower = result['label'].lower()
        desc_lower = result['description'].lower()
        is_general = 'muscle' in label_lower and 'esophag' in (label_lower + desc_lower)
        is_not_specific = not any(term in label_lower or term in desc_lower 
                                 for term in ['mucosa', 'propria', 'externa', 'inner', 'outer', 'circular', 'longitudinal'])
        
        if is_general and is_not_specific:
            print(f"   🎯 This appears more appropriately general!")
    
    # Also check some specific UBERON IDs that might be relevant
    print(f"\n🔍 CHECKING RELATED ESOPHAGEAL TERMS:")
    related_terms = [
        "UBERON:0001043",  # esophagus
        "UBERON:0006562",  # pharyngoesophageal muscle
        "UBERON:0004648",  # current mapping
    ]
    
    for uberon_id in related_terms:
        term_info = lookup_uberon_term(uberon_id)
        if term_info['valid']:
            print(f"\n{uberon_id}: {term_info['label']}")
            print(f"  Description: {term_info['description'][:100]}...")

if __name__ == "__main__":
    main()