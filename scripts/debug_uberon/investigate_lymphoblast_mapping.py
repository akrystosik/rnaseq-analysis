#!/usr/bin/env python3
"""
Investigate Lymphoblast mapping - CL:0000542
Check if this correctly represents lymphoblasts vs lymphocytes
And verify developmental stage appropriateness
"""

import requests
import json
from urllib.parse import quote

def lookup_cl_term(cl_id):
    """Look up a specific Cell Ontology term"""
    try:
        iri = f"http://purl.obolibrary.org/obo/{cl_id.replace(':', '_')}"
        double_encoded_iri = quote(quote(iri, safe=''), safe='')
        url = f"https://www.ebi.ac.uk/ols4/api/ontologies/cl/terms/{double_encoded_iri}"
        
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return {
                'valid': True,
                'cl_id': cl_id,
                'label': data.get('label', 'Unknown'),
                'description': data.get('description', [''])[0] if data.get('description') else '',
                'synonyms': data.get('synonyms', [])
            }
        else:
            return {'valid': False, 'cl_id': cl_id, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'valid': False, 'cl_id': cl_id, 'error': str(e)}

def search_lymphocyte_terms():
    """Search for lymphoblast and lymphocyte related terms"""
    search_terms = [
        'lymphoblast',
        'lymphocyte', 
        'immature lymphocyte',
        'mature lymphocyte',
        'lymphocyte precursor'
    ]
    
    results = []
    for term in search_terms:
        try:
            url = "https://www.ebi.ac.uk/ols4/api/search"
            params = {
                'q': term,
                'ontology': 'cl',
                'rows': 8,
                'format': 'json'
            }
            
            print(f"Searching Cell Ontology: '{term}'")
            response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                docs = data.get('response', {}).get('docs', [])
                
                for doc in docs:
                    if doc.get('ontology_name') == 'cl':
                        iri = doc.get('iri', '')
                        if 'CL_' in iri:
                            cl_id = iri.split('CL_')[-1]
                            cl_id = f"CL:{cl_id}"
                            
                            result = {
                                'cl_id': cl_id,
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
    print("Investigating Lymphoblast Mapping")
    print("=" * 35)
    
    # Check current mapping
    current_cl = "CL:0000542"
    print(f"Current mapping: {current_cl}")
    
    current_term = lookup_cl_term(current_cl)
    if current_term['valid']:
        print(f"  Label: {current_term['label']}")
        print(f"  Description: {current_term['description']}")
        if current_term['synonyms']:
            print(f"  Synonyms: {', '.join(current_term['synonyms'][:5])}...")
    else:
        print(f"  Error: {current_term['error']}")
    
    print()
    
    # Analyze if this is actually lymphoblast or lymphocyte
    if current_term['valid']:
        label_lower = current_term['label'].lower()
        desc_lower = current_term['description'].lower()
        
        is_lymphoblast = 'lymphoblast' in label_lower or 'lymphoblast' in desc_lower
        is_lymphocyte = 'lymphocyte' in label_lower or 'lymphocyte' in desc_lower
        is_immature = any(term in label_lower or term in desc_lower 
                         for term in ['immature', 'precursor', 'progenitor', 'blast'])
        is_mature = 'mature' in label_lower or 'mature' in desc_lower
        
        print("🔍 DEVELOPMENTAL STAGE ANALYSIS:")
        print(f"  Current term label: '{current_term['label']}'")
        print(f"  Contains 'lymphoblast': {is_lymphoblast}")
        print(f"  Contains 'lymphocyte': {is_lymphocyte}")
        print(f"  Indicates immature/precursor: {is_immature}")
        print(f"  Indicates mature: {is_mature}")
        
        if is_lymphoblast and not is_mature:
            print("  ✅ Current term correctly represents lymphoblast (immature)")
        elif is_lymphocyte and is_mature:
            print("  ⚠️  Current term represents mature lymphocyte, not lymphoblast")
        elif is_lymphocyte and not is_mature:
            print("  ❓ Current term represents lymphocyte but maturity unclear")
        else:
            print("  ❓ Current term classification unclear")
    
    print()
    
    # Search for alternatives and compare
    print("🔎 SEARCHING FOR LYMPHOID CELL TYPES:")
    search_results = search_lymphocyte_terms()
    
    # Remove duplicates
    unique_results = {}
    for result in search_results:
        cl_id = result['cl_id']
        if cl_id not in unique_results:
            unique_results[cl_id] = result
    
    # Sort by relevance and group by type
    sorted_results = sorted(unique_results.values(), key=lambda x: x['search_term'])
    
    lymphoblast_terms = []
    lymphocyte_terms = []
    other_terms = []
    
    for result in sorted_results:
        label_lower = result['label'].lower()
        if 'lymphoblast' in label_lower:
            lymphoblast_terms.append(result)
        elif 'lymphocyte' in label_lower:
            lymphocyte_terms.append(result)
        else:
            other_terms.append(result)
    
    print(f"\n📋 LYMPHOBLAST TERMS FOUND:")
    for i, result in enumerate(lymphoblast_terms[:5], 1):
        current_indicator = "👉 CURRENT" if result['cl_id'] == current_cl else ""
        print(f"{i}. {result['cl_id']}: {result['label']} {current_indicator}")
        print(f"   Description: {result['description'][:80]}...")
        
        # Check if this is appropriate for tissue samples
        desc_lower = result['description'].lower()
        if 'immature' in desc_lower or 'precursor' in desc_lower:
            print(f"   ✅ Confirmed immature/precursor cell")
        elif 'mature' in desc_lower:
            print(f"   ⚠️  Mature cell (unexpected for lymphoblast)")
    
    print(f"\n📋 LYMPHOCYTE TERMS FOUND:")
    for i, result in enumerate(lymphocyte_terms[:5], 1):
        print(f"{i}. {result['cl_id']}: {result['label']}")
        print(f"   Description: {result['description'][:80]}...")
        
        desc_lower = result['description'].lower()
        label_lower = result['label'].lower()
        if 'mature' in desc_lower or 'mature' in label_lower:
            print(f"   ✅ Confirmed mature lymphocyte")
        elif 'immature' in desc_lower:
            print(f"   ⚠️  Immature lymphocyte")
    
    print(f"\n📊 ASSESSMENT:")
    print("For tissue samples labeled 'Lymphoblast':")
    print("• Should represent immature lymphoid precursor cells")
    print("• Should NOT represent mature, functional lymphocytes") 
    print("• Term should explicitly indicate developmental stage")
    
    if current_term['valid']:
        current_label = current_term['label'].lower()
        current_desc = current_term['description'].lower()
        
        if 'lymphoblast' in current_label:
            print("✅ Current mapping label correctly indicates lymphoblast")
        elif 'lymphocyte' in current_label and 'mature' not in current_desc:
            print("❓ Current mapping uses lymphocyte terminology - check appropriateness")
        elif 'lymphocyte' in current_label and 'mature' in current_desc:
            print("❌ Current mapping represents mature lymphocyte, not lymphoblast")
        
    print(f"\n💡 CONTEXT:")
    print("In RNA-seq experiments:")
    print("• Lymphoblast samples often come from cell lines or cultured precursor cells")
    print("• These may be EBV-transformed lymphoblastoid cell lines")
    print("• Different from primary mature lymphocytes circulating in blood")

if __name__ == "__main__":
    main()