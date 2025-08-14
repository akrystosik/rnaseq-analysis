#!/usr/bin/env python3
"""
Systematic search for proper UBERON IDs using OLS4 API
Compare results with current mappings to identify discrepancies
"""

import json
import requests
import time
from urllib.parse import quote
import re

def clean_tissue_text(tissue_text):
    """Clean tissue text for better search matching"""
    clean_text = tissue_text.lower()
    
    # Remove common prefixes
    clean_text = re.sub(r'^(brain - |cells - |artery - |colon - |esophagus - |heart - |kidney - |skin - |cervix - )', '', clean_text)
    
    # Remove parenthetical info but keep important context
    # Keep basal ganglia context but remove BA codes
    if 'basal ganglia' in clean_text:
        clean_text = re.sub(r'\s*\(ba\d+\)', '', clean_text)  # Remove BA24 etc
    else:
        clean_text = re.sub(r'\s*\([^)]*\)', '', clean_text)  # Remove all parentheticals
    
    clean_text = clean_text.strip()
    return clean_text

def generate_search_variants(tissue_text):
    """Generate search term variants for better matching"""
    variants = []
    clean_text = clean_tissue_text(tissue_text)
    variants.append(clean_text)
    
    # Specific mappings for known problematic terms
    if 'visceral' in clean_text and 'adipose' in clean_text:
        variants.extend([
            'intra-abdominal adipose tissue',
            'visceral adipose tissue', 
            'omental adipose tissue',
            'mesenteric fat'
        ])
    elif 'subcutaneous' in clean_text and 'adipose' in clean_text:
        variants.extend(['subcutaneous adipose tissue', 'subcutaneous fat'])
    elif 'caudate' in clean_text:
        variants.extend(['caudate nucleus', 'striatum'])
    elif 'frontal cortex' in clean_text:
        variants.extend(['frontal cortex', 'prefrontal cortex', 'brodmann area 9'])
    elif 'anterior cingulate cortex' in clean_text:
        variants.extend(['anterior cingulate cortex', 'brodmann area 24'])
    elif 'cortex' in clean_text and 'brain' not in tissue_text.lower():
        variants.extend(['cerebral cortex'])
    elif 'mammary' in clean_text:
        variants.extend(['mammary gland'])
    elif 'tibial' in clean_text and 'artery' in tissue_text.lower():
        variants.extend(['posterior tibial artery'])
    elif 'tibial' in clean_text and 'nerve' in tissue_text.lower():
        variants.extend(['tibial nerve'])
    
    return variants[:4]  # Limit to avoid too many API calls

def search_uberon_for_tissue(tissue_text, max_results=5):
    """Search UBERON ontology for matching terms"""
    search_variants = generate_search_variants(tissue_text)
    all_results = []
    
    for variant in search_variants:
        try:
            url = "https://www.ebi.ac.uk/ols4/api/search"
            params = {
                'q': variant,
                'ontology': 'uberon',
                'rows': max_results,
                'start': 0,
                'format': 'json'
            }
            
            print(f"    Searching: '{variant}'")
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
                                'search_variant': variant,
                                'score': doc.get('score', 0)
                            }
                            all_results.append(result)
            
            time.sleep(0.3)  # Be respectful to API
            
        except Exception as e:
            print(f"      Error searching '{variant}': {e}")
    
    # Remove duplicates and sort by score
    unique_results = {}
    for result in all_results:
        uberon_id = result['uberon_id']
        if uberon_id not in unique_results or result['score'] > unique_results[uberon_id]['score']:
            unique_results[uberon_id] = result
    
    return sorted(unique_results.values(), key=lambda x: x['score'], reverse=True)

def validate_current_mapping(uberon_id):
    """Validate current UBERON ID by looking it up"""
    if uberon_id.startswith('CL:'):
        return {'valid': True, 'label': 'Cell Ontology term', 'note': 'Skipping CL validation'}
    
    try:
        iri = f"http://purl.obolibrary.org/obo/{uberon_id.replace(':', '_')}"
        double_encoded_iri = quote(quote(iri, safe=''), safe='')
        url = f"https://www.ebi.ac.uk/ols4/api/ontologies/uberon/terms/{double_encoded_iri}"
        
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return {
                'valid': True,
                'label': data.get('label', 'Unknown'),
                'description': data.get('description', [''])[0] if data.get('description') else ''
            }
        else:
            return {'valid': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'valid': False, 'error': str(e)}

def main():
    print("Systematic UBERON Search and Validation")
    print("=" * 60)
    
    # Load current mappings
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json', 'r') as f:
        current_mappings = json.load(f)
    
    # Get unique tissue names (avoiding duplicate case variations)
    unique_tissues = {}
    for tissue_name in current_mappings.keys():
        normalized = tissue_name.lower()
        if normalized not in unique_tissues:
            unique_tissues[normalized] = tissue_name
    
    tissues_to_check = list(unique_tissues.values())
    print(f"Checking {len(tissues_to_check)} unique tissue types\n")
    
    results = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'total_tissues': len(tissues_to_check),
        'discrepancies': [],
        'validated_mappings': [],
        'search_failures': []
    }
    
    for i, tissue_name in enumerate(tissues_to_check, 1):
        print(f"[{i:2d}/{len(tissues_to_check)}] {tissue_name}")
        current_uberon = current_mappings[tissue_name]
        print(f"  Current: {current_uberon}")
        
        # Validate current mapping
        current_validation = validate_current_mapping(current_uberon)
        if current_validation['valid']:
            print(f"    ✓ Valid: {current_validation['label']}")
        else:
            print(f"    ✗ Invalid: {current_validation['error']}")
        
        # Search for alternatives
        search_results = search_uberon_for_tissue(tissue_name)
        
        if search_results:
            best_match = search_results[0]
            print(f"  Best match: {best_match['uberon_id']} - {best_match['label']}")
            
            # Check for discrepancy
            if current_uberon != best_match['uberon_id']:
                print(f"    ⚠️  DISCREPANCY FOUND!")
                
                discrepancy = {
                    'tissue_name': tissue_name,
                    'current_mapping': {
                        'uberon_id': current_uberon,
                        'validation': current_validation
                    },
                    'suggested_mapping': best_match,
                    'all_search_results': search_results[:3]
                }
                results['discrepancies'].append(discrepancy)
            else:
                print(f"    ✓ Mapping confirmed")
                results['validated_mappings'].append({
                    'tissue_name': tissue_name,
                    'uberon_id': current_uberon,
                    'validation': current_validation
                })
        else:
            print(f"    ❌ No search results found")
            results['search_failures'].append({
                'tissue_name': tissue_name,
                'current_uberon': current_uberon
            })
        
        print()
        time.sleep(0.5)  # Be respectful to API
        
        # Limit for initial run to avoid overwhelming API
        if i >= 30:
            print("Limiting to first 30 tissues for initial analysis...")
            break
    
    # Summary
    print("=" * 60)
    print("ANALYSIS SUMMARY")
    print(f"Tissues analyzed: {len(results['validated_mappings']) + len(results['discrepancies']) + len(results['search_failures'])}")
    print(f"Confirmed mappings: {len(results['validated_mappings'])}")
    print(f"Discrepancies found: {len(results['discrepancies'])}")
    print(f"Search failures: {len(results['search_failures'])}")
    
    if results['discrepancies']:
        print(f"\nDISCREPANCIES REQUIRING REVIEW:")
        for disc in results['discrepancies']:
            print(f"\n• {disc['tissue_name']}")
            print(f"  Current:   {disc['current_mapping']['uberon_id']}")
            if disc['current_mapping']['validation']['valid']:
                print(f"    Label: {disc['current_mapping']['validation']['label']}")
            print(f"  Suggested: {disc['suggested_mapping']['uberon_id']}")
            print(f"    Label: {disc['suggested_mapping']['label']}")
    
    # Save detailed results
    output_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/debug_uberon/systematic_search_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Detailed results saved to: {output_file}")

if __name__ == "__main__":
    main()