#!/usr/bin/env python3
"""
Investigate the partner's specific finding about UBERON:0016529 for visceral adipose tissue
"""

import json
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
            return {
                'valid': False,
                'uberon_id': uberon_id,
                'error': f'HTTP {response.status_code}'
            }
    except Exception as e:
        return {
            'valid': False,
            'uberon_id': uberon_id,
            'error': str(e)
        }

def main():
    print("Investigating Partner's UBERON:0016529 Issue")
    print("=" * 50)
    
    # Load current mapping
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json', 'r') as f:
        current_mappings = json.load(f)
    
    visceral_tissue = "Adipose - Visceral (Omentum)"
    current_mapping = current_mappings.get(visceral_tissue)
    
    print(f"Tissue: {visceral_tissue}")
    print(f"Current mapping: {current_mapping}")
    print()
    
    # Check what UBERON:0016529 actually is
    print("🔍 Investigating UBERON:0016529 (partner's concern):")
    uberon_0016529 = lookup_uberon_term("UBERON:0016529")
    
    if uberon_0016529['valid']:
        print(f"  Label: {uberon_0016529['label']}")
        print(f"  Description: {uberon_0016529['description']}")
        if uberon_0016529['synonyms']:
            print(f"  Synonyms: {', '.join(uberon_0016529['synonyms'])}")
    else:
        print(f"  ❌ Invalid: {uberon_0016529['error']}")
    
    print()
    
    # Check our current mapping UBERON:0014454
    print("✅ Checking our current mapping UBERON:0014454:")
    uberon_0014454 = lookup_uberon_term("UBERON:0014454")
    
    if uberon_0014454['valid']:
        print(f"  Label: {uberon_0014454['label']}")
        print(f"  Description: {uberon_0014454['description']}")
        if uberon_0014454['synonyms']:
            print(f"  Synonyms: {', '.join(uberon_0014454['synonyms'])}")
    else:
        print(f"  ❌ Invalid: {uberon_0014454['error']}")
    
    print()
    
    # Also check a few related adipose tissue terms
    related_terms = [
        "UBERON:0001013",  # adipose tissue
        "UBERON:0002190",  # subcutaneous adipose tissue
        "UBERON:0003427",  # abdominal adipose tissue
    ]
    
    print("📚 Related adipose tissue terms for context:")
    for term in related_terms:
        result = lookup_uberon_term(term)
        if result['valid']:
            print(f"  {term}: {result['label']}")
            if 'visceral' in result['description'].lower() or 'omental' in result['description'].lower():
                print(f"    🎯 Contains visceral/omental: {result['description'][:100]}...")
        else:
            print(f"  {term}: ❌ {result['error']}")
    
    print()
    
    # Analysis
    print("📊 ANALYSIS:")
    if uberon_0016529['valid'] and uberon_0014454['valid']:
        print(f"Partner concern: UBERON:0016529 = '{uberon_0016529['label']}'")
        print(f"Our current:     UBERON:0014454 = '{uberon_0014454['label']}'")
        
        # Check if 0016529 is appropriate for visceral adipose tissue
        desc_0016529 = uberon_0016529['description'].lower()
        label_0016529 = uberon_0016529['label'].lower()
        
        if any(term in desc_0016529 or term in label_0016529 for term in ['visceral', 'omental', 'intra-abdominal', 'peritoneal']):
            print("⚠️  UBERON:0016529 appears related to visceral/omental anatomy")
        else:
            print("❌ UBERON:0016529 does NOT appear suitable for visceral adipose tissue")
            
        desc_0014454 = uberon_0014454['description'].lower()
        label_0014454 = uberon_0014454['label'].lower()
        
        if any(term in desc_0014454 or term in label_0014454 for term in ['visceral', 'omental', 'intra-abdominal', 'peritoneal']):
            print("✅ UBERON:0014454 IS appropriate for visceral adipose tissue")
        else:
            print("⚠️  UBERON:0014454 relationship to visceral adipose needs verification")
    
    # Save investigation results
    results = {
        'investigation_date': '2025-08-14',
        'partner_concern': 'UBERON:0016529 incorrectly used for visceral adipose tissue',
        'tissue_name': visceral_tissue,
        'current_mapping': current_mapping,
        'uberon_0016529': uberon_0016529,
        'uberon_0014454': uberon_0014454,
        'related_terms': {term: lookup_uberon_term(term) for term in related_terms}
    }
    
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/debug_uberon/partner_issue_investigation.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Investigation results saved to reports/debug_uberon/partner_issue_investigation.json")

if __name__ == "__main__":
    main()