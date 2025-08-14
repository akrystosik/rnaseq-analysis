#!/usr/bin/env python3
"""
Validate the proposed esophagus muscularis correction
From UBERON:0004648 (muscularis mucosa) to UBERON:0003832 (esophagus muscle)
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
    print("Validating Esophagus Muscularis Correction")
    print("=" * 45)
    
    current_uberon = "UBERON:0004648"
    proposed_uberon = "UBERON:0003832"
    
    print("❌ CURRENT (too specific):")
    current_term = lookup_uberon_term(current_uberon)
    if current_term['valid']:
        print(f"  {current_uberon}: {current_term['label']}")
        print(f"  Description: {current_term['description']}")
        if current_term['synonyms']:
            print(f"  Synonyms: {', '.join(current_term['synonyms'])}")
    
    print()
    
    print("✅ PROPOSED (appropriately general):")
    proposed_term = lookup_uberon_term(proposed_uberon)
    if proposed_term['valid']:
        print(f"  {proposed_uberon}: {proposed_term['label']}")
        print(f"  Description: {proposed_term['description']}")
        if proposed_term['synonyms']:
            print(f"  Synonyms: {', '.join(proposed_term['synonyms'])}")
    
    print()
    
    # Analysis
    if current_term['valid'] and proposed_term['valid']:
        print("📊 CORRECTION RATIONALE:")
        print(f"• Current term is TOO SPECIFIC: '{current_term['label']}'")
        print("  - Refers only to muscularis mucosa (thin layer within mucosa)")
        print("  - Does not represent general esophageal muscle tissue")
        
        print(f"• Proposed term is APPROPRIATELY GENERAL: '{proposed_term['label']}'")
        print("  - Covers any muscle organ that is part of esophagus") 
        print("  - Includes both muscularis mucosa AND muscularis propria layers")
        print("  - Better matches 'Esophagus - Muscularis' tissue samples")
        
        print("\n✅ RECOMMENDATION: Change to UBERON:0003832")
        
        return {
            'recommendation': 'CHANGE_RECOMMENDED',
            'current': current_term,
            'proposed': proposed_term,
            'rationale': 'Current term too specific (only muscularis mucosa), proposed covers general esophageal muscle tissue'
        }
    
    return {'recommendation': 'UNABLE_TO_VALIDATE', 'error': 'Could not retrieve term information'}

if __name__ == "__main__":
    result = main()
    
    # Save validation result
    import json
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/debug_uberon/esophagus_correction_validation.json', 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"\n✅ Validation saved to reports/debug_uberon/esophagus_correction_validation.json")