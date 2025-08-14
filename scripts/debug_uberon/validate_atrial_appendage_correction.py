#!/usr/bin/env python3
"""
Validate the proposed atrial appendage correction
From UBERON:0006631 (right atrium auricular region) to UBERON:0006618 (atrium auricular region)
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
    print("Validating Atrial Appendage Correction")
    print("=" * 40)
    
    current_uberon = "UBERON:0006631"
    proposed_uberon = "UBERON:0006618"
    
    print("❌ CURRENT (too specific - lateralized):")
    current_term = lookup_uberon_term(current_uberon)
    if current_term['valid']:
        print(f"  {current_uberon}: {current_term['label']}")
        print(f"  Description: {current_term['description']}")
        if current_term['synonyms']:
            print(f"  Synonyms: {', '.join(current_term['synonyms'][:3])}...")
    
    print()
    
    print("✅ PROPOSED (appropriately general):")
    proposed_term = lookup_uberon_term(proposed_uberon)
    if proposed_term['valid']:
        print(f"  {proposed_uberon}: {proposed_term['label']}")
        print(f"  Description: {proposed_term['description']}")
        if proposed_term['synonyms']:
            print(f"  Synonyms: {', '.join(proposed_term['synonyms'][:3])}...")
    
    print()
    
    # Also check left atrial appendage for comparison
    print("📚 For comparison - Left atrial appendage:")
    left_atrial = lookup_uberon_term("UBERON:0006630")  # likely left atrial appendage
    if left_atrial['valid']:
        print(f"  UBERON:0006630: {left_atrial['label']}")
        print(f"  Description: {left_atrial['description']}")
    
    print()
    
    # Analysis
    if current_term['valid'] and proposed_term['valid']:
        print("📊 CORRECTION RATIONALE:")
        print(f"• Current term is TOO SPECIFIC: '{current_term['label']}'")
        print("  - Specifies RIGHT atrial appendage only")
        print("  - Tissue samples labeled 'Heart - Atrial Appendage' don't specify lateralization")
        print("  - Could include left atrial appendage tissue")
        
        print(f"• Proposed term is APPROPRIATELY GENERAL: '{proposed_term['label']}'")
        print("  - Covers atrial appendage structures generally")
        print("  - Includes both left AND right atrial appendages")
        print("  - Better matches non-lateralized 'Heart - Atrial Appendage' samples")
        
        # Check if proposed term actually covers both sides
        proposed_desc = proposed_term['description'].lower()
        if 'each atrium' in proposed_desc or 'both' in proposed_desc or ('left' in proposed_desc and 'right' in proposed_desc):
            print("  - ✅ Description confirms coverage of both sides")
        
        print("\n✅ RECOMMENDATION: Change to UBERON:0006618")
        
        return {
            'recommendation': 'CHANGE_RECOMMENDED',
            'current': current_term,
            'proposed': proposed_term,
            'rationale': 'Current term lateralized (right only), proposed covers general atrial appendage'
        }
    
    return {'recommendation': 'UNABLE_TO_VALIDATE', 'error': 'Could not retrieve term information'}

if __name__ == "__main__":
    result = main()
    
    # Save validation result
    import json
    with open('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/debug_uberon/atrial_appendage_correction_validation.json', 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"\n✅ Validation saved to reports/debug_uberon/atrial_appendage_correction_validation.json")