#!/usr/bin/env python3
"""
Test the ENSG extraction logic for complex ENCODE IDs.
"""

def test_ensg_extraction():
    """Test the ENSG extraction logic."""
    
    test_cases = [
        "ENSG00000000003.14;ENSG00000000003.10",
        "ENSG00000000419.8;ENSG00000000419.12", 
        "ENSG00000000457.13;ENSG00000000457.9",
        "ENSG00000001234.5",
        "12345",
        "gSpikein_ERCC_123"
    ]
    
    for original_id_source in test_cases:
        print(f"\nTesting: {original_id_source}")
        
        # Simulate the extraction logic
        target_id = ''
        mapping_source = 'unmapped'
        
        # Handle complex versioned Ensembl IDs
        if not target_id.startswith(('ENSG', 'ENTREZ:', 'gSpikein')) and 'ENSG' in original_id_source:
            # Try to extract clean ENSG ID from complex format
            ensg_parts = []
            if ';' in original_id_source:
                # Split by semicolon and extract ENSG IDs
                parts = original_id_source.split(';')
                for part in parts:
                    if 'ENSG' in part:
                        clean_ensg = part.split('.')[0].strip()  # Remove version
                        if clean_ensg.startswith('ENSG'):
                            ensg_parts.append(clean_ensg)
            elif 'ENSG' in original_id_source:
                # Single ENSG, possibly with version
                clean_ensg = original_id_source.split('.')[0].strip()
                if clean_ensg.startswith('ENSG'):
                    ensg_parts.append(clean_ensg)
            
            # Use the first valid ENSG ID found
            if ensg_parts:
                candidate_ensg = ensg_parts[0]
                target_id = candidate_ensg
                mapping_source = 'ensg_extraction'
            else:
                target_id = "fallback"
                mapping_source = 'unmapped'
        
        # Handle other cases
        elif original_id_source.startswith('gSpikein'):
            target_id = original_id_source
            mapping_source = 'spike_in'
        elif original_id_source.isdigit():
            target_id = f"ENTREZ:{original_id_source}"
            mapping_source = 'entrez_id_fallback'
        else:
            target_id = original_id_source
            mapping_source = 'unmapped'
        
        print(f"  Result: {target_id} (source: {mapping_source})")

if __name__ == '__main__':
    test_ensg_extraction()