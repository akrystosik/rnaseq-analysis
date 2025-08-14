# Tissue-to-UBERON Mapping Audit

## Critical Findings

### 1. Confirmed Mapping Error - Visceral Adipose Tissue

**Issue**: Lines 3 and 68 in `tissue_to_uberon.json`
```json
"Adipose - Visceral (Omentum)": "UBERON:0016529",
"adipose - visceral (omentum)": "UBERON:0016529",
```

**Problem**: UBERON:0016529 is incorrect for visceral/omental adipose tissue

**Correct Mapping**: 
- **UBERON:0014454** - intra-abdominal adipose tissue (visceral adipose tissue)
- Definition: "Adipose tissue located inside the peritoneal cavity, packed in between internal organs and torso"

**Evidence**: 
- Web search confirmed UBERON:0014454 is the proper term for intra-abdominal/visceral adipose tissue
- This term specifically covers omental adipose tissue located in the peritoneal cavity

### 2. Impact Assessment
- Affects both capitalized and lowercase versions of visceral adipose tissue mappings
- Partner using run_20250502_211754 would have incorrect ontology annotations
- All downstream analyses using these mappings have incorrect tissue classifications

## Next Steps
1. Systematic audit of all remaining mappings
2. Check partner data impact scope
3. Create corrected mapping file
4. Version control and documentation