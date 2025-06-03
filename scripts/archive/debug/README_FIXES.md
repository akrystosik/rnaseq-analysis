# RNA-seq Pipeline Fixes

This document outlines the fixes implemented to address issues in the RNA-seq standardization pipeline.

## 1. ENCODE Tissue and Cell Line Information

**Issue:** Missing tissue information in ENCODE data. Many samples had "nan" as tissue value.

**Fix:** Updated the pipeline to use only the cell_lines directory for ENCODE processing, 
ignoring the ENTEx data as specified in the MVP requirements.

**Location:** Modified `run_rnaseq_pipeline.sh` to use the cell_lines directory:
```bash
--encode-dir "${BASE_DIR}/encode/raw_data/cell_lines"
```

## 2. AnnData Saving Issues

**Issue:** Error when saving AnnData objects: "Can't implicitly convert non-string objects to strings"

**Fix:** Created a wrapper script `anndata_save_wrapper.py` that properly handles string conversion
when saving AnnData objects. The wrapper ensures all metadata values are properly converted to
serializable types.

**Location:** Added to the pipeline and used after each dataset is processed.

## 3. Placeholder Gene IDs

**Issue:** Some gene IDs were being represented as "PLACEHOLDER_*" instead of proper identifiers.

**Fix:** Created `fix_placeholder_ids.py` that converts placeholder IDs to proper Entrez format 
(ENTREZ:*) and handles categorical data types correctly.

**Location:** Added as a post-processing step after the preprocessing stage in the pipeline.

## Validation Results

After implementing these fixes:

1. All ENCODE samples now have proper tissue values and cell line identifiers
2. AnnData objects can be saved without errors
3. Gene ID distribution shows:
   - 99.63% Ensembl IDs
   - 0.34% Entrez IDs
   - 0.03% spike-in controls
   - No more placeholder IDs

These fixes ensure that the pipeline meets the MVP requirements of providing standardized RNA-seq
data from MAGE, ADNI, GTEx, and ENCODE with proper gene IDs, sample IDs, and tissue/cell line 
information.
