# Changelog

## [1.0.0] - 2025-04-30

### Added
- Initial release of the RNA-seq standardization pipeline
- Support for ENCODE, GTEx, MAGE, and ADNI datasets
- Standardization to GENCODE v24 gene annotations and hg38 reference genome
- Combined dataset creation with sparse matrix representation

### Fixed

#### 1. ENCODE Tissue and Cell Line Information
- **Issue**: Missing tissue information in ENCODE data (80.98% of samples had "nan" as tissue value)
- **Fix**: Updated the pipeline to use only the cell_lines directory for ENCODE processing
- **File**: Modified `run_rnaseq_pipeline.sh` to use the correct directory path
- **Impact**: All ENCODE samples now have proper tissue values and cell line identifiers

#### 2. AnnData Saving Issues
- **Issue**: Error when saving AnnData objects: "Can't implicitly convert non-string objects to strings"
- **Fix**: Created a wrapper script `anndata_save_wrapper.py` that properly handles string conversion
- **File**: Added new script and integrated it into the pipeline
- **Impact**: AnnData objects can now be saved without errors

#### 3. Placeholder Gene IDs
- **Issue**: Some gene IDs were being represented as "PLACEHOLDER_*" instead of proper identifiers (0.34%)
- **Fix**: Created `fix_placeholder_ids.py` that converts placeholder IDs to proper Entrez format (ENTREZ:*)
- **File**: Added new script and integrated it as a post-processing step
- **Impact**: No more placeholder IDs in the final datasets

#### 4. Numeric Gene IDs
- **Issue**: The `gene_id` column in preprocessed datasets contained sequential numeric indices
- **Fix**: Modified the preprocessing to use Ensembl IDs as the values for the `gene_id` column
- **File**: Updated `preprocess_dataset_gene_ids.py`
- **Impact**: The `gene_id` column now contains proper biological identifiers (99.63% Ensembl IDs)

### Improved
- Enhanced memory efficiency with sparse matrix representation
- Preserved version information in original_gene_id column
- Added comprehensive metadata validation
- Improved error handling throughout the pipeline
