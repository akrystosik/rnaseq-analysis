# Changelog

## [0.1.1] - 2025-04-30

### Fixed
- Fixed parsing of tab-delimited files in ADNI dataset processing
- Added handling for escaped tabs in ADNI CSV files
- Fixed categorical data comparison issues in preprocessed datasets
- Added proper application of metadata from JSON files to all datasets
- Ensured preservation of Ensembl version in original_gene_id field
- Fixed placeholder ID handling in preprocessed datasets

### Added
- New utility script `fix_categorical_columns.py` to fix categorical column comparison issues
- New utility script `apply_dataset_metadata.py` to apply metadata from JSON files
- Improved error handling and reporting throughout the pipeline

### Changed
- Updated standardize_datasets.py to better handle mixed file formats
- Modified processes to preserve original gene IDs with version information
- Improved gene_id handling in preprocessed datasets

## 2025-04-30
### Fixed
- Fixed metadata serialization issues in standardize_datasets.py
- Simplified metadata structure to ensure proper serialization to h5ad files
- Preserved essential metadata fields while removing complex nested structures
- properly parse tab-delimited files via fix to /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py
- categorical column comparison issue with the preprocessed ENCODE dataset


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

## [0.1.2] - 2025-04-30

### Fixed
- Added defensive processing for ADNI files with escaped tabs
- Fixed categorical data comparison issues in preprocessed datasets
- Improved metadata serialization to handle complex nested structures
- Added fix_placeholder_ids.py script to replace placeholder IDs with proper Entrez IDs
- Integrated placeholder ID fix into the pipeline workflow

### Changed
- Made metadata serialization more robust by converting numpy types
- Enhanced categorical column handling to prevent comparison errors
- Updated preprocessing to ensure consistent gene ID representation

## [2025-04-30]
### Fixed
- Fixed gene ID mapping in ENCODE preprocessing to achieve 100% mapping rate
- Modified `preprocess_encode_dataset` function in `preprocess_dataset_gene_ids.py` to use Ensembl IDs as gene identifiers
- Improved extraction of Ensembl IDs from compound identifiers
- Successfully standardized gene IDs across all datasets with:
  - 99.63% proper Ensembl IDs
  - 0.34% Entrez IDs
  - 0.03% spike-in controls
- Fixed placeholder ID handling to ensure all genes have proper identifiers
[2025-04-30] Full Pipeline Test Results
Achievements

Successfully implemented gene ID mapping across all datasets:

ENCODE: 99.63% Ensembl IDs, 0.34% Placeholder IDs, 0.03% Spike-in IDs
MAGE: 100% Ensembl IDs
ADNI: 100% Ensembl IDs


Created combined dataset with 1,388 samples and 40,926 genes
Sparse matrix representation achieved 11% memory savings

Validation Issues Identified

ENCODE and MAGE datasets failed validation due to metadata issues:

Missing or incorrect tissue ontology mappings
Data type validation failures
Missing cell type information


224 placeholder IDs remain in ENCODE dataset (0.34%)

Next Steps

Improve tissue-to-ontology mappings
Add missing cell type annotations
Standardize data type fields across datasets

