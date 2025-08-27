# RNA-seq Analysis Pipeline

This directory contains scripts for processing and analyzing RNA-seq data from multiple sources (ENCODE, GTEx, MAGE, ADNI, and ENTEx).

## Recent Fixes

The pipeline has been fixed to address several issues:

1. **ENTEx Metadata Loading**: Now properly loads all samples from the metadata file (242 samples vs 63 previously)
2. **Variable Naming Consistency**: Fixed variable name mismatches in `rnaseq_utils.py`
3. **Missing Function Implementation**: Added implementation for `map_tissue_to_ontology()`
4. **Path Construction Fix**: Resolved issues with undefined variables in file paths

## Key Scripts

- `standardize_datasets.py`: Processes RNA-seq data into standardized AnnData format
- `rnaseq_utils.py`: Utility functions for data processing and standardization
- `run_rnaseq_pipeline.sh`: Main pipeline runner
- `optimized_gene_expression.py`: Functions for querying gene expression (designed for model training)
- `optimized_gene_expression_query_examples.py`: Examples of how to use the query functions
- `create_combined_dataset.py`: Creates a combined dataset excluding GTEx
- `generate_entex_metadata.py`: Generates ENTEx metadata
- `fix_entex_metadata_format.py`: Fixes metadata format issues

## Data Processing

The pipeline standardizes data from multiple sources:
- ENCODE: 276 samples, 67,151 genes
- GTEx: 19,616 samples, 58,988 genes
- MAGE: 731 samples, 19,428 genes
- ADNI: 650 samples, 17,991 genes
- ENTEx: 203 samples

## Output

The standardized data is saved in AnnData (.h5ad) format in the output directory specified in the pipeline script.
