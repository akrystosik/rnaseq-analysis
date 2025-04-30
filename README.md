# RNA-seq Standardization Pipeline

A pipeline for standardizing RNA-seq gene expression data from multiple sources (ENCODE, GTEx, MAGE, ADNI) into a consistent format for cross-dataset analysis.

## Overview

This pipeline takes diverse RNA-seq datasets and processes them into standardized AnnData objects with:
- Consistent gene identifiers (GENCODE v24)
- Reference genome standardization (hg38)
- Unified metadata structure
- Common expression units

## Features

- **Multi-dataset support**: Processes ENCODE, GTEx, MAGE, and ADNI data
- **GENCODE v24 standardization**: All gene IDs mapped to consistent Ensembl IDs
- **Original ID preservation**: Version information is maintained in original_gene_id column
- **Metadata harmonization**: Consistent metadata fields across all datasets
- **Combined analysis**: Enables creation of a unified dataset with all genes and samples

## Directory Structure

```
/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/
├── scripts/             # Core pipeline scripts
│   ├── standardize_datasets.py            # Initial data conversion
│   ├── preprocess_dataset_gene_ids.py     # Gene ID standardization
│   ├── fix_placeholder_ids.py             # Fix for placeholder gene IDs
│   ├── anndata_save_wrapper.py            # AnnData saving utility
│   └── run_rnaseq_pipeline.sh             # Main pipeline script
├── metadata/            # Metadata and mapping files
│   └── json/            # JSON configuration files
├── standardized_data/   # Output directory for standardized datasets
└── logs/                # Pipeline logs
```

## Key Data Structure

The pipeline produces AnnData objects with the following key components:

1. **obs**: Sample metadata
   - tissue
   - cell_line (for ENCODE)
   - subject_id
   - sex
   - age
   - data_type
   - expression_unit

2. **var**: Gene metadata
   - gene_id: Ensembl ID without version (GENCODE v24)
   - original_gene_id: Original ID with version preserved
   - ensembl_id: Same as gene_id
   - gene_name: Gene symbol
   - gene_type: Gene type (protein_coding, lincRNA, etc.)
   - chromosome: Chromosome location
   - mapping_source: Source of gene ID mapping

3. **X**: Expression matrix (samples × genes)

## Running the Pipeline

Basic usage:

```bash
./scripts/run_rnaseq_pipeline.sh [OUTPUT_DIR]
```

For more detailed usage instructions, see [USAGE.md](USAGE.md).

## Recent Fixes

This pipeline has been thoroughly debugged and improved to ensure reliable processing of multiple RNA-seq datasets. See [CHANGELOG.md](CHANGELOG.md) for details of recent fixes and improvements.

## Requirements

- Python 3.8+
- scanpy
- pandas
- numpy
- anndata

## License

[Include your license information here]
