# RNA-seq Standardization Pipeline Usage Guide

This document provides detailed instructions for running the RNA-seq standardization pipeline.

## Quick Start

```bash
# Run the complete pipeline with default settings
./scripts/run_rnaseq_pipeline.sh

# Run the pipeline with a specific output directory
./scripts/run_rnaseq_pipeline.sh /path/to/output/directory

# Force regeneration of mapping files
./scripts/run_rnaseq_pipeline.sh --force-mapping
```

## Input Data Requirements

The pipeline expects data in the following locations:

- **ENCODE**: `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/raw_data/cell_lines`
- **GTEx**: `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz`
- **MAGE**: `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage`
- **ADNI**: `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/adni_microarray`

## Pipeline Stages

The pipeline runs through the following stages:

1. **Generate Entrez to Ensembl Mapping**: Creates reference mappings between different gene ID systems
2. **Initial Data Conversion**: Processes each dataset into standardized format
3. **Generate Gene ID Reference Mapping**: Creates a comprehensive gene ID reference
4. **Enhanced Metadata Standardization**: Standardizes metadata across datasets
5. **Dataset Preprocessing**: Applies consistent gene IDs across datasets
6. **Fix Placeholder IDs**: Converts any placeholder IDs to proper identifiers
7. **Create Combined Dataset**: Creates datasets with common genes and all genes

## Output Files

The pipeline generates the following outputs in the specified output directory:

- `encode_standardized_v1.h5ad`: Standardized ENCODE dataset
- `gtex_standardized_v1.h5ad`: Standardized GTEx dataset
- `mage_standardized_v1.h5ad`: Standardized MAGE dataset
- `adni_standardized_v1.h5ad`: Standardized ADNI dataset
- `preprocessed_data/`: Folder containing preprocessed datasets with consistent gene IDs
- `combined_all_genes_sparse_standardized.h5ad`: Combined dataset with all genes
- `validation_report_*.json`: Validation report for all standardized datasets

## Examining Results

You can examine the standardized datasets using Python and scanpy:

```python
import scanpy as sc

# Load a standardized dataset
adata = sc.read_h5ad("path/to/standardized_data/encode_standardized_v1.h5ad")

# Check sample metadata
print(adata.obs)

# Check gene information
print(adata.var)

# Check expression values
print(adata.X)
```

## Troubleshooting

If you encounter issues with the pipeline, check the log files in the `logs/` directory. Each pipeline run creates a timestamped log file with detailed information about each step.

Common issues:

1. **Missing input data**: Ensure all required input datasets are available
2. **Memory errors**: For large datasets, increase available memory
3. **File permissions**: Ensure write permissions in the output directory

## Advanced Usage

### Running Individual Steps

You can run individual parts of the pipeline using the component scripts:

```bash
# Just standardize ENCODE data
python scripts/standardize_datasets.py --encode-dir "/path/to/encode/data" --output-dir "/path/to/output"

# Just preprocess a dataset
python scripts/preprocess_dataset_gene_ids.py --data-dir "/path/to/data" --output-dir "/path/to/output"
```

### Customizing the Pipeline

The pipeline behavior can be customized by modifying the JSON configuration files in the `metadata/json/` directory:

- `encode_metadata.json`: ENCODE-specific metadata
- `gene_id_reference_mapping.csv`: Gene ID mapping reference
