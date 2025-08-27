# RNA-seq Analysis Pipeline

Multi-omics RNA-seq standardization and analysis pipeline integrating ADNI, ENCODE, GTEx, and MAGE datasets with comprehensive metadata validation.

## Overview

This repository contains a standardized RNA-seq analysis pipeline designed for multi-dataset integration and comparative analysis across diverse biological contexts. The pipeline harmonizes expression data from major genomics consortiums into a unified analytical framework.

## Datasets

### Integrated RNA-seq Datasets
- **ADNI** - Alzheimer's Disease Neuroimaging Initiative expression data
- **ENCODE** - Encyclopedia of DNA Elements cell line expression
- **GTEx** - Genotype-Tissue Expression project tissue samples  
- **ENTEx** - Encyclopedia of DNA Elements tissue expression
- **MAGE** - Multi-ancestry Analysis of Gene Expression

## Pipeline Features

### Data Standardization
- **Gene ID harmonization** across platforms (Ensembl, Entrez, Symbol)
- **Metadata standardization** with controlled vocabularies
- **Quality control** and batch effect correction
- **Cross-platform normalization**

### Multi-omics Integration
- Paired with **WGS variant data** for comprehensive analysis
- **Expression quantitative trait loci (eQTL)** analysis ready
- **Population genetics** integration capabilities

### Single Cell Analysis
- **scRNA-seq processing** pipeline included
- **Genotype-expression matching** for personalized genomics
- **Cell type annotation** and analysis tools

## Data Processing

### Standardized Outputs
All datasets processed to consistent formats:
- **AnnData (.h5ad)** format for Python analysis
- **Standardized metadata** with ontology terms
- **Gene expression matrices** with unified gene identifiers
- **Sample metadata** with demographic and clinical annotations

### Quality Assurance
- Comprehensive validation reports
- Cross-dataset compatibility checks  
- Metadata schema compliance verification
- Statistical quality metrics

## Quick Start

```bash
# Clone repository  
git clone https://github.com/akrystosik/rnaseq-analysis.git
cd rnaseq-analysis

# Install dependencies
pip install -r requirements.txt

# Run standardization pipeline
bash scripts/run_rnaseq_pipeline.sh

# Generate combined dataset
python scripts/create_combined_dataset.py
```

## Structure

```
rnaseq-analysis/
├── data/                    # Processed datasets
│   ├── preprocessed_data/   # Standardized .h5ad files
│   └── standardized_data/   # Quality-controlled data
├── scripts/                 # Processing and analysis scripts
├── metadata/               # Gene mappings and JSON metadata  
├── docs/                   # Comprehensive documentation
├── reports/                # Validation and QC reports
└── results/                # Analysis outputs and figures
```

## Documentation

- [Methods](docs/RNASEQ_STANDARDIZATION_METHODS.md)
- [Pipeline Usage](docs/USAGE.md)
- [Metadata Schema](docs/METADATA_PROGRESS_VALIDATION_PLAN.md)
- [Version History](docs/VERSION.md)

## Validation

All processed data includes:
- ✅ **Metadata validation** against schema requirements
- ✅ **Cross-platform gene ID** consistency checks  
- ✅ **Statistical quality control** metrics
- ✅ **Reproducibility testing** with documented versions

## Citation

If you use this pipeline or processed datasets, please cite:

```bibtex
@software{rnaseq_analysis_pipeline,
  title = {Multi-omics RNA-seq Standardization Pipeline},
  author = {Krystosik, Amy},
  year = {2025}, 
  url = {https://github.com/akrystosik/rnaseq-analysis}
}
```


---
**Note**: This pipeline specializes in RNA-seq datasets with paired genomic data for integrated multi-omics analysis.