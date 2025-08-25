# Croissant Metadata Generation for H5AD Genomics Datasets

## Overview

This directory contains Python scripts and generated Croissant JSON-LD metadata files for four standardized genomics datasets. Croissant is the MLCommons standard for machine learning dataset metadata, designed to make datasets more discoverable and interoperable.

## Generated Files

### Scripts
- `croissant_generator.py` - Main generator script with dataset-specific metadata
- `examine_h5ad_structure.py` - Utility to examine h5ad file structures  
- `validate_croissant.py` - Validation script using mlcroissant library
- `h5ad_structure_summary.json` - Detailed structural analysis of all h5ad files

### Croissant Metadata Files
- `adni_standardized_preprocessed_croissant.jsonld` - ADNI Alzheimer's dataset
- `encode_standardized_preprocessed_croissant.jsonld` - ENCODE functional genomics dataset
- `gtex_standardized_preprocessed_croissant.jsonld` - GTEx tissue expression dataset  
- `mage_standardized_preprocessed_croissant.jsonld` - MAGE population genetics dataset

## Dataset Information

### GTEx (Genotype-Tissue Expression)
- **Samples**: 19,616 samples across 49 tissues from 838 individuals
- **Genes**: 58,988 genes
- **Key Publications**: 
  - Aguet, F. et al. Science 369, 1318-1330 (2020)
  - GTEx Consortium. Nature 550, 204-213 (2017)
- **License**: dbGaP Authorized Access
- **Description**: Multi-tissue RNA-seq data for understanding genetic regulation

### ENCODE (Encyclopedia of DNA Elements) 
- **Samples**: 7 samples from various cell lines
- **Genes**: 65,586 genes
- **Key Publications**:
  - ENCODE Consortium. Nature 583, 699-710 (2020)
  - ENCODE Consortium. Nature 489, 57-74 (2012)
- **License**: CC0
- **Description**: Functional genomics data from cell lines and tissues

### ADNI (Alzheimer's Disease Neuroimaging Initiative)
- **Samples**: 650 samples from brain tissue
- **Genes**: 17,991 genes  
- **Key Publications**:
  - Petersen RC, et al. Neurology 74, 201-209 (2010)
- **License**: ADNI Data Use Agreement
- **Description**: Microarray data from Alzheimer's disease research

### MAGE (Multi-ancestry Analysis of Gene Expression)
- **Samples**: 731 samples from lymphoblastoid cell lines
- **Genes**: 19,428 genes
- **Key Publications**:
  - Taylor DJ, et al. Nature (2024). DOI: 10.1038/s41586-024-07708-2
- **License**: CC0  
- **Description**: Population genetics study across 26 globally-distributed populations

## Croissant Metadata Structure

Each generated metadata file follows the MLCommons Croissant 1.0 specification and includes:

### Core Metadata
- **@context**: Official Croissant JSON-LD context
- **@type**: sc:Dataset (Schema.org Dataset)
- **conformsTo**: http://mlcommons.org/croissant/1.0
- **name**: Human-readable dataset name
- **description**: Detailed description with sample/gene counts
- **url**: Official dataset homepage
- **license**: Dataset license information
- **creator**: Dataset creator/organization
- **citation**: Key publications
- **keywords**: Search keywords

### Genomics-Specific Metadata
- **sc:organism**: "Homo sapiens" 
- **sc:assayType**: Assay technology (RNA sequencing, Microarray)
- **sc:tissueType**: Tissue types included

### Distribution Information
- **distribution**: FileObject describing the h5ad file
  - **name**: Filename
  - **description**: File description  
  - **contentUrl**: File location URL
  - **encodingFormat**: "application/x-hdf5"
  - **contentSize**: File size in bytes

### Data Structure Description
- **recordSet**: Describes the internal structure
  - **gene_expression_matrix**: Main expression data
  - **sample_metadata**: Sample-level annotations
  - **gene_metadata**: Gene-level annotations

### Processing Information
- **sc:processingInformation**: Processing pipeline details
  - **standardization**: "Cell by Gene schema v3.0.0"
  - **preprocessing**: Processing description
  - **processing_date**: Processing timestamp
  - **reference_genome**: Reference genome version
  - **gencode_version**: Gene annotation version

### Dataset Statistics  
- **sc:datasetStatistics**: Key statistics
  - **numberOfSamples**: Total sample count
  - **numberOfGenes**: Total gene count
  - **dataShape**: Matrix dimensions [samples, genes]
  - **tissueTypes**: List of tissue types (limited to top 10)
  - **subjectDemographics**: Demographic summary

## Usage

### Generate Croissant Metadata
```bash
python3 croissant_generator.py
```

### Examine H5AD Structure
```bash
python3 examine_h5ad_structure.py
```

### Validate Generated Files
```bash
python3 validate_croissant.py
```

## Technical Implementation

### Architecture
The generator uses a modular design with:
- `DatasetInfo` class containing literature-derived metadata for each dataset
- `CroissantGenerator` class with methods for h5ad parsing and JSON-LD generation
- Automated extraction of metadata from AnnData objects
- Type-safe JSON serialization handling numpy/pandas types

### Design Principles
- **Literature-based**: Dataset information sourced from original publications
- **Cell by Gene Schema**: Follows CZI single-cell curation standards  
- **First principles**: Deep understanding of each dataset's provenance
- **Reproducible**: Version-controlled scripts with deterministic output
- **Validated**: Schema compliance checking with official tools

### Key Features
- Handles multiple genomics dataset types (RNA-seq, microarray)
- Extracts comprehensive metadata from h5ad files
- Generates compliant Croissant JSON-LD with official context
- Includes genomics-specific vocabulary extensions
- Provides detailed field descriptions for all metadata columns
- Handles tissue ontology corrections and standardization

## Validation Notes

The generated files are structurally valid JSON-LD but may have validation warnings for:
- Missing file hashes (md5/sha256) - acceptable for metadata-only descriptions
- Missing `source` fields - the files describe existing h5ad datasets rather than providing loading instructions
- Minor context differences - uses standard Croissant context with minimal extensions

These are design choices appropriate for describing existing preprocessed datasets rather than providing full data loading pipelines.

## Future Enhancements

1. **File hashing**: Add SHA256 hash computation for file integrity
2. **Source mapping**: Add source field mappings for full data loading compliance
3. **Responsible AI**: Integrate Croissant RAI vocabulary extensions  
4. **Validation**: Address minor schema compliance issues
5. **Automation**: Integrate into preprocessing pipeline

## References

- [Croissant Format Specification](https://docs.mlcommons.org/croissant/docs/croissant-spec.html)
- [MLCommons Croissant Repository](https://github.com/mlcommons/croissant)
- [Cell by Gene Schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md)

Generated on: 2025-08-25 using Python 3.10 and mlcroissant 1.0.22