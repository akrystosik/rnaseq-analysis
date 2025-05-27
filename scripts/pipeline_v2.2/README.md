# RNA-seq Standardization Pipeline v2.2

## Overview

The RNA-seq Standardization Pipeline v2.2 is a comprehensive, production-ready system for processing multi-omics datasets (ENCODE, GTEx, MAGE, ADNI, ENTEx) into harmonized AnnData (.h5ad) format. This pipeline standardizes gene IDs to GENCODE v24/Ensembl, applies consistent metadata schemas, integrates controlled-access data, and creates unified datasets for downstream analysis.

## Key Features

### ✅ **Core Pipeline (Stages 0-4)**
- **Gene Mapping**: GENCODE v24/Ensembl standardization with 99%+ coverage
- **Data Standardization**: Raw data → standardized H5AD conversion
- **Metadata Enhancement**: Ontology mapping (UBERON, CL, EFO, NCBITaxon)
- **Validation**: Comprehensive dataset validation and quality control

### 🚀 **v2.2 Enhancements (Stages 2.8-4.2)**
- **Controlled-Access Integration**: GTEx demographics from dbGaP
- **Developmental Stage Mapping**: Age → HsapDv ontology (96.5% coverage)
- **Technical Metadata**: MAGE RIN scores and quality metrics
- **Partner Deliverables**: CZI schema-compliant exports
- **Advanced Analytics**: Gene overlap analysis and enhanced validation

## Pipeline Architecture

```
Stage 0    : Entrez-to-Ensembl mapping (NCBI source)
Stage A    : Primary gene reference map generation
Stage 1    : Raw data standardization
Stage 1.6  : ENCODE-specific ID mappings
Stage 2    : Metadata enhancement
Stage 2.5  : Gene ID preprocessing
Stage 2.6  : Placeholder ID fixes
Stage 2.7  : ENCODE mapping quality analysis
Stage 2.8  : Controlled-access metadata integration ⭐
Stage 2.9  : Developmental stage mapping ⭐
Stage 3    : Combined dataset creation
Stage 3.5  : MAGE technical metadata integration ⭐
Stage 4    : Validation
Stage 4.1  : Gene overlap analysis ⭐
Stage 4.2  : Partner deliverables generation ⭐
```

⭐ = New in v2.2

## Quick Start

### Prerequisites
- Python 3.8+ with required packages (scanpy, pandas, numpy, etc.)
- Access to raw data directories
- Optional: GTEx controlled-access files for complete metadata

### Basic Usage

```bash
# Full pipeline run with all enhancements
./run_rnaseq_pipeline.sh

# Force regeneration of all outputs
./run_rnaseq_pipeline.sh --force

# Skip metadata enhancement stages (basic pipeline only)
./run_rnaseq_pipeline.sh --skip-metadata-enhancement

# Show help and available options
./run_rnaseq_pipeline.sh --help
```

### Configuration

Edit the configuration section in `run_rnaseq_pipeline.sh`:

```bash
# Input and output directories
INPUT_BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
OUTPUT_BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq"

# Raw data paths
ENCODE_RAW_DATA_DIR="${INPUT_BASE_DIR}/encode/raw_data"
GTEX_RAW_FILE_INPUT="${INPUT_BASE_DIR}/gtex/raw_data/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"
MAGE_RAW_DATA_DIR="${INPUT_BASE_DIR}/mage"
ADNI_RAW_DATA_DIR="${INPUT_BASE_DIR}/adni_microarray"
```

## Data Requirements

### Core Data Files
- **ENCODE**: Cell line RNA-seq data (`/encode/raw_data/cell_lines/`)
- **GTEx**: Gene expression TPM file (`GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz`)
- **MAGE**: 1000 Genomes RNA-seq data (`/mage/`)
- **ADNI**: Microarray data (`/adni_microarray/`)

### Metadata Files
- **GENCODE v24 GTF**: Gene annotations (`gencode.v24.annotation.gtf.gz`)
- **ADNI Demographics**: Subject demographics (`PTDEMOG_25Apr2025.csv`)
- **ADNI Diagnosis**: Diagnosis summaries (`DXSUM_30Apr2025.csv`)
- **MAGE Population**: 1000 Genomes PED file (`20130606_g1k_3202_samples_ped_population.txt`)

### Optional Controlled-Access Files
- **GTEx Phenotypes**: `phs000424.v10.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz`
- **MAGE Technical**: `sample.metadata.MAGE.v1.0.txt`

## Output Structure

```
output_directory/
├── standardized_data/run_TIMESTAMP/
│   ├── *_standardized_v2.h5ad           # Enhanced H5AD files
│   ├── combined_dataset_all_genes_sparse.h5ad  # Combined dataset
│   ├── validation_reports/               # Validation results
│   └── partner_deliverables/            # CZI schema exports
├── preprocessed_data/run_TIMESTAMP/
│   └── *_standardized_preprocessed.h5ad # Gene ID standardized files
├── metadata/
│   ├── json_generated/                  # Generated metadata
│   └── gene_mapping_generated/          # Gene mapping files
└── logs/pipeline_runs/                  # Execution logs
```

## Key Outputs

### Dataset Files
- **Individual H5ADs**: `{dataset}_standardized_preprocessed.h5ad`
- **Combined Dataset**: `combined_dataset_all_genes_sparse.h5ad` (21,004 samples × 68,621 unique genes)

### Validation Reports
- **Main Report**: `validation_report.json` (100% pass rate)
- **Gene Overlap**: Gene intersection analysis across datasets
- **Metadata Coverage**: Completeness assessment

### Partner Deliverables
- **Ethnicity Mapping**: `subject_ethnicity_mapping_with_ontology.csv` (99.8% coverage)
- **CZI Schema Files**: Schema v3.0.0 compliant outputs
- **Technical Reports**: RNA quality assessments

## Metadata Standards

### Ontology Compliance
- **Tissues**: UBERON terms (e.g., `UBERON:0000178` for blood)
- **Cell Types**: CL terms (e.g., `CL:0000542` for lymphoblast)
- **Assays**: EFO terms (e.g., `EFO:0009922` for RNA-seq)
- **Species**: NCBITaxon terms (`NCBITaxon:9606` for human)
- **Ethnicity**: HANCESTRO terms (e.g., `HANCESTRO:0005` for European)
- **Development**: HsapDv terms (e.g., `HsapDv:0000087` for adult stage)

### CZI Schema Compliance
Full compliance with CZI single-cell curation schema v3.0.0:
- Required fields: `assay_ontology`, `cell_type_ontology_term_id`, `tissue_ontology_term_id`
- Recommended fields: `self_reported_ethnicity_ontology_term_id`, `sex`, `organism_ontology_term_id`
- Optional fields: `development_stage_ontology_term_id`, `disease_ontology_term_id`

## Quality Metrics

### Gene ID Standardization
- **ADNI**: 100% Ensembl coverage (17,991 genes)
- **ENCODE**: 98.3% Ensembl coverage (64,499/65,586 genes)
- **GTEx**: 100% Ensembl coverage (58,988 genes)
- **MAGE**: 100% Ensembl coverage (19,428 genes)

### Metadata Completeness
- **Ethnicity**: 99.8% coverage with HANCESTRO ontology
- **Developmental Stage**: 96.5% coverage with HsapDv ontology
- **Sex**: 100% coverage across all datasets
- **Age**: 100% coverage where privacy allows (3/4 datasets)

### RNA Quality Assessment
- **GTEx RIN**: Mean 7.3 (range 2.3-10.0)
- **MAGE RIN**: Mean 9.7 (range 7.9-10.0), 96.6% excellent quality

## Advanced Features

### Gene Overlap Analysis
- **Core genes**: 13,244 genes present in all datasets
- **Unique genes**: 68,621 total across all datasets
- **Redundancy factor**: 2.36x (shows high overlap)
- **RNA-seq similarity**: GTEx ∩ ENCODE = 81.6% Jaccard similarity

### Controlled-Access Integration
- **GTEx Demographics**: Sex and ethnicity from dbGaP phenotype files
- **ADNI Age Calculation**: Birth date to visit date calculation
- **Privacy Protection**: Compliant with data use agreements

## Troubleshooting

### Common Issues

1. **Missing Input Files**
   ```bash
   # Check file paths in configuration section
   ls -la $GTEX_RAW_FILE_INPUT
   ```

2. **Permission Errors**
   ```bash
   # Ensure write permissions
   chmod -R 755 $OUTPUT_BASE_DIR
   ```

3. **Memory Issues**
   ```bash
   # Monitor memory usage
   ./monitor_v2.2.sh
   ```

### Debug Mode
```bash
# Enable verbose logging
export PYTHONPATH="${SCRIPTS_DIR}:${PYTHONPATH}"
python -u script_name.py --verbose
```

### Validation Failures
Check validation reports for specific issues:
```bash
cat $OUTPUT_DIR/validation_reports/validation_report.json
```

## Development

### Testing
```bash
# Individual component tests
python test_celltype_mapping.py
python test_ensg_extraction.py

# Validation tests
python validate_ethnicity_mapping.py
```

### Adding New Datasets
1. Create standardization script following existing patterns
2. Add to `standardize_datasets.py`
3. Update metadata configuration files
4. Add validation rules
5. Test with pilot data

## Documentation

### Key Files
- **CLAUDE.md**: Development guidance and architecture
- **VERSION.md**: Version history and changes
- **METADATA_PROGRESS_VALIDATION_PLAN.md**: Validation strategy
- **v2.2_partner_presentation.ipynb**: Results showcase

### Partner Resources
- **ANCESTRY_REPORT_FOR_PARTNER.md**: Ethnicity analysis
- **FINAL_DELIVERY_STATUS.md**: Completion summary
- **SAMPLE_ETHNICITY_SUMMARY.md**: Sample-level ethnicity

## Citation

When using this pipeline, please cite:

```
RNA-seq Standardization Pipeline v2.2
Multi-omics Data Harmonization with Comprehensive Metadata Integration
[Institution/Lab Name], 2025
```

## License

[Specify license terms]

## Support

For issues and questions:
- Check troubleshooting section above
- Review log files in `$OUTPUT_DIR/logs/`
- Examine validation reports for specific errors

## Version History

- **v2.2**: Metadata enhancement, controlled-access integration, partner deliverables
- **v2.1**: Enhanced validation, gene mapping improvements
- **v2.0**: Core pipeline with standardization and validation
- **v1.0**: Initial implementation

---

**Pipeline v2.2 Status**: ✅ Production Ready  
**Validation**: 100% pass rate (4/4 datasets)  
**Coverage**: 21,004 samples, 68,621 unique genes  
**Partner Ready**: Complete with deliverables and documentation