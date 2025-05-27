# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Pipeline Version: v2.2 (Development)

**Status:** Development branch focused on critical validation fixes  
**Based on:** v2.1 production pipeline  
**Target:** Address MAGE tissue mapping, ADNI data_type validation, and ENCODE gene ID format detection issues

## Repository Overview

This is an RNA-seq standardization pipeline that processes multi-omics datasets (ENCODE, GTEx, MAGE, ADNI, ENTEx) into harmonized AnnData (.h5ad) format. The pipeline standardizes gene IDs to GENCODE v24/Ensembl, applies consistent metadata schemas, and creates unified datasets for downstream analysis.

## Pipeline Architecture

The pipeline follows a multi-stage processing workflow:

**Stage 0**: Gene mapping preparation (Entrez-to-Ensembl, GTF processing)
**Stage A**: Primary gene annotation reference map generation
**Stage 1**: Raw data to standardized v1 H5AD conversion
**Stage 1.6**: ENCODE-specific ID mapping generation
**Stage 2**: Enhanced metadata standardization (v1 to v2 H5ADs)
**Stage 2.5**: Gene ID preprocessing and standardization
**Stage 2.6**: Placeholder gene ID fixes
**Stage 2.7**: ENCODE mapping quality analysis
**Stage 3**: Combined dataset creation (sparse format)
**Stage 4**: Validation and reporting

## Key Commands

### Main Pipeline Execution
```bash
# Full pipeline run
./run_rnaseq_pipeline.sh

# Force regeneration of all outputs
./run_rnaseq_pipeline.sh --force

# Force only gene mapping regeneration
./run_rnaseq_pipeline.sh --force-mapping
```

### Environment Setup
```bash
# Initial setup (run once)
./setup_workspace.sh
```

### Pipeline Monitoring
```bash
# Monitor active pipeline run
./monitor_pipeline.sh

# Run in screen session for long-running processes
screen -S rnaseq_pipeline ./run_rnaseq_pipeline.sh
```

### Testing and Validation
```bash
# Validate standardized datasets
python validate_standardized_datasets.py --input-dir /path/to/data --output-file validation_report.json

# Individual component tests
python test_celltype_mapping.py
python test_ensg_extraction.py
```

### Single-Cell Processing (Optional)
```bash
# Process GTEx single-cell data
./single_cell/run_process_gtex_single_cell.sh
```

## Core Components

### Data Processing Scripts
- `standardize_datasets.py`: Main data conversion (raw → v1 H5AD)
- `standardize_metadata.py`: Metadata enhancement (v1 → v2 H5AD)
- `preprocess_dataset_gene_ids.py`: Gene ID standardization
- `create_combined_dataset_all_genes_sparse.py`: Multi-dataset integration

### Gene Mapping Infrastructure
- `gene_id_mapping_reference.py`: Primary reference map generation
- `entrez-to-ensembl-mapping.py`: NCBI-based gene mapping
- `generate_encode_mapping.py`: ENCODE-specific mappings
- `fix_placeholder_ids.py`: Placeholder ID resolution

### Utilities and Validation
- `rnaseq_utils.py`: Shared utility functions and constants
- `validate_standardized_datasets.py`: Dataset validation framework
- `anndata_save_wrapper.py`: Safe H5AD saving with error handling

## Configuration System

### Metadata Configuration
JSON files in `/metadata/json/` define dataset-specific processing:
- `encode_metadata.json`: ENCODE cell line configurations
- `gtex_metadata.json`: GTEx tissue mappings and protocols
- `mage_metadata.json`: 1000 Genomes population mappings
- `tissue_to_uberon.json`: Tissue ontology mappings

### Path Configuration
Input/output paths are configured in `run_rnaseq_pipeline.sh`:
- `INPUT_BASE_DIR`: Source data location
- `OUTPUT_BASE_DIR`: Generated outputs and logs
- Separation allows clean input preservation and organized output structure

## Gene ID Standardization

The pipeline uses a hierarchical gene mapping strategy:
1. **Primary Reference**: GENCODE v24 GTF → Ensembl ID mapping
2. **NCBI Backup**: Entrez → Ensembl mappings for unmapped genes
3. **Dataset-Specific**: ENCODE custom mappings for remaining genes
4. **Placeholder Handling**: Final cleanup of unmapped entries

Target: >99% Ensembl ID mapping rate across all datasets.

## Data Flow Patterns

### Input Processing
Raw data → Standardized v1 → Enhanced v2 → Preprocessed → Combined

### Metadata Enhancement
Basic fields → Ontology terms → Validation → Integration

### Quality Control
Each stage includes validation checkpoints and mapping quality metrics.

## Development Notes

### Error Handling
- Scripts use `anndata_save_wrapper.py` for safe H5AD writes
- Pipeline includes comprehensive logging and error recovery
- Each stage can be run independently with appropriate flags

### Memory Management
- Large datasets processed in chunks where possible
- Sparse matrix formats used for memory efficiency
- Temporary files cleaned up automatically

### Debugging
- Log files include timestamps and stage markers
- Individual test scripts for component validation
- Pipeline monitor for real-time status tracking