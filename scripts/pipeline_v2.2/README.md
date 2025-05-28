# RNA-seq Pipeline v2.2 Scripts

This directory contains the core pipeline scripts organized by functionality.

## Core Pipeline Scripts (Main Directory)

**Main Pipeline:**
- `run_rnaseq_pipeline.sh` - Main pipeline orchestration script
- `rnaseq_utils.py` - Shared utility functions and constants

**Data Processing (Pipeline Stages):**
- `standardize_datasets.py` - Stage 1: Raw data to standardized H5AD
- `standardize_metadata.py` - Stage 2: Metadata enhancement  
- `preprocess_dataset_gene_ids.py` - Stage 2.5: Gene ID preprocessing
- `create_combined_dataset_all_genes_sparse.py` - Stage 3: Combined dataset creation
- `validate_standardized_datasets.py` - Stage 4: Validation

**Gene Mapping Infrastructure:**
- `entrez-to-ensembl-mapping.py` - Stage 0: NCBI gene mapping
- `gene_id_mapping_reference.py` - Stage A: Primary reference map
- `generate_encode_mapping.py` - Stage 1.6: ENCODE-specific mappings
- `fix_placeholder_ids.py` - Stage 2.6: Placeholder ID resolution

**Metadata Enhancement (v2.2 Stages):**
- `integrate_missing_metadata.py` - Stage 2.8: Controlled-access metadata
- `map_developmental_stages.py` - Stage 2.9: Developmental stage mapping
- `integrate_mage_technical_metadata.py` - Stage 3.5: MAGE technical metadata

**Optional Extensions:**
- `process_gtex_single_cell.py` - GTEx single-cell processing

## Supporting Scripts (Subdirectories)

### `analysis/`
Partner deliverables and advanced analysis scripts:
- `analyze_gene_overlap.py` - Stage 4.1: Gene overlap analysis
- `create_subject_level_ethnicity_mapping.py` - Stage 4.2: Ethnicity mappings
- `create_czi_schema_compliant_mapping.py` - Stage 4.2: CZI schema mappings
- Other ethnicity processing scripts

### `validation/`
Testing and validation utilities:
- `test_celltype_mapping.py` - Cell type mapping validation
- `test_ensg_extraction.py` - Gene ID extraction testing
- `check_*.py` - Dataset field validation scripts
- `validate_ethnicity_mapping.py` - Ethnicity mapping validation

### `utilities/`
Infrastructure and monitoring tools:
- `anndata_save_wrapper.py` - Safe H5AD file saving
- `setup_workspace.sh` - Environment setup
- `monitor_*.sh` - Pipeline monitoring scripts

## Configuration Files

- `CLAUDE.md` - Claude Code guidance for this pipeline
- `.gitignore` - Git exclusion rules
- `file_list.md` / `metadata_file_list.md` - Reference documentation

## Usage

Run the complete pipeline:
```bash
./run_rnaseq_pipeline.sh
```

See `./run_rnaseq_pipeline.sh --help` for all available options.