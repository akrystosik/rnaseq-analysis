# RNA-seq Data Standardization Pipeline - Methods

**Last Updated: ** May 22, 2025

## 1. Overview

This document outlines the methods used to process and standardize diverse RNA-sequencing (RNA-seq) and microarray expression datasets into a harmonized collection of AnnData objects. The goal is to create analysis-ready datasets with consistent gene identifiers, rich ontology-mapped metadata, and clear provenance, suitable for downstream integrative analyses, including linkage with Whole Genome Sequencing (WGS) data.

The pipeline handles bulk RNA-seq, single-cell RNA-seq (via pseudobulking), and microarray data from multiple public and internal sources.

## 2. Data Sources and Types

The pipeline currently processes data from the following sources:

| Dataset          | Original Data Type(s)             | Assay Type(s)        | Processed As | Key Characteristics                                                                |
| :--------------- | :-------------------------------- | :------------------- | :----------- | :--------------------------------------------------------------------------------- |
| **ENCODE**       | Bulk RNA-seq (TSV)                | RNA-seq (polyA+, total) | Bulk RNA-seq | Cell lines; GRCh38/GENCODE v29 originally.                                         |
| **ENTEx**        | Bulk RNA-seq (TSV)                | RNA-seq (polyA+, total) | Bulk RNA-seq | Primary tissues from ENCODE donors; GRCh38/GENCODE v29 originally.                   |
| **GTEx**         | Bulk RNA-seq (GCT), snRNA-seq (H5AD) | RNA-seq (polyA+)     | Bulk, SC PB  | Multiple tissues from many donors; GRCh38/GENCODE v39 (bulk), v26 (SC) originally. |
| **MAGE**         | Bulk RNA-seq (CSV)                | RNA-seq (total)      | Bulk RNA-seq | Lymphoblastoid Cell Lines (LCLs) from 1000 Genomes Project; GRCh38/GENCODE v38.      |
| **ADNI**         | Microarray (TXT/CSV)              | Microarray           | Microarray   | Blood samples; Longitudinal study; Affymetrix Human Genome U219 Array.               |

*SC PB = Single-Cell Pseudobulk*

## 3. Metadata Standards and Ontologies

A core goal is to harmonize metadata to common standards and map relevant fields to established ontologies. Dataset-specific metadata details are often specified in JSON configuration files (e.g., `gtex_metadata.json`, `adni_metadata.json`) located in `metadata/json/`.

| Metadata Field                                | Standard / Ontology                                     | Primary Script(s) for Standardization                               | Notes                                                                                                                               |
| :-------------------------------------------- | :------------------------------------------------------ | :------------------------------------------------------------------ | :---------------------------------------------------------------------------------------------------------------------------------- |
| **Dataset Identifier** (`dataset`)            | Dataset short name (e.g., 'encode', 'gtex')             | `standardize_datasets.py`, `rnaseq_utils.standardize_metadata`    | Categorical.                                                                                                                        |
| **Sample Identifier** (`sample_id`)           | Unique string per sample                                | `standardize_datasets.py` (per dataset logic)                       | Serves as `.obs.index`.                                                                                                             |
| **Subject/Donor Identifier** (`subject_id`)   | Consistent ID for linking across assays/datasets        | `standardize_datasets.py` (per dataset logic)                       | Aim for stability. Links to WGS data.                                                                                               |
| **Tissue** (`tissue`, `tissue_ontology`)      | UBERON                                                  | `standardize_datasets.py`, `rnaseq_utils.standardize_metadata` (via `TissueOntologyMapper`), `metadata/json/tissue_to_uberon.json` | `tissue` is standardized string, `tissue_ontology` is UBERON ID.                                                                    |
| **Cell Type** (`cell_type`, `cell_type_ontology_term_id`) | Cell Ontology (CL)                                      | (Primarily for SC/Cell Lines) `standardize_metadata.py` (planned for `celltype_to_cl.json`) | For cell lines and single-cell derived data.                                                                                        |
| **Assay Type** (`data_type`, `assay_ontology`) | Experimental Factor Ontology (EFO)                      | `standardize_datasets.py`, `rnaseq_utils.standardize_metadata`, `metadata/json/assay_to_efo.json` | `data_type` is standardized string (e.g., "RNA-seq"), `assay_ontology` is EFO ID.                                                     |
| **Expression Unit** (`expression_unit`)       | Standardized string (e.g., "TPM", "Counts", "Normalized intensity") | `standardize_datasets.py` (per dataset logic)                       | Categorical.                                                                                                                        |
| **Age** (`age`, `developmental_stage_ontology`) | Age as string, HsapDv for developmental stage         | `standardize_datasets.py`, `rnaseq_utils.standardize_metadata`, `metadata/json/age_to_hsapdv.json` | `age` is string representation.                                                                                                     |
| **Sex** (`sex`)                               | "male", "female", "unknown" (aligns with CellXGene)   | `standardize_datasets.py`, `rnaseq_utils.standardize_metadata`, `metadata/json/sex_standardization.json` | Categorical.                                                                                                                        |
| **Species** (`species`, `species_ontology`)   | NCBI Taxonomy (e.g., `NCBITaxon:9606` for human)      | `standardize_datasets.py`, `rnaseq_utils.standardize_metadata`, `metadata/json/species_to_taxon.json` | Typically "human".                                                                                                                  |
| **Self-Reported Ethnicity** (`self_reported_ethnicity`, `self_reported_ethnicity_ontology_term_id`) | Standard string labels, HANCESTRO ontology terms (or "multiethnic", "unknown") | `rnaseq_utils.standardize_metadata`, `metadata/json/ethnicity_to_hancestro.json` | Derived from source `race` and `is_hispanic_or_latino` fields.                                                                        |
| **Race** (`race`)                             | NIH categories (lowercase, e.g., "white")             | `rnaseq_utils.standardize_metadata`                                 | Standardized string.                                                                                                                |
| **Is Hispanic or Latino** (`is_hispanic_or_latino`) | "hispanic or latino", "not hispanic or latino", "unknown or not reported" | `rnaseq_utils.standardize_metadata`                                 | Standardized string.                                                                                                                |
| **Disease/Diagnosis** (`diagnosis`, `diagnosis_code`) | Dataset-specific, aiming for MONDO where applicable     | `standardize_datasets.py` (ADNI specific logic), `adni_metadata.json` | ADNI has `diagnosis_code` and `diagnosis` label. Longitudinal history in `.uns`. Cell lines have disease info in `encode_metadata.json`. |
| **Gene Identifiers**                          | Ensembl Gene IDs (ENSG), GENCODE v24 (harmonized)     | See Section 4.2                                                     | `.var.index` and `var['gene_id']` are primary. `var['ensembl_id']` stores the clean Ensembl ID.                                      |
| **Reference Genome**                          | GRCh38/hg38 (harmonized)                                | Stored in `.uns`                                                    | `harmonized_reference_genome` and `original_reference_genome`.                                                                      |
| **GENCODE Version**                           | v24 (harmonized)                                        | Stored in `.uns`                                                    | `harmonized_gencode_version` and `original_gencode_version`.                                                                        |

## 4. Transformations and Processing

The pipeline is executed via `run_rnaseq_pipeline.sh` and involves several Python scripts.

### 4.1. Gene Identifier Mapping and Annotation Resources

**4.1.1. Entrez to Ensembl Mapping (NCBI Source)**
*   **Script:** `entrez-to-ensembl-mapping.py`
*   **Purpose:** Downloads and processes the `gene2ensembl.gz` file from NCBI FTP to create a direct mapping between Entrez Gene IDs and Ensembl Gene IDs for human.
*   **Output:** `metadata/json/entrez_to_ensembl_mapping.csv`
*   **Usage:** This mapping is a key input for generating the primary gene annotation reference.

**4.1.2. GENCODE v24 GTF Preparation**
*   **Source:** `gencode.v24.annotation.gtf.gz` (manual download assumed to be in `metadata/gene_mapping/`).
*   **Action:** The pipeline script unzips this file if the unzipped version is not present or outdated.
*   **Output:** `metadata/gene_mapping/gencode.v24.annotation.gtf`
*   **Usage:** This GTF file is parsed to extract gene information for GENCODE v24, forming the backbone of our harmonized gene annotations.

**4.1.3. Primary Gene Annotation & Reference Map Generation**
*   **Script:** `gene_id_mapping_reference.py`
*   **Purpose:** Creates a comprehensive gene ID reference mapping database (`gencode_v24_complete_mapping.csv`) that serves as the single source of truth for gene identifiers and basic annotations (gene name, type, chromosome).
*   **Inputs:**
    *   GENCODE v24 GTF file (parsed gene information).
    *   `entrez_to_ensembl_mapping.csv` (from step 4.1.1).
    *   Raw data directories for ENCODE and ENTEx (to identify numeric/Entrez IDs present in those datasets that might need mapping or placeholder creation).
*   **Output:**
    *   `metadata/gene_mapping/gencode_v24_complete_mapping.csv`: The primary reference file. Contains Ensembl IDs (versionless), original Ensembl IDs (with version from GTF), gene names, types, chromosome, start, end, strand, and mapped numeric/Entrez IDs where available.
    *   A JSON version is also saved.
*   **Usage:** This map is critically used by:
    *   `standardize_datasets.py` (via `rnaseq_utils.load_gencode_mapping`) for initial gene annotation in Stage 1.
    *   `preprocess_dataset_gene_ids.py` (as `--reference-mapping`) for final gene ID standardization and detailed `.var` population in Stage 2.5.

**4.1.4. ENCODE-specific ID to Ensembl Mapping (Auxiliary)**
*   **Script:** `generate_encode_mapping.py`
*   **Purpose:** Processes ENCODE gene IDs from raw ENCODE TSV files (which can be complex, e.g., `gene_id|transcript_id|...` or numeric) and attempts to extract standard Ensembl IDs (versionless) or Entrez IDs.
*   **Output:** `metadata/gene_mapping/encode_specific_mappings/encode_id_to_ensembl_mapping.csv`.
*   **Usage:** Used by `preprocess_dataset_gene_ids.py` specifically for ENCODE data to help resolve its original complex gene identifiers to more standard forms before final mapping against the primary reference.

### 4.2. Bulk RNA-seq and Microarray Processing (Iterative Stages)

**Stage 1: Initial Data Conversion (Raw Data -> `*_standardized_v1.h5ad`)**

*   **Script:** `standardize_datasets.py`
*   **Core Actions:**
    1.  **Data Ingestion:** Reads raw expression data files (TSV for ENCODE/ENTEx, GCT for GTEx, CSV for MAGE, TXT/CSV for ADNI).
    2.  **Metadata Extraction:** Extracts basic sample metadata from filenames, directory structures, or dedicated metadata files (e.g., GTEx attribute files, ADNI demographics/diagnosis files).
    3.  **Expression Matrix Unification:** For each dataset, combines expression data from multiple files/samples into a unified matrix (genes x samples). Handles gene ID aggregation (e.g., if multiple original IDs map to the same standardized Ensembl ID, `max` expression is typically taken).
    4.  **Initial Gene ID Standardization:** Uses `rnaseq_utils.standardize_ensembl_id` to remove version numbers from Ensembl IDs.
    5.  **GENCODE Annotation:** Adds basic gene annotations (gene name, type, chromosome) to the variable metadata (`.var`) by looking up standardized Ensembl IDs in the primary gene annotation reference (`gencode_v24_complete_mapping.csv`) via `rnaseq_utils.load_gencode_mapping` and `rnaseq_utils.add_gencode_annotations`.
    6.  **Initial Metadata Standardization:** Applies general metadata standardization rules using `rnaseq_utils.standardize_metadata`. This includes:
        *   Ensuring core metadata fields exist.
        *   Standardizing values (e.g., sex to "male"/"female"/"unknown").
        *   Mapping tissue, assay, age, species to ontology terms using JSON mapping files (`tissue_to_uberon.json`, etc.).
        *   Deriving `race`, `is_hispanic_or_latino`, and `self_reported_ethnicity` fields.
    7.  **Dataset-Specific JSON Configuration:** Loads and applies dataset-specific configurations from `metadata/json/<dataset_name>_metadata.json` using `rnaseq_utils.load_dataset_specific_metadata` and `rnaseq_utils.apply_dataset_specific_metadata`. This can override/supplement general metadata and add dataset-specific `.uns` entries.
    8.  **AnnData Object Creation:** Creates an AnnData object (`.X` stores expression, `.obs` stores sample metadata, `.var` stores gene metadata).
    9.  **Saving:** Saves the processed data as `<dataset_name>_standardized_v1.h5ad` using a robust `save_anndata` function that handles serialization of complex types in `.uns`.

*   **Dataset-Specific Logic in Stage 1:**
    *   **ENCODE/ENTEx:** `process_encode_data` handles both. Reads TSV files. Uses `entex_metadata.json` for detailed ENTEx sample/donor info and `encode_metadata.json` for cell line info and defaults.
    *   **GTEx (Bulk):** `process_gtex_data` reads GCT file. `load_gtex_metadata` merges sample and subject attribute files.
    *   **MAGE:** `process_mage_data` reads CSV files. Integrates metadata from the 1000 Genomes PED file for sex and population/ethnicity.
    *   **ADNI:** `process_adni_data` reads TXT/CSV expression files. Integrates `PTDEMOG` (demographics) and `DXSUM` (longitudinal diagnosis) files. Calculates `age_at_expression_sampling`. Stores full longitudinal diagnosis history in `.uns['longitudinal_diagnoses_adni']`.

**Stage 2: Enhanced Metadata Standardization (`*_standardized_v1.h5ad` -> `*_standardized_v2.h5ad`)**

*   **Script:** `standardize_metadata.py` (this is the main script for Stage 2, distinct from the util function).
*   **Purpose:** To further refine and ensure consistency of metadata, primarily by re-applying dataset-specific JSON configurations and ensuring all core ontology mappings and harmonized versions are correctly set in the AnnData objects produced by Stage 1.
*   **Core Actions:**
    1.  Loads `<dataset_name>_standardized_v1.h5ad`.
    2.  **Re-applies Dataset-Specific JSON:** Uses `load_dataset_specific_metadata` and `apply_dataset_specific_metadata` from `rnaseq_utils.py` to ensure `.uns` and `.obs` reflect the latest JSON configurations.
    3.  **Ensures Core Fields & Harmonization Standards:** Explicitly sets/validates `harmonized_gencode_version` (to "v24") and `harmonized_reference_genome` (to "hg38") in `.uns`. Ensures core `.obs` columns like `dataset`, `sample_id`, `species`, `sex`, `tissue`, `data_type`, `expression_unit` exist and have appropriate default values or dtypes.
    4.  **Tissue Ontology Mapping:** Applies `TissueOntologyMapper` (defined within `standardize_metadata.py`) to populate/update `tissue_ontology` and `tissue_ontology_confidence` in `.obs` based on `metadata/json/tissue_to_uberon.json`.
    5.  **Saving:** Saves the enhanced dataset as `<dataset_name>_standardized_v2.h5ad` using the same robust `save_anndata` logic.

**Stage 2.5: Dataset Gene ID Preprocessing (`*_standardized_v2.h5ad` -> `*_standardized_preprocessed.h5ad`)**

*   **Script:** `preprocess_dataset_gene_ids.py`
*   **Purpose:** To rigorously standardize gene identifiers in the `.var` DataFrame of each dataset, ensuring they align with the primary gene annotation reference and that `.var` is comprehensively populated.
*   **Core Actions:**
    1.  Loads `<dataset_name>_standardized_v2.h5ad`.
    2.  **Reference Mapping:** Uses `gencode_v24_complete_mapping.csv` (generated in Stage A) as the definitive reference.
    3.  **Gene ID Standardization & Annotation:**
        *   Iterates through `adata.var_names` (original gene IDs from Stage 1/2).
        *   **ENCODE Special Handling (`preprocess_encode_dataset`):** Uses `encode_id_to_ensembl_mapping.csv` (from Stage 1.6) to first resolve ENCODE's complex/numeric original IDs to Ensembl or Entrez IDs. Then, these resolved IDs are mapped against the primary reference. The `.var` index is set to the final mapped ID (Ensembl, ENTREZ:, or gSpikein), made unique.
        *   **ENTEx Special Handling (`preprocess_entex_dataset`):** Maps original IDs (Ensembl with/without version, numeric/Entrez, gSpikein) to standard Ensembl IDs (versionless) or keeps them as ENTREZ:/gSpikein if appropriate, using the primary reference.
        *   **Other Datasets (GTEx, MAGE, ADNI - `preprocess_other_dataset`):** Assumes these datasets already use Ensembl-like IDs. Strips versions and maps against the primary reference.
        *   Populates `.var` columns:
            *   `gene_id`: The final standardized gene ID used as the new index (e.g., `ENSG00000XXXXXX`, `ENTREZ:12345`, `gSpikein_XYZ`). This often matches `ensembl_id`.
            *   `ensembl_id`: The clean, versionless Ensembl ID (ENSG), or an ENTREZ ID if Ensembl mapping failed but Entrez was available, or the original ID if it's a spike-in or unmapped.
            *   `original_id_from_input` (ENCODE specific): The numeric string from the input AnnData's `.var_names`.
            *   `original_id_from_source` (ENCODE specific): The complex ID string from the `original_ids` column of the Stage 1/2 AnnData.
            *   `original_gene_id` (ENTEx, others): The original gene ID from the input AnnData's `.var_names`.
            *   `gene_name`, `gene_type`, `chromosome`: From the primary reference map.
            *   `mapping_source`, `mapping_confidence`: Indicate how the mapping was achieved.
    4.  **Index Update:** Sets `adata.var.index` to the newly standardized `gene_id` values, ensuring uniqueness (e.g., using `anndata.utils.make_index_unique`).
    5.  **Saving:** Saves the dataset as `<dataset_name>_standardized_preprocessed.h5ad`. For ENCODE, it uses the `anndata_save_wrapper.py` script for potentially safer saving due to complex `.uns` structures.

**Stage 2.6: Fix Placeholder Gene IDs**

*   **Script:** `fix_placeholder_ids.py`
*   **Purpose:** To convert any gene IDs in `.var['gene_id']` or `.var['ensembl_id']` that were temporarily marked as `PLACEHOLDER_<numeric_id>` (typically from `gene_id_mapping_reference.py` if an Entrez ID didn't have a direct Ensembl match in GENCODE v24) into a proper `ENTREZ:<numeric_id>` format.
*   **Action:** Reads each `*_standardized_preprocessed.h5ad`, modifies the relevant `.var` columns, and overwrites the file (after creating a `.fixed` temporary file).

**Stage 2.7: Analyze ENCODE Gene Mapping Quality (Informational)**

*   **Action:** An inline Python script within `run_rnaseq_pipeline.sh` that loads the ENCODE preprocessed H5AD and prints statistics about the gene ID mapping quality (e.g., percentage of genes mapped, counts by mapping source).
*   **Output:** Console output and `encode_mapping_stats.json`.

### 4.3. Single-Cell RNA-seq Processing (GTEx Example - Conceptual)

*This is a planned extension and current scripts primarily focus on bulk data. The GTEx single-cell atlas is already available as an H5AD.*

*   **Input:** GTEx Single-Cell RNA-seq Atlas (e.g., `GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad`).
*   **Key Steps (Conceptual - to be implemented in e.g., `process_gtex_single_cell.py`):**
    1.  **Load Data:** Read the public GTEx snRNA-seq H5AD.
    2.  **Pseudobulking:**
        *   For each unique combination of `subject_id` and `cell_type` (e.g., from `public_obs.Granular_cell_type`), aggregate gene expression across all cells belonging to that group. This typically involves taking the sum or mean of raw counts (or normalized counts if appropriate).
        *   The resulting expression matrix will have rows representing `subject_id`_`cell_type` (or a multi-index) and columns representing genes.
    3.  **Metadata Integration:**
        *   Create an `.obs` DataFrame for the pseudobulk samples. Each row corresponds to a `subject_id`_`cell_type` combination.
        *   Populate with `subject_id`, `cell_type`, and other relevant metadata inherited from the original GTEx subject/sample attributes (e.g., `tissue`, `sex`, `age`, ethnicity fields). This will require merging with the metadata loaded by `load_gtex_metadata`.
    4.  **Gene ID Standardization:** Apply gene ID standardization to the pseudobulk expression matrix's gene dimension (columns) using the primary gene annotation reference (`gencode_v24_complete_mapping.csv`), similar to Stage 2.5 for bulk data. Populate `.var` accordingly.
    5.  **Output:** Save as `gtex_scrnaseq_pseudobulk_standardized.h5ad`.

### 4.4. Combined Dataset Creation (All Genes, Sparse)

*   **Script:** `create_combined_dataset_all_genes_sparse.py`
*   **Purpose:** To merge all processed bulk/pseudobulk datasets (`*_standardized_preprocessed.h5ad`) into a single, large AnnData object. This object will include all unique genes across all datasets, with missing genes for a particular sample/dataset represented as NaN or 0 in a sparse matrix.
*   **Key Actions:**
    1.  **Identify Datasets:** Finds all `*_standardized_preprocessed.h5ad` files.
    2.  **Load Reference Mapping:** Uses `gencode_v24_complete_mapping.csv`.
    3.  **Unified Gene Index:** Creates a comprehensive list of all unique gene IDs (standardized Ensembl IDs, ENTREZ IDs, spike-ins) present across all input datasets.
    4.  **Sparse Matrix Construction:**
        *   Iterates through each dataset.
        *   For each sample, maps its gene expression values to the unified gene index. If a gene from the unified index is not present in the current dataset for that sample, its value will implicitly be zero (or NaN if explicitly handled, though sparse matrices typically store only non-zero values).
        *   Builds a `scipy.sparse.csr_matrix` to store the combined expression data efficiently.
    5.  **Combine Metadata:** Concatenates `.obs` DataFrames from all datasets.
    6.  **Create Unified `.var`:** Creates a new `.var` DataFrame based on the unified gene index, populating it with annotations from the reference mapping and information about which datasets each gene is present in.
    7.  **Save Combined AnnData:** Saves the final combined sparse AnnData object.

### 4.5. WGS Data Processing and Linkage (Future Work)

*   **Current Status:** WGS data processing is handled by a separate pipeline.
*   **Goal for Integration:** To link `subject_id` from the RNA-seq AnnData objects to the corresponding sample IDs in the WGS VCF files.
*   **Planned Steps:**
    1.  **Establish Linking Key:** Confirm the exact `subject_id` format to be used for linkage.
    2.  **Create External Mapping File:** Generate a mapping file: `rna_seq_subject_id` -> `wgs_vcf_sample_id`. This may require manual curation or leveraging existing project metadata.
    3.  **RNA-seq AnnData Enhancement:** Add a `wgs_vcf_sample_id` column to the `.obs` of relevant RNA-seq AnnData objects.
    4.  **VCF Header Annotation (Separate Process):** Annotate WGS VCF file headers with `##META` lines indicating the paired RNA-seq sample ID(s) and subject ID, as per VCF specifications. This is outside the scope of the current RNA-seq pipeline scripts but is a related downstream step for full integration.

## 5. Validation

*   **Script:** `validate_standardized_datasets.py`
*   **Purpose:** To perform automated checks on the final standardized AnnData objects (`*_standardized_v2.h5ad` or `*_standardized_preprocessed.h5ad`).
*   **Checks Performed:**
    *   Presence and correctness of harmonized GENCODE version (`v24`) and reference genome (`hg38`) in `.uns`.
    *   Presence and completeness of core metadata fields in `.obs` (e.g., `tissue`, `sex`, `species`, `data_type`).
    *   Validity of ontology mappings (e.g., `tissue_ontology` starts with `UBERON:`).
    *   Format of gene identifiers in `.var_names` (expects Ensembl format `ENSG...`).
*   **Output:** A JSON and HTML validation report summarizing the status for each dataset.

## 6. Workspace Setup

*   **Script:** `setup_workspace.sh`
*   **Purpose:** To prepare the environment (primarily for Kubernetes or similar containerized setups) by installing necessary system dependencies and Python packages.
*   **Actions:**
    *   Installs system packages like `build-essential`, `gcc`, `python3-dev`.
    *   Installs Python packages listed in the script (numpy, pandas, scanpy, anndata, etc.) using pip.
    *   Verifies installations.

## 7. Known Limitations and Future Work

*   **MAGE Age Data:** Age data for MAGE donors is not currently integrated; requires sourcing from 1000 Genomes supplementary information if available.
*   **Cell Type Ontology for ENCODE/MAGE:** While cell line names are present, a formal `cell_type_ontology_term_id` (CL terms) column needs to be systematically added for ENCODE cell lines and MAGE LCLs by creating and using a `celltype_to_cl.json` mapping.
*   **ADNI "Worst Diagnosis":** Logic to derive a "worst diagnosis over time" for ADNI subjects from the longitudinal data in `.uns` and add it to `.obs` or a donor summary is pending.
*   **GTEx Single-Cell Pseudobulking:** The full processing pipeline for GTEx single-cell data to create pseudobulk H5ADs needs to be implemented and integrated.
*   **WGS Linkage:** Full implementation of WGS data linkage.
*   **Schema Validation Tooling:** Explore and integrate more formal AnnData schema validation tools if available/recommended by CZI.
---