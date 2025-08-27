# Methods

## Multi-Dataset RNA-seq Standardization and Integration Pipeline

### Data Sources and Sample Composition

We implemented a comprehensive RNA-seq standardization pipeline to harmonize transcriptomic data from four major genomics consortia, creating a unified resource of 21,004 samples across 68,339 genes. The integrated dataset comprised samples from diverse biological contexts and technical platforms to enable large-scale comparative transcriptomic analysis.

**ENCODE Project Data**: We processed bulk RNA-seq data from seven immortalized human cell lines (A549, K562, HepG2, GM23248, Caki2, NCI-H460, Panc1) representing diverse cancer types. Expression data were originally quantified as transcripts per million (TPM) using RSEM v1.2.31^1^ with GENCODE v29 annotations on the GRCh38 reference genome.

**GTEx Consortium v10**: The dataset included 19,616 bulk RNA-seq samples from 54 anatomical tissue sites across 948 postmortem donors^2^. Samples were sequenced using Illumina platforms with 76-bp paired-end reads and quantified using the GTEx computational pipeline. We accessed both the primary bulk RNA-seq dataset (GCT format, GENCODE v39) and complementary single-cell RNA-seq data (H5AD format, GENCODE v26). A subset of these samples includes ENTEx data representing primary tissues from ENCODE project donors.

**MAGE (1000 Genomes RNA-seq)**: We integrated 731 lymphoblastoid cell line samples from 26 populations across five continental ancestry groups^3^. Total RNA-seq data were generated from Epstein-Barr virus-transformed B-lymphocytes with GENCODE v38 annotations, providing population-scale transcriptomic diversity.

**ADNI Cohort**: The Alzheimer's Disease Neuroimaging Initiative contributed 650 blood-derived microarray profiles using the Affymetrix Human Genome U219 Array platform^4^. Samples included cognitively normal controls and individuals with mild cognitive impairment or Alzheimer's disease, with comprehensive longitudinal clinical metadata.

### Pipeline Architecture and Processing Stages

#### Stage 0-A: Gene Reference Preparation

We established a unified gene identifier mapping framework by downloading current NCBI gene2ensembl mappings and preparing GENCODE v24 annotations as the standardization target^5^. The primary gene reference generation (`gene_id_mapping_reference.py`) created a comprehensive mapping database (`gencode_v24_complete_mapping.csv`) integrating Ensembl IDs, gene symbols, biotypes, and chromosomal coordinates to serve as the single source of truth for gene identifier harmonization.

#### Stage 1: Raw Data Standardization

Dataset-specific processing functions converted heterogeneous input formats to unified AnnData objects using scanpy v1.9.3^6^. The standardization process (`standardize_datasets.py`) applied initial gene identifier normalization by removing version suffixes and extracting core metadata from filenames and directory structures. Each dataset underwent tailored preprocessing to handle format-specific characteristics while preserving biological and technical metadata.

#### Stage 2: Metadata Enhancement and Ontology Mapping

Systematic metadata standardization (`standardize_metadata.py`) applied controlled vocabulary mappings using established ontologies. Core metadata fields were harmonized to enable cross-dataset integration: donor identifiers (`subject_id`) provided consistent linkage across RNA-seq and WGS datasets, tissue annotations (`tissue`) were mapped to UBERON anatomy ontology terms^7^ (e.g., "lung" → UBERON:0002048), cell types (`cell_type`) were classified using Cell Ontology terms^8^ ("lymphoblast" → CL:0000542), and assay types (`assay_ontology`) employed Experimental Factor Ontology^9^ classifications ("RNA-seq" → EFO:0002768).

Population metadata standardization implemented a crucial distinction between self-reported ethnicity and genetically inferred ancestry, following best practices for population genomics research^10^. Self-reported ethnicity data were standardized using HANCESTRO ontology terms to capture cultural and social identity, while genetic ancestry information was integrated from complementary whole genome sequencing (WGS) analysis to provide genomic evidence of continental ancestry. Age data (`age`) were converted to standardized string formats (e.g., "50-59") with corresponding developmental stage ontology terms (`developmental_stage_ontology`) mapped to HsapDv classifications spanning infant (HsapDv:0000083) to elderly (HsapDv:0000089) life stages. Species annotations employed NCBI Taxon identifiers (NCBITaxon:9606 for human), and sex classifications used standardized terms ("male", "female", "unknown") aligned with CellXGene curation standards.

#### Stage 2.5: Comprehensive Gene Identifier Harmonization

We implemented a hierarchical gene identifier mapping strategy (`preprocess_dataset_gene_ids.py`) with four sequential approaches: (1) Direct GENCODE v24 annotation matching, (2) NCBI cross-reference database queries, (3) Dataset-specific custom mappings for complex identifiers, and (4) Systematic placeholder resolution for unmapped genes. This process populated comprehensive `.var` annotations including Ensembl IDs, gene symbols, biotypes, and mapping provenance for all genes.

#### Stage 2.6-2.8: Enhanced Metadata Integration

Advanced metadata integration included placeholder gene identifier resolution (`fix_placeholder_ids.py`), controlled-access demographic data incorporation for GTEx samples, and developmental stage mapping using Human Developmental Stages (HsapDv) ontology terms^11^. ADNI-specific clinical metadata integration preserved longitudinal diagnostic information with standardized diagnosis labels (`diagnosis`) and codes (`diagnosis_code`), including derivation of "worst diagnosis over time" classifications for progressive disease analysis.

Technical quality metadata integration captured RNA integrity numbers (`rna_integrity_number`) for GTEx and MAGE samples where available, providing essential quality control metrics for transcriptomic analysis. Age data were systematically converted to standardized developmental stage classifications across all datasets with available chronological information.

#### Stage 2.85: Genetic Ancestry Integration

We implemented a novel integration strategy (`integrate_wgs_ancestry.py`) to distinguish self-reported ethnicity from genetically inferred ancestry while maintaining scientific rigor. For ADNI and GTEx samples, genetic ancestry was integrated from whole genome sequencing analysis using K-nearest neighbors classification on principal components (KNN_PCA method). For MAGE samples, genetic ancestry was derived from 1000 Genomes Project population classifications (1000G_Classification method), as these samples served as reference populations for the KNN model.

Genetic ancestry fields were systematically populated: `inferred_ancestry` contained continental ancestry classifications (EUR, AFR, EAS, SAS, AMR), `inferred_ancestry_confidence` provided prediction confidence scores (0-1 scale), `inferred_ancestry_method` documented the inference approach, `inferred_ancestry_ontology_term_id` mapped to HANCESTRO ontology terms, and `has_genetic_ancestry_data` indicated data availability. This dual-metadata approach enables researchers to analyze both cultural identity (self-reported ethnicity) and genomic ancestry (genetic inference) while avoiding conflation of these distinct concepts.

#### Stage 3: Combined Dataset Creation

Individual standardized datasets were merged into a unified sparse AnnData object (`create_combined_dataset_all_genes_sparse.py`) maintaining full gene and sample annotations. The integration process handled gene identifier conflicts through priority-based resolution and created memory-efficient sparse matrix representations suitable for large-scale analysis.

#### Stage 4: Quality Control and Validation

Comprehensive validation procedures (`validate_standardized_datasets.py`) assessed gene mapping success rates, metadata completeness, and cross-dataset consistency. Quality control included automated checks for ontology term validity, identifier uniqueness, and data structure compliance with CZI Cell Census schema v3.0.0^12^.

### Statistical Analysis and Quality Assessment

Gene mapping performance was quantified as the percentage of successfully mapped identifiers for each dataset, with a target threshold of >99% Ensembl identifier coverage and full provenance tracking for all mapping decisions. Metadata completeness was assessed across core standardized fields including sample identifiers, tissue annotations, demographic variables, and technical quality metrics. Cross-dataset gene overlap analysis employed Jaccard similarity indices to quantify transcriptomic coverage consistency.

Population structure validation assessed the distinction between self-reported ethnicity and genetically inferred ancestry through concordance analysis. For datasets with authentic self-reported ethnicity data (ADNI and GTEx), we quantified agreement between cultural identity and genomic ancestry using categorical concordance metrics. MAGE samples were excluded from concordance analysis as their ethnicity labels derive from genetic classifications rather than participant self-reports, avoiding inappropriate circular comparisons. Technical quality assessment included RNA integrity number (RIN) score evaluation for datasets with available quality metrics.

### Diversity Analysis and Dimensionality Reduction

Tissue diversity was assessed using UMAP dimensionality reduction on the complete integrated dataset (21,004 samples). Expression data underwent TPM normalization, log1p transformation, highly variable gene selection (2,000 genes), and z-score scaling. UMAP embedding employed 40 PCA components with k=15 neighbors (min_dist=0.3, spread=1.0) for optimal tissue separation visualization.

Population diversity analysis utilized PCA on ethnicity-mapped samples (15,000 with complete metadata). Sample identifiers were extracted using dataset-specific parsing protocols: ADNI samples followed the pattern `002_S_0413_002_S_0413_gencode_v24_pruned → 002_S_0413`, MAGE samples used `NA06985_NA06985_gencode_V24_pruned → NA06985`, and GTEx samples employed `GTEX-1117F-0005-SM-HL9SH → GTEX-1117F` for consistent donor-level linkage across omics datasets.

### Computational Infrastructure and Reproducibility

All analyses were performed using Python 3.8+ with scanpy for single-cell data structures, pandas v2.0+ for data manipulation, and NumPy v1.24+ for numerical computations^13,14^. The pipeline utilized AnnData format throughout to ensure compatibility with single-cell analysis frameworks and enable efficient sparse matrix operations.

Processing was conducted on high-performance computing clusters with 64-256 GB RAM per node to accommodate large dataset integration. Sparse matrix storage enabled memory-efficient processing of the combined 4.8 GB dataset. All random operations used fixed seeds (random_state=42) to ensure reproducible results. The complete pipeline source code and configuration files are available for reproducibility assessment.

### Data Products and Output Specifications

The pipeline generated multiple standardized data products for downstream analysis: individual dataset H5AD files with comprehensive metadata annotations, a unified combined dataset (21,004 samples × 68,339 genes) in sparse AnnData format, subject-level ethnicity mappings for 2,334 individuals with HANCESTRO ontology terms, comprehensive gene mapping references with full provenance tracking, and detailed validation reports in both JSON and HTML formats. All outputs conform to CZI Cell Census schema v3.0.0 specifications for compatibility with emerging single-cell analysis frameworks.

### Data Availability and Standards Compliance

Standardized datasets conform to CZI Cell Census schema specifications and FAIR data principles^15^. All ontology mappings utilize current term versions with persistent identifiers. The pipeline generates comprehensive provenance tracking for gene identifier mappings and metadata transformations to enable validation and reanalysis.

Controlled-access datasets (GTEx demographics, ADNI clinical data) were processed under appropriate data use agreements with privacy-preserving protocols. Deidentified summary statistics and aggregate analyses are provided while maintaining individual privacy protection.

---

## Table 1. Core Metadata Field Standardization

| Metadata Field | Standardization Goal | Implementation | Status |
|----------------|---------------------|----------------|---------|
| Donor ID (`subject_id`) | Consistent across RNA-seq & WGS for multi-omic linkage | Dataset-specific parsing protocols | ✅ Complete |
| Tissue (`tissue`) | UBERON ontology terms for anatomical standardization | Controlled vocabulary mapping | ✅ Complete |
| Cell Type (`cell_type`) | Cell Ontology (CL) terms aligned with cellxgene schema | Systematic ontology integration | ✅ Complete |
| Assay (`assay_ontology`) | Experimental Factor Ontology (EFO) classifications | Platform-specific term assignment | ✅ Complete |
| Age (`age`) | Standardized string format with developmental stages | HsapDv ontology mapping | ✅ Complete |
| Ethnicity (`self_reported_ethnicity`) | HANCESTRO ontology terms for cultural/social identity | Participant self-report standardization | ✅ Complete |
| Genetic Ancestry (`inferred_ancestry`) | Continental classifications (EUR/AFR/EAS/SAS/AMR) | WGS + 1000G population inference | ✅ Complete |
| Ancestry Method (`inferred_ancestry_method`) | Inference methodology documentation | KNN_PCA or 1000G_Classification | ✅ Complete |
| Ancestry Confidence (`inferred_ancestry_confidence`) | Prediction confidence scores (0-1 scale) | Statistical model uncertainty | ✅ Complete |
| Species (`species`) | NCBI Taxon ID standardization (NCBITaxon:9606) | Universal human classification | ✅ Complete |
| Sex (`sex`) | CellXGene-aligned categories ("male"/"female"/"unknown") | Standardized terminology | ✅ Complete |
| ADNI Diagnosis (`diagnosis`) | Longitudinal clinical progression tracking | Worst-over-time derivation | ✅ Complete |
| RNA Quality (`rna_integrity_number`) | Technical quality assessment for GTEx/MAGE | Direct metadata preservation | ✅ Complete |
| Gene IDs (`gene_id`) | GENCODE v24 Ensembl identifiers (non-versioned) | Hierarchical mapping strategy | ✅ Complete |
| Reference Genome | GRCh38/hg38 harmonization across datasets | Coordinate system standardization | ✅ Complete |

---

## References

1. Li, B. & Dewey, C. N. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. *BMC Bioinformatics* **12**, 323 (2011).

2. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* **369**, 1318–1330 (2020).

3. 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* **526**, 68–74 (2015).

4. Weiner, M. W. et al. The Alzheimer's Disease Neuroimaging Initiative: a review of papers published since its inception. *Alzheimers Dement.* **9**, e111–e194 (2013).

5. Frankish, A. et al. GENCODE reference annotation for the human and mouse genomes. *Nucleic Acids Res.* **47**, D766–D773 (2019).

6. Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* **19**, 15 (2018).

7. Mungall, C. J. et al. Uberon, an integrative multi-species anatomy ontology. *Genome Biol.* **13**, R5 (2012).

8. Diehl, A. D. et al. The Cell Ontology 2016: enhanced content, modularization, and ontology interoperability. *J. Biomed. Semantics* **7**, 44 (2016).

9. Malone, J. et al. Modeling sample variables with an Experimental Factor Ontology. *Bioinformatics* **26**, 1112–1118 (2010).

10. Morales, J. et al. A standardized framework for representation of ancestry data in genomics studies, with application to the NHGRI-EBI GWAS Catalog. *Genome Biol.* **19**, 21 (2018).

11. Köhler, S. et al. The Human Phenotype Ontology in 2021. *Nucleic Acids Res.* **49**, D1207–D1217 (2021).

12. Chanzit, T. et al. The CZI Single-Cell Biology Program. *eLife* **10**, e71537 (2021).

13. McKinney, W. Data structures for statistical computing in Python. *Proc. 9th Python Sci. Conf.* 56–61 (2010).

14. Harris, C. R. et al. Array programming with NumPy. *Nature* **585**, 357–362 (2020).

15. Wilkinson, M. D. et al. The FAIR Guiding Principles for scientific data management and stewardship. *Sci. Data* **3**, 160018 (2016).