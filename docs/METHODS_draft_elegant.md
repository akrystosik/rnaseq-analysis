# Methods

## Multi-Consortium RNA-seq Harmonization and Population Genomics Integration

### Dataset Compilation and Standardization

We assembled a comprehensive transcriptomic resource of 21,004 human samples by harmonizing data from four major genomics consortia. The integrated dataset comprised 19,616 GTEx tissue samples (93.4%) spanning 54 anatomical sites from 948 donors, 731 MAGE lymphoblastoid cell lines (3.5%) representing 26 global populations from the 1000 Genomes Project, 650 ADNI blood samples (3.1%) from Alzheimer's disease research participants, and 7 ENCODE immortalized cell lines (0.03%) modeling diverse cancer types.

Technical harmonization addressed platform heterogeneity across Illumina RNA-seq (GTEx, MAGE, ENCODE) and Affymetrix microarray (ADNI) technologies. All expression data were standardized to GENCODE v24 annotations on the GRCh38 reference genome through hierarchical gene identifier mapping. We implemented a four-tier mapping strategy: direct GENCODE matching, NCBI cross-reference resolution, dataset-specific custom mappings, and systematic placeholder assignment for unmapped genes. This approach achieved >99% mapping success across all datasets, creating a unified 68,339-gene expression matrix.

### Population Metadata Curation and Ancestry Integration

We implemented a rigorous distinction between self-reported ethnicity and genetically inferred ancestry, addressing a critical gap in population genomics research. Self-reported ethnicity data captured cultural and social identity using HANCESTRO ontology standardization, while genetic ancestry provided independent genomic evidence through integration with whole genome sequencing analysis.

For genetic ancestry inference, we applied dataset-specific methodologies optimized for data provenance. ADNI and GTEx samples received ancestry assignments from K-nearest neighbors classification on principal components of genome-wide SNP data (KNN_PCA method, mean confidence = 0.92). MAGE samples utilized their established 1000 Genomes population classifications (1000G_Classification method, confidence = 1.0) to avoid circular validation, as these samples served as reference populations for the KNN model.

Demographic metadata underwent systematic ontology mapping: tissue annotations employed UBERON anatomy terms, cell types used Cell Ontology classifications, developmental stages followed HsapDv ontology, and technical quality metrics preserved RNA integrity scores where available. This standardization enabled cross-consortium comparative analysis while maintaining dataset-specific biological signatures.

### Statistical Validation and Quality Assessment

We validated the ethnicity-ancestry distinction through concordance analysis on samples with authentic participant self-reports. Among 20,203 samples with both genetic ancestry predictions and self-reported ethnicity data, we observed 84.1% concordance (ADNI: 85.5%, GTEx: 83.6%), reflecting the scientifically appropriate distinction between cultural identity and genomic ancestry. MAGE samples were excluded from concordance analysis as their ethnicity labels derive from genetic classifications rather than participant reports.

Gene expression quality control employed standard metrics including mapping rates (>99% target), RNA integrity assessment (mean RIN = 8.2 for available samples), and cross-dataset consistency validation through Jaccard similarity analysis. Principal component analysis confirmed that biological factors (tissue type, ancestry, disease state) explained substantially more variance than technical factors (platform, batch effects), validating successful harmonization.

### Computational Implementation and Reproducibility

All analyses employed Python-based frameworks optimized for large-scale genomics: scanpy for expression data structures, pandas for metadata manipulation, and scikit-learn for ancestry classification. The pipeline generated AnnData-formatted outputs compatible with single-cell analysis ecosystems and CZI Cell Census schema v3.0.0. Processing utilized high-performance computing infrastructure with 64-256 GB RAM nodes, enabling memory-efficient sparse matrix operations on the 4.8 GB integrated dataset.

### Data Products and Dissemination

The standardization pipeline produced multiple research-ready outputs: individual consortium datasets with harmonized annotations, a unified multi-consortium expression atlas (21,004 samples × 68,339 genes), comprehensive ancestry-ethnicity mapping for 20,934 samples with dual metadata, and detailed validation reports documenting data quality and processing provenance. All outputs maintain FAIR data principles with persistent ontology identifiers and comprehensive metadata tracking.

---

## Extended Data Table 1. Dataset Integration Summary

| Consortium | Samples | Platform | Tissues/Conditions | Ancestry Coverage | Key Features |
|------------|---------|----------|-------------------|------------------|--------------|
| GTEx v10 | 19,616 | Illumina RNA-seq | 54 tissue sites | 99.7% (19,553) | Population-scale tissue atlas |
| MAGE | 731 | Illumina RNA-seq | Lymphoblastoid cells | 100% (731) | Global population diversity |
| ADNI | 650 | Affymetrix U219 | Blood samples | 100% (650) | Alzheimer's disease cohort |
| ENCODE | 7 | Illumina RNA-seq | Immortalized cell lines | 0% (0) | Cancer cell models |
| **Integrated** | **21,004** | **Multi-platform** | **Comprehensive** | **99.7% (20,934)** | **Unified transcriptomic resource** |

---

## Extended Data Table 2. Population Genomics Integration

| Ancestry | Samples | Proportion | Method | Ethnicity Concordance | Representative Populations |
|----------|---------|------------|--------|--------------------|---------------------------|
| European (EUR) | 17,127 | 81.5% | KNN_PCA + 1000G | 84.2% | CEU, GBR, FIN, IBS, TSI |
| African (AFR) | 2,522 | 12.0% | KNN_PCA + 1000G | 82.1% | YRI, LWK, GWD, MSL, ESN |
| Admixed American (AMR) | 689 | 3.3% | KNN_PCA + 1000G | 79.3% | MXL, PUR, CLM, PEL |
| East Asian (EAS) | 384 | 1.8% | KNN_PCA + 1000G | 88.7% | CHB, JPT, CHS, CDX, KHV |
| South Asian (SAS) | 212 | 1.0% | KNN_PCA + 1000G | 91.2% | GIH, PJL, BEB, STU, ITU |
| **Total** | **20,934** | **99.7%** | **Dual-method** | **84.1%** | **26 global populations** |

Note: Concordance reflects distinction between self-reported ethnicity (cultural identity) and genetic ancestry (genomic evidence).