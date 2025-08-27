# Methods

## Multi-Consortium RNA-seq Harmonization and Population Genomics Integration

### Dataset Compilation and Standardization

We assembled a comprehensive transcriptomic resource of 21,004 human samples by harmonizing data from four major genomics consortia. The integrated dataset comprised 19,616 GTEx tissue samples (93.4%) spanning 54 anatomical sites from 948 donors, 731 MAGE lymphoblastoid cell lines (3.5%) representing 26 global populations from the 1000 Genomes Project, 650 ADNI blood samples (3.1%) from Alzheimer's disease research participants, and 7 ENCODE immortalized cell lines (0.03%) modeling diverse cancer types.

Technical harmonization addressed platform heterogeneity across Illumina RNA-seq (GTEx, MAGE, ENCODE) and Affymetrix microarray (ADNI) technologies. All expression data were standardized to GENCODE v24 annotations on the GRCh38 reference genome through hierarchical gene identifier mapping. We implemented a four-tier mapping strategy: direct GENCODE matching, NCBI cross-reference resolution, dataset-specific custom mappings, and systematic placeholder assignment for unmapped genes. This approach achieved >99% mapping success across all datasets, creating a unified 68,339-gene expression matrix.

### Population Metadata Curation and Ancestry Integration

We implemented a rigorous distinction between self-reported ethnicity and genetically inferred ancestry, addressing a critical gap in population genomics research. Self-reported ethnicity data captured cultural and social identity using HANCESTRO ontology standardization, while genetic ancestry provided independent genomic evidence through integration with whole genome sequencing analysis.

For genetic ancestry inference, we applied dataset-specific methodologies optimized for data provenance. ADNI and GTEx samples received ancestry assignments from K-nearest neighbors classification on principal components of genome-wide SNP data (KNN_PCA method). MAGE samples utilized their established 1000 Genomes population classifications (1000G_Classification method) to avoid circular validation, as these samples served as reference populations for the KNN model.

Demographic metadata underwent systematic ontology mapping: tissue annotations employed UBERON anatomy terms, cell types used Cell Ontology classifications, developmental stages followed HsapDv ontology, and technical quality metrics preserved RNA integrity scores where available. This standardization enabled cross-consortium comparative analysis while maintaining dataset-specific biological signatures.

### Statistical Validation and Quality Assessment

We validated the ethnicity-ancestry distinction through concordance analysis on samples with authentic participant self-reports. Concordance analysis demonstrated scientifically appropriate agreement between cultural identity and genomic ancestry, with MAGE samples excluded as their ethnicity labels derive from genetic classifications rather than participant reports.

Gene expression quality control employed standard metrics including mapping rates, RNA integrity assessment, and cross-dataset consistency validation through Jaccard similarity analysis. Principal component analysis confirmed that biological factors (tissue type, ancestry, disease state) explained substantially more variance than technical factors (platform, batch effects), validating successful harmonization.

### Computational Implementation and Reproducibility

All analyses employed Python-based frameworks optimized for large-scale genomics: scanpy for expression data structures, pandas for metadata manipulation, and scikit-learn for ancestry classification. The pipeline generated AnnData-formatted outputs compatible with single-cell analysis ecosystems and CZI Cell Census schema v3.0.0. Processing utilized high-performance computing infrastructure with 64-256 GB RAM nodes, enabling memory-efficient sparse matrix operations on the 4.8 GB integrated dataset.

### Data Products and Dissemination

The standardization pipeline produced multiple research-ready outputs: individual consortium datasets with harmonized annotations, a unified multi-consortium expression atlas, comprehensive ancestry-ethnicity mapping with dual metadata, and detailed validation reports documenting data quality and processing provenance. All outputs maintain FAIR data principles with persistent ontology identifiers and comprehensive metadata tracking.

