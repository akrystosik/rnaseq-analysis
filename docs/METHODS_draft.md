# Methods

## Multi-Consortium RNA-seq Harmonization and Metadata Curation

### Data Acquisition and Harmonization

We developed a comprehensive RNA-seq standardization pipeline to integrate transcriptomic data from four major repositories: ENCODE Project (1), GTEx v10 (2), 1000 Genomes MAGE lymphoblastoid lines (3), and ADNI microarray datasets (4). The pipeline harmonized 21,004 samples across 68,339 genes to GENCODE v24/GRCh38 annotations.

ENCODE data comprised RNA-seq profiles from seven cell lines (A549, K562, HepG2, GM23248, Caki2, NCI-H460, Panc1) quantified as TPM using RSEM v1.2.31. GTEx provided 19,616 tissue samples from 54 anatomical sites sequenced on Illumina platforms with 76-bp paired-end reads. MAGE contributed 731 lymphoblastoid samples representing 26 populations from the 1000 Genomes Project. ADNI supplied 650 microarray profiles from the Affymetrix Human Genome U219 Array with comprehensive clinical metadata.

Gene identifiers were harmonized through hierarchical mapping: GENCODE v24 direct annotation (primary) and NCBI cross-reference resolution (secondary). This strategy achieved >99% Ensembl identifier coverage with full provenance tracking for all mapping decisions.

### Metadata Standardization Schema

Metadata harmonization employed controlled vocabularies and standardized ontologies (5-8) to enable cross-consortium integration. Tissue annotations mapped to UBERON anatomy terms (e.g., "lung" → UBERON:0002048), cell types to Cell Ontology classifications ("lymphoblast" → CL:0000542), and assays to Experimental Factor Ontology terms ("RNA-seq" → EFO:0002768).

Population metadata implemented a rigorous distinction between self-reported ethnicity and genetically inferred ancestry. Self-reported ethnicity utilized HANCESTRO classifications (9) with three-tier mapping from NIH demographic categories through 1000 Genomes population codes to continental ancestry classifications. Genetic ancestry data were integrated from complementary whole genome sequencing analysis to provide independent genomic evidence, with dataset-specific methodologies: KNN classification on principal components for ADNI/GTEx samples, and established 1000 Genomes population classifications for MAGE samples to avoid circular validation.

Age-based developmental stages employed Human Developmental Stages (HsapDv) ontology with life-stage classifications from infant (HsapDv:0000083) to elderly (HsapDv:0000089). Technical quality metrics preserved RNA integrity numbers where available, and clinical metadata maintained longitudinal diagnostic information for ADNI samples.

### Computational Pipeline Implementation

Processing employed Python 3.8+ with scanpy (11), pandas, and NumPy using AnnData data structures (12). The pipeline comprised sequential stages: gene mapping preparation, raw data standardization, metadata enhancement, identifier preprocessing, dataset combination, and comprehensive validation.

Quality control monitored gene mapping rates (>99% target), metadata completeness across all standardized fields, and cross-dataset consistency through Jaccard similarity analysis. Technical validation ensured CZI Cell Census schema v3.0.0 compliance (13) with proper ontology term assignment and categorical data formatting.

Sample identifiers were standardized using dataset-specific parsing protocols to enable multi-omic linkage: ADNI samples (002_S_0413_002_S_0413_gencode_v24_pruned → 002_S_0413), MAGE samples (NA06985_NA06985_gencode_V24_pruned → NA06985), and GTEx samples (GTEX-1117F-0005-SM-HL9SH → GTEX-1117F).

### Data Products and Validation

Output products included individual standardized H5AD files with comprehensive metadata annotations, a combined 4.8-GB multi-consortium dataset (21,004 samples × 68,339 genes), ancestry-ethnicity mapping for samples with dual metadata, gene mapping references with full provenance tracking, and detailed validation reports documenting quality control metrics.

Diversity analysis employed UMAP dimensionality reduction on the complete integrated dataset. Expression data underwent TPM normalization, log1p transformation, highly variable gene selection (2,000 genes), and z-score scaling. UMAP embedding utilized 40 PCA components with k=15 neighbors (min_dist=0.3, spread=1.0) for optimal biological structure visualization.

All analyses used fixed random seeds (seed=42) for reproducibility. Sparse matrix storage enabled memory-efficient processing of the combined dataset while maintaining full sample and gene annotations.

---

## Extended Data Table 1. Metadata Curation Standards

| Field | Standard/Goal | Implementation | 
|-------|---------------|----------------|
| Donor ID (`subject_id`) | Consistent across RNA-seq & WGS for multi-omic linkage | Dataset-specific parsing protocols |
| Tissue (`tissue`) | UBERON ontology terms for anatomical standardization | Controlled vocabulary mapping |
| Cell Type (`cell_type`) | Cell Ontology terms aligned with CellXGene schema | Systematic ontology integration |
| Assay (`assay_ontology`) | Experimental Factor Ontology classifications | Platform-specific term assignment |
| Age (`age`) | Standardized format with developmental stage mapping | HsapDv ontology integration |
| Self-Reported Ethnicity (`self_reported_ethnicity`) | HANCESTRO terms for cultural/social identity | Three-tier hierarchical mapping |
| Genetic Ancestry (`inferred_ancestry`) | Continental classifications with confidence scoring | WGS integration with dual methodology |
| Species (`species`) | NCBI Taxon ID standardization | Universal human classification |
| Sex (`sex`) | CellXGene-aligned categories | Standardized terminology |
| Gene IDs (`gene_id`) | GENCODE v24 Ensembl identifiers | Hierarchical mapping strategy |
| Reference Genome | GRCh38/hg38 harmonization | Coordinate system standardization |

---

## References

1. ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. *Nature* **583**, 699–710 (2020).

2. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* **369**, 1318–1330 (2020).

3. 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* **526**, 68–74 (2015).

4. Weiner, M. W. et al. The Alzheimer's Disease Neuroimaging Initiative: progress report and future plans. *Alzheimers Dement.* **6**, 202–211 (2010).

5. Mungall, C. J. et al. Uberon, an integrative multi-species anatomy ontology. *Genome Biol.* **13**, R5 (2012).

6. Diehl, A. D. et al. The Cell Ontology 2016: enhanced content, modularization, and ontology interoperability. *J. Biomed. Semantics* **7**, 44 (2016).

7. Malone, J. et al. Modeling sample variables with an Experimental Factor Ontology. *Bioinformatics* **26**, 1112–1118 (2010).

8. Morales, J. et al. A standardized framework for representation of ancestry data in genomics studies, with application to the NHGRI-EBI GWAS Catalog. *Genome Biol.* **19**, 21 (2018).

9. Popejoy, A. B. & Fullerton, S. M. Genomics is failing on diversity. *Nature* **538**, 161–164 (2016).

10. NIH Policy and Guidelines on The Inclusion of Women and Minorities as Subjects in Clinical Research. NIH Guide **23**, 11 (1994).

11. Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* **19**, 15 (2018).

12. Virshup, I. et al. The scverse project provides a computational ecosystem for single-cell omics data analysis. *Nat. Biotechnol.* **41**, 604–606 (2023).

13. Chanzit, T. et al. Single-cell RNA sequencing technologies and bioinformatics pipelines. *Exp. Mol. Med.* **50**, 96 (2018).