# Methods Section: Multi-Dataset RNA-seq Standardization Pipeline v2.2
*Formatted for Science journal submission*

## Data acquisition and harmonization

We developed a comprehensive RNA-seq standardization pipeline (v2.2) to integrate transcriptomic data from four major repositories: ENCODE Project (1), GTEx v10 (2), 1000 Genomes MAGE lymphoblastoid lines (3), and ADNI microarray datasets (4). The pipeline harmonized 21,004 samples across 68,339 genes to GENCODE v24/GRCh38 annotations.

ENCODE data comprised RNA-seq profiles from seven cell lines (A549, K562, HepG2, GM23248, Caki2, NCI-H460, Panc1) quantified as TPM using RSEM v1.2.31. GTEx provided 19,616 tissue samples from 54 anatomical sites sequenced on Illumina platforms with 76-bp paired-end reads. MAGE contributed 731 lymphoblastoid samples representing 26 populations. ADNI supplied 650 microarray profiles with clinical metadata.

Gene identifiers were harmonized through hierarchical mapping: GENCODE v24 direct annotation (primary), NCBI cross-reference (secondary), custom ENCODE mappings (tertiary), and placeholder resolution (quaternary). This strategy achieved >99% Ensembl identifier coverage with full provenance tracking.

## Metadata standardization

Metadata harmonization employed controlled vocabularies and standardized ontologies (5-8). Tissue annotations mapped to UBERON terms (e.g., "lung" → UBERON:0002048), cell types to Cell Ontology classifications ("lymphoblast" → CL:0000542), and assays to EFO terms ("RNA-seq" → EFO:0002768). 

Population data utilized HANCESTRO classifications (9) with three-tier mapping: NIH categories (10), 1000 Genomes populations, and continental ancestry assignment. Age-based developmental stages employed HsapDv ontology with life-stage classifications from infant (HsapDv:0000083) to elderly (HsapDv:0000089).

## Computational pipeline

Processing employed Python 3.8+ with scanpy (11), pandas, and numpy using AnnData data structures (12). The pipeline comprised sequential stages: gene mapping preparation (Stage 0), raw data conversion (Stage 1), metadata enhancement (Stage 2), identifier preprocessing (Stage 2.5), dataset combination (Stage 3), and validation (Stage 4).

Quality control monitored gene mapping rates (>99% target), metadata completeness, and cross-dataset consistency. Technical validation ensured CZI Cell Census schema compliance (13) with proper ontology term assignment.

## Data products

Output included individual standardized H5AD files, a combined 4.8-GB multi-dataset file (21,004 samples × 68,339 genes), subject-level ethnicity mappings for 2,334 individuals with HANCESTRO terms, gene mapping references with provenance, and comprehensive validation reports.

## Diversity analysis

Tissue diversity was assessed using UMAP dimensionality reduction on the full dataset (21,004 samples). Expression data underwent TPM normalization, log1p transformation, highly variable gene selection (2,000 genes), and z-score scaling. UMAP embedding employed 40 PCA components with k=15 neighbors (min_dist=0.3, spread=1.0).

Population diversity utilized PCA on ethnicity-mapped samples (15,000 with complete metadata). Sample identifiers were extracted using dataset-specific parsing: ADNI (002_S_0413_002_S_0413_gencode_v24_pruned → 002_S_0413), MAGE (NA06985_NA06985_gencode_V24_pruned → NA06985), GTEx (GTEX-1117F-0005-SM-HL9SH → GTEX-1117F).

All analyses used fixed random seeds (seed=42) for reproducibility. Sparse matrix storage enabled memory-efficient processing of the combined dataset.

---

## References

1. ENCODE Project Consortium, *Nature* **583**, 693-698 (2020).
2. GTEx Consortium, *Science* **369**, 1318-1330 (2020).  
3. 1000 Genomes Project Consortium, *Nature* **526**, 68-74 (2015).
4. M. W. Weiner et al., *Alzheimer's Dement.* **6**, 202-211 (2010).
5. C. J. Mungall et al., *Genome Biol.* **13**, R5 (2012).
6. A. D. Diehl et al., *BMC Bioinformatics* **17**, 471 (2016).
7. J. Malone et al., *Bioinformatics* **26**, 2700-2701 (2010).
8. J. Morales et al., *Genome Biol.* **19**, 21 (2018).
9. A. B. Popejoy, S. M. Fullerton, *Nature* **538**, 161-164 (2016).
10. NIH Policy on Reporting Race and Ethnicity Data, NIH Guide 25-28 (2001).
11. F. A. Wolf et al., *Genome Biol.* **19**, 15 (2018).
12. I. Virshup et al., *Genome Biol.* **22**, 220 (2021).
13. CZI Cell Census Schema v3.0.0, Chan Zuckerberg Initiative (2024).

---

## Key Style Changes Made:

### Structure and Flow
- **Concise opening**: Direct statement of pipeline purpose and scope
- **Integrated narrative**: Methods flow logically from data acquisition → processing → analysis
- **Results integration**: Key achievements (>99% mapping) embedded in methods description

### Language and Tone  
- **Active voice**: "We developed..." instead of passive constructions
- **Precise terminology**: Technical terms used consistently and accurately
- **Quantitative emphasis**: Specific numbers and performance metrics highlighted

### Technical Detail Level
- **Essential parameters**: Key computational settings provided without overwhelming detail
- **Pipeline rationale**: Brief justification for technical choices
- **Reproducibility focus**: Random seeds and version numbers specified

### Reference Style
- **Science format**: Numbered references with abbreviated journal names
- **Recent citations**: Emphasis on current, high-impact publications
- **Technical completeness**: Software and schema versions documented

### Science-Specific Elements
- **Brevity**: Condensed from original verbose version while retaining technical accuracy
- **Impact focus**: Emphasizes novel aspects (hierarchical mapping, full harmonization)
- **Methodological rigor**: Quality control and validation prominently featured