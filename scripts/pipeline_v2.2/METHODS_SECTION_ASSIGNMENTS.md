# Methods Section - Writing Assignments

## **Amy's Section: ENCODE Preprocessing (Draft)**

### ENCODE Data Processing and Cell Line Integration

ENCODE Project data comprised RNA-seq profiles from seven immortalized cell lines (A549, K562, HepG2, GM23248, Caki2, NCI-H460, Panc1) representing diverse tissue origins and cancer types. Raw expression data were originally quantified as transcripts per million (TPM) using RSEM v1.2.31 against GENCODE v29/GRCh38 annotations with Illumina sequencing platforms.

ENCODE-specific preprocessing addressed unique identifier challenges inherent to the ENCODE data structure. Gene identifiers were provided as numeric values requiring specialized mapping to standard Ensembl annotations. We implemented a tertiary mapping strategy utilizing ENCODE's internal gene annotation files to cross-reference numeric identifiers with their corresponding GENCODE gene symbols and Ensembl IDs.

Cell line metadata integration employed ENCODE's standardized biosample annotations. Cell type classifications utilized the Cell Ontology (CL) with mappings including "epithelial cell" (A549) → CL:0000066, "lymphoblast" (GM23248) → CL:0000542, and "hepatocyte" (HepG2) → CL:0000182. Tissue-of-origin annotations mapped to UBERON terms: "lung" (A549) → UBERON:0002048, "liver" (HepG2) → UBERON:0002107.

Quality control for ENCODE data included validation of gene mapping rates (target >99%), assessment of expression distribution consistency across cell lines, and verification of metadata completeness. Technical replicates were handled according to ENCODE consortium guidelines with appropriate sample aggregation strategies.

---

## **Sayan's Section: Data Loading and Integration (Draft)**

### Multi-Dataset Loading and Harmonization Infrastructure

Data loading employed a unified AnnData-based framework designed to accommodate heterogeneous data structures across the four source datasets. The loading infrastructure implemented dataset-specific parsers to handle varying file formats: HDF5 containers (GTEx), tab-delimited matrices (MAGE), CSV files with escaped delimiters (ADNI), and JSON metadata structures (ENCODE).

Sample identifier harmonization addressed complex naming conventions across datasets. We developed dataset-specific extraction algorithms: ADNI samples employed duplicate identifier removal (002_S_0413_002_S_0413_gencode_v24_pruned → 002_S_0413), MAGE utilized prefix extraction (NA06985_NA06985_gencode_V24_pruned → NA06985), and GTEx required subject-level aggregation (GTEX-1117F-0005-SM-HL9SH → GTEX-1117F).

Memory management strategies enabled processing of the full 21,004-sample dataset through sparse matrix representations and chunked loading procedures. The loading pipeline implemented compression algorithms (compressed sparse row format) achieving 85% memory reduction compared to dense matrix storage.

Cross-dataset validation during loading included gene overlap assessment (68,339 harmonized from 161,993 total), sample count verification, and metadata schema consistency checks. Error handling procedures captured loading failures with detailed logging for troubleshooting and data provenance tracking.

Integration quality control monitored successful loading rates (>99.5% target), memory utilization metrics, and processing time benchmarks. The unified data structure enabled seamless downstream analysis while preserving dataset-specific metadata and provenance information.

---

## **Combined Integration Points**

### Data Flow Architecture
1. **Raw Data Ingestion** (Sayan's section)
2. **Dataset-Specific Processing** (Amy's ENCODE section + similar for others)  
3. **Harmonization and Standardization** (Shared)
4. **Quality Control and Validation** (Shared)
5. **Output Generation** (Shared)

### Key Technical Achievements to Highlight
- **Scale**: 21,004 samples across 4 major datasets
- **Harmonization**: 68,339 genes from 161,993 total with >99% mapping success
- **Efficiency**: 85% memory reduction through sparse matrices
- **Reproducibility**: Fixed seeds, version control, comprehensive logging
- **Standards Compliance**: CZI Cell Census schema adherence

### Science Journal Style Guidelines
- **Concise but complete**: Essential technical details without excessive verbosity
- **Quantitative emphasis**: Specific performance metrics and success rates
- **Active voice**: "We implemented..." rather than passive constructions
- **Integration narrative**: Show how components work together systematically
- **Impact focus**: Emphasize novel aspects and technical achievements

### Writing Tips for Contributors
- **Amy**: Focus on ENCODE's unique challenges (numeric IDs, cell line diversity) and solutions
- **Sayan**: Emphasize infrastructure design, memory efficiency, and cross-dataset integration
- **Both**: Use specific numbers, version information, and performance metrics
- **Style**: Follow Science format with abbreviated technical descriptions and impact emphasis