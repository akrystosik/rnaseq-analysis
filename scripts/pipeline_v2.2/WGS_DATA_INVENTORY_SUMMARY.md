# WGS Data Inventory Summary

## Overview
This document provides a comprehensive inventory of all WGS sample and subject IDs available in the pipeline, covering ENCODE cell lines, GTEx subjects, and legacy ENCODE reference data.

## Total Data Summary
- **Total WGS samples**: 995
- **Primary datasets**: 3 (ENCODE, GTEx, ENCODE_legacy)
- **Data types**: Germline variants, somatic variants, alignment data

## Dataset Breakdown

### 1. ENCODE Cell Line WGS Data (9 samples)
**Primary ENCODE cell lines with complete WGS processing:**

| Cell Line | Sample ID | Data Types | File Types Available |
|-----------|-----------|------------|---------------------|
| A549 | A549_WGS_pair1 | germline+somatic | BAM, VCF_germline, VCF_somatic |
| Caki2 | Caki2_WGS_pair1 | germline+somatic | VCF_germline, VCF_somatic |
| GM23248 | GM23248_WGS_pair1 | germline+somatic | BAM, VCF_germline, VCF_somatic |
| HepG2 | HepG2_WGS_pair1 | germline+somatic | BAM, VCF_germline, VCF_somatic |
| K562 | K562_WGS_pair1 | germline+somatic | BAM, VCF_germline, VCF_somatic |
| NCI-H460 | NCI-H460_WGS_pair1 | germline+somatic | BAM, VCF_germline, VCF_somatic |
| Panc1 | Panc1_WGS_pair1 | germline+somatic | BAM, VCF_germline, VCF_somatic |
| T47D | T47D_WGS_pair1 | germline+somatic | BAM, VCF_germline, VCF_somatic |
| sknmc | sknmc_WGS_pair1 | germline+somatic | BAM, VCF_germline, VCF_somatic |

**Key characteristics:**
- All samples follow the naming pattern: `{CellLine}_WGS_pair1`
- 8 samples have complete BAM + VCF data, 1 sample (Caki2) has VCF-only
- Both germline (DeepVariant) and somatic (DeepSomatic) variants available
- All use GRCh38 reference genome

### 2. GTEx WGS Data (982 samples)
**Large-scale population-based WGS from GTEx project:**

| Metric | Value |
|--------|-------|
| Total samples | 982 |
| ID format | GTEX-{4-5 character suffix} |
| Data type | Germline variants only |
| File type | VCF_germline |
| Source VCF | GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf |

**ID pattern analysis:**
- Most common prefixes: 13* (107 samples), 14* (82 samples), 11* (58 samples)
- ID suffix lengths: 4-5 characters
- Examples: GTEX-16NPV, GTEX-1H1DG, GTEX-R3RS

**Sample coverage statistics:**
- Average variants per sample: ~3.6M total variants
- Average SNPs: ~3.3M per sample  
- Average Indels: ~260K per sample

### 3. ENCODE Legacy/Reference Data (4 samples)
**Historical ENCODE VCF files lifted from hg19 to hg38:**

| ENCODE File ID | Subject ID | Variants | Notes |
|----------------|------------|----------|-------|
| ENCFF752OAX | ENCODE_ref_ENCFF752OAX | 3,775,298 | Primary large variant set |
| ENCFF785JVR | ENCODE_ref_ENCFF785JVR | 3,544 | Structural variants |
| ENCFF574MDJ | ENCODE_ref_ENCFF574MDJ | 301 | Small variant set |
| ENCFF863MPP | ENCODE_ref_ENCFF863MPP | 296 | Small variant set |

**Key characteristics:**
- All lifted from hg19 to hg38 using liftOver
- Primarily contain indels and structural variants
- One large comprehensive set (ENCFF752OAX) with 3.8M variants

## Identifier Patterns and Conventions

### ENCODE Cell Lines
- **Sample naming**: `{CellLine}_WGS_pair1`
- **Subject naming**: `{CellLine}` (e.g., K562, HepG2)
- **Pair format**: All samples use `pair1` designation
- **Cell line types**: Cancer cell lines, immortalized cell lines

### GTEx Subjects  
- **ID format**: `GTEX-{suffix}`
- **Suffix patterns**: 4-5 alphanumeric characters
- **Subject = Sample**: GTEx uses same ID for subject and sample
- **Population**: Diverse human population from GTEx tissue donors

### Data File Organization
```
/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/
├── bam/                           # ENCODE alignment files
│   ├── A549/A549_WGS_pair1.marked_duplicates.bam
│   ├── K562/K562_WGS_pair1.marked_duplicates.bam
│   └── ...
├── variants/
│   ├── deepvariant/              # Germline variants (ENCODE)
│   ├── deepsomatic/              # Somatic variants (ENCODE)  
│   └── encode_lifted/            # Legacy hg38 lifted variants
└── qc/
    └── gtex_counts/              # GTEx variant summaries
```

## Data Types and Analysis Ready Status

### Germline Variants
- **ENCODE**: 9 samples (DeepVariant called)
- **GTEx**: 982 samples (joint-called population VCF)
- **ENCODE_legacy**: 4 reference samples
- **Total**: 995 samples

### Somatic Variants  
- **ENCODE**: 9 samples (DeepSomatic called)
- **Analysis**: Tumor vs normal comparison using cell line replicates

### Alignment Data
- **ENCODE**: 8 samples with BAM files
- **Quality**: Marked duplicates, indexed, analysis-ready

## Computational Processing Pipeline

### ENCODE Processing
1. **Alignment**: BWA-MEM to GRCh38
2. **Post-processing**: Picard MarkDuplicates  
3. **Germline calling**: DeepVariant
4. **Somatic calling**: DeepSomatic (tumor vs normal mode)
5. **QC**: FastQC, samtools stats, variant stats

### GTEx Processing
- **Source**: Pre-processed joint-called VCF from GTEx v9 release
- **Genome**: GRCh38/hg38
- **Processing**: Population-scale variant calling across 953 individuals

## Integration with RNA-seq Pipeline

The WGS data complements the RNA-seq standardization pipeline by providing:
- **Variant context**: Genetic variants affecting gene expression
- **Population structure**: Diverse genetic backgrounds in GTEx
- **Cell line variants**: Somatic alterations in ENCODE cell lines
- **Cross-dataset integration**: Shared sample identifiers enable multi-omic analysis

## File Inventory Output
Complete sample inventory saved to: `wgs_sample_inventory.csv`

**CSV Structure:**
- `sample_id`: Unique sample identifier
- `subject_id`: Subject/cell line identifier  
- `dataset`: Source dataset (ENCODE/GTEx/ENCODE_legacy)
- `data_type`: Variant types available (germline/somatic/germline+somatic)
- `file_types_available`: Data file types (BAM/VCF_germline/VCF_somatic/VCF_lifted_hg38)

This inventory provides the foundation for cross-dataset variant analysis and integration with the standardized RNA-seq expression data.