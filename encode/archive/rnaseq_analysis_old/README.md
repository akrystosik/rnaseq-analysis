# RNA-seq Analysis Pipeline

## Directory Structure
- raw_data/: Original unprocessed data
  - fastq/: Raw FASTQ files
  - sra/: SRA format files (if downloading from NCBI)
- processed_data/: All processed data files
  - bam/: Aligned BAM files
  - counts/: Raw count matrices
  - normalized/: Normalized expression data
- scripts/: Analysis and processing scripts
  - download/: Scripts for data download
  - process/: Data processing scripts
  - qc/: Quality control scripts
  - analysis/: Analysis scripts
- qc/: Quality control outputs
  - fastqc/: FastQC reports
  - mapping_qc/: Alignment quality metrics
  - expression_qc/: Expression QC metrics
- metadata/: Sample and experiment metadata
  - encode/: ENCODE metadata
  - sample_sheets/: Sample information sheets
- docs/: Documentation and notes

## Pipeline Overview
This pipeline processes RNA-seq data from ENCODE K562 cell line experiments.
See the docs/ directory for detailed documentation.

## Usage
[Instructions for running the pipeline will be added here]

## Contact
Amy Krystosik
akrystosik@chanzuckerberg.com