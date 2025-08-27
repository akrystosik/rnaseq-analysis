# GTEx Dataset Clarification

## Two Dataset Scopes

When working with the GTEx data in this project, it's important to understand that there are two different scopes of data being referenced:

1. **Full GTEx Metadata**: 
   - Contains information about all 48,231 samples collected from 981 donors
   - Represents the complete GTEx collection across 70 tissue types
   - Referenced in metadata files: `GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt` and `GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt`

2. **Processed Expression Data**:
   - Contains only the 19,616 samples from 946 donors that were processed into the TPM matrix
   - Represents the subset of samples that passed quality control and were included in the final expression dataset
   - Referenced in the TPM file: `GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct`
   - This is the data that was converted to AnnData format for analysis

## Why the Difference?

The difference between these two datasets is due to quality control filtering and processing requirements in the GTEx project. Not all collected samples were ultimately included in the processed expression data due to factors such as:

- Sample quality issues
- RNA quality thresholds
- Sequencing depth requirements
- Other technical considerations

## In Our Analysis

Throughout our documentation and code:
- When we reference "samples in the metadata," we're referring to the full 48,231 samples.
- When we reference "the TPM matrix" or "processed data," we're referring to the 19,616 samples.
- Our tissue-specific analyses based on the AnnData files use only the processed data.

This distinction is important when interpreting sample counts, donor distributions, and other metrics throughout the project.