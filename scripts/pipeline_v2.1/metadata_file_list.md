# Metadata File Contents

The following metadata files were used or referenced by the pipeline, along with their contents:

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/encode_metadata.json`

```json
{
  "label": "ENCODE RNA-seq Expression Data (Processed)",
  "dataset_version": "Various ENCODE experiments",
  "version": "1.0",
  "gencode_version": "v29", 
  "harmonized_gencode_version": "v24",
  "reference_genome": "GRCh38", 
  "harmonized_reference_genome": "hg38",
  "doi": [
    "10.1038/s41586-020-2493-4" 
  ],
  "license_terms": {
    "license_url": "https://www.encodeproject.org/about/data-use-policy/",
    "download_date": "2024-11-01", 
    "dataset_owner": "ENCODE Project Consortium"
  },
  "data_type": "RNA-seq",
  "platform": "Illumina",
  "genome_info": {
    "genome_version": "GRCh38", 
    "confidence": "high",
    "notes": "ENCODE data standards typically use GRCh38 reference genome."
  },
  "rna_seq_protocol": {
    "protocol_type": "RNA-seq mixed",
    "protocol_confidence": "medium", 
    "selection_method": "Variable - includes polyA plus RNA-seq and total RNA-seq",
    "sequencing_platform": "Illumina (various models)",
    "read_configuration": "Typically paired-end",
    "strand_specificity": "Typically strand-specific (reverse); based on standard ENCODE processing pipelines." 
  },
  "assay_ontology": "EFO:0009922", 
  "obs_columns": { 
    "data_type": "RNA-seq",
    "assay_ontology": "EFO:0009922",
    "species": "human",
    "species_ontology": "NCBITaxon:9606",
    "expression_unit": "TPM"
  },
  "metadata_source": {
    "citation": "ENCODE Project Consortium. (2020). Perspectives on ENCODE. Nature, 583(7818), 699-710.",
    "url": "https://www.encodeproject.org/",
    "notes": "ENCODE project data, originally using GRCh38/GENCODE v29, harmonized to GENCODE v24/hg38 for this project. Includes encode cell lines."
  },
  "extras": [
    {"key": "data_product_type", "value": "h5ad"},
    {"key": "dataset_name_long", "value": "ENCODE RNA-seq Expression Data (Encode Cell Lines, Processed by CZI DNA2Cell Pipeline)"},
    {"key": "project_description_short", "value": "Processed ENCODE RNA-seq TPM data (cell lines), harmonized to GENCODE v24/hg38."},
    {"key": "assay_label", "value": "RNA-seq"},
    {"key": "organism_label", "value": "Homo sapiens"}
  ],
  "dataset_info": {
    "project": "ENCODE",
    "cell_lines": {
      "A549": {
        "description": "Human lung carcinoma cell line derived from a 58 year old European male with lung adenocarcinoma",
        "tissue": "lung",
        "tissue_ontology": "UBERON:0002048",
        "disease": "adenocarcinoma",
        "disease_ontology": "MONDO:0004992",
        "cell_type": "epithelial",
        "cell_ontology": "CL:0000066",
        "sex": "male",
        "sex_ontology": "PATO:0000384",
        "organism": "human",
        "organism_ontology": "NCBITaxon:9606",
        "age": "58",
        "age_units": "year",
        "development_stage": "adult",
        "development_stage_ontology": "HsapDv:0000087",
        "ethnicity": "European",
        "ethnicity_ontology": "HANCESTRO:0005",
        "subject_id": "ENCDO000AAZ",
        "external_ids": [
          "GEO:SAMN05733878",
          "Cellosaurus:CVCL_0023"
        ],
        "assay": "polyA plus RNA-seq",
        "assay_ontology": "EFO:0009922",
        "assay_term_id": "OBI:0002571",
        "nucleic_acid_type": "polyadenylated mRNA",
        "nucleic_acid_term_id": "SO:0000871",
        "depletion": "rRNA",
        "depleted_in_term_id": ["SO:0000252"],
        "size_range": ">200",
        "strand_specificity": "reverse",
        "experiment_accession": "ENCSR000CON",
        "experiment_description": "The libraries contained in this experiment come from the whole cell fraction of independent growths of cell line A549. They are stranded PE76 Illumina GAIIx RNA-Seq libraries from rRNA-depleted Poly-A+ RNA > 200 nucleotides in size.",
        "experiment_date_released": "2011-10-17",
        "replication_type": "isogenic",
        "bio_replicate_count": 2,
        "file": {
          "accession": "ENCFF244DNJ",
          "format": "tsv",
          "type": "gene quantifications",
          "assembly": "GRCh38",
          "genome_annotation": "V29",
          "download_url": "https://www.encodeproject.org/files/ENCFF244DNJ/@@download/ENCFF244DNJ.tsv",
          "biological_replicate": 2,
          "technical_replicate": "2_1",
          "date_created": "2021-04-05T23:46:45.413132+00:00"
        },
        "pipeline": {
          "version": "1.2.3",
          "analysis_id": "ENCAN400RTQ",
          "title": "Bulk RNA-seq",
          "accession": "ENCPL862USL",
          "software": "RSEM",
          "software_version": "1.2.31",
          "award_rfas": ["ENCODE4"]
        },
        "biosample_ontology": {
          "term_id": "EFO:0001086",
          "term_name": "A549",
          "classification": "cell line",
          "dbxrefs": ["Cellosaurus:CVCL_0023"],
          "organ_slims": ["lung"],
          "cell_slims": ["cancer cell"],
          "developmental_slims": ["endoderm"],
          "system_slims": ["respiratory system"],
          "synonyms": ["A549 cell", "A-549"]
        }
      },
      "K562": {
        "description": "Human chronic myelogenous leukemia cell line",
        "tissue": "bone marrow",
        "tissue_ontology": "UBERON:0002371",
        "disease": "chronic myelogenous leukemia",
        "disease_ontology": "MONDO:0005059",
        "cell_type": "lymphoblast",
        "cell_ontology": "CL:0000945",
        "sex": "female",
        "sex_ontology": "PATO:0000383",
        "organism": "human",
        "organism_ontology": "NCBITaxon:9606",
        "age": "53",
        "age_units": "year",
        "development_stage": "adult",
        "development_stage_ontology": "HsapDv:0000087",
        "ethnicity": "unknown",
        "subject_id": "ENCDO000AAD",
        "external_ids": [
          "GEO:SAMN04284550",
          "Cellosaurus:CVCL_0004"
        ],
        "assay": "RNA-seq",
        "assay_ontology": "EFO:0009922",
        "assay_term_id": "OBI:0001271",
        "nucleic_acid_type": "RNA",
        "nucleic_acid_term_id": "SO:0000356",
        "assay_category": "total RNA-seq",
        "depletion": "rRNA",
        "depleted_in_term_id": ["SO:0000252"],
        "size_range": ">200",
        "experiment_accession": "ENCSR000AEL",
        "replication_type": "isogenic",
        "bio_replicate_count": 2,
        "file": {
          "accession": "ENCFF171FQU",
          "format": "tsv",
          "type": "gene quantifications",
          "assembly": "GRCh38",
          "genome_annotation": "V29",
          "download_url": "https://www.encodeproject.org/files/ENCFF171FQU/@@download/ENCFF171FQU.tsv",
          "biological_replicate": 2,
          "technical_replicate": "2_1",
          "date_created": "2021-04-06T00:00:56.723705+00:00"
        },
        "pipeline": {
          "version": "1.2.3",
          "analysis_id": "ENCAN653KGT",
          "title": "Bulk RNA-seq",
          "accession": "ENCPL862USL",
          "software": "RSEM",
          "software_version": "1.2.31",
          "award_rfas": ["ENCODE4"]
        },
        "biosample_ontology": {
          "term_id": "EFO:0002067",
          "term_name": "K562",
          "classification": "cell line",
          "dbxrefs": ["Cellosaurus:CVCL_0004"],
          "organ_slims": ["blood", "bodily fluid"],
          "cell_slims": ["cancer cell", "leukocyte", "hematopoietic cell"],
          "developmental_slims": ["mesoderm"],
          "system_slims": ["immune system"],
          "synonyms": ["K-562", "K-562 cell", "K562 cell"]
        }
      },
      "HepG2": {
        "description": "Human liver hepatocellular carcinoma cell line",
        "tissue": "liver",
        "tissue_ontology": "UBERON:0002107",
        "disease": "hepatocellular carcinoma",
        "disease_ontology": "MONDO:0007256",
        "cell_type": "epithelial",
        "cell_ontology": "CL:0000066",
        "sex": "male",
        "sex_ontology": "PATO:0000384",
        "organism": "human",
        "organism_ontology": "NCBITaxon:9606",
        "age": "15",
        "age_units": "year",
        "development_stage": "juvenile",
        "development_stage_ontology": "HsapDv:0000082",
        "ethnicity": "European",
        "subject_id": "ENCDO000AAC",
        "external_ids": [
          "GEO:SAMN04284581",
          "Cellosaurus:CVCL_0027"
        ],
        "assay": "RNA-seq",
        "assay_ontology": "EFO:0009922",
        "assay_term_id": "OBI:0001271",
        "nucleic_acid_type": "RNA",
        "nucleic_acid_term_id": "SO:0000356",
        "assay_category": "total RNA-seq",
        "depletion": "rRNA",
        "depleted_in_term_id": ["SO:0000252"],
        "size_range": ">200",
        "experiment_accession": "ENCSR245ATJ",
        "file": {
          "accession": "ENCFF863QWG",
          "format": "tsv",
          "type": "gene quantifications",
          "assembly": "GRCh38",
          "genome_annotation": "V29",
          "size": 11067419,
          "md5sum": "37dee2d5f5a694c1fa014541b0129650",
          "download_url": "https://www.encodeproject.org/files/ENCFF863QWG/@@download/ENCFF863QWG.tsv",
          "s3_uri": "s3://encode-public/2020/10/30/432cbe0e-fa24-47bb-9e89-4b082b18066d/ENCFF863QWG.tsv",
          "azure_uri": "https://datasetencode.blob.core.windows.net/dataset/2020/10/30/432cbe0e-fa24-47bb-9e89-4b082b18066d/ENCFF863QWG.tsv?sv=2019-10-10&si=prod&sr=c&sig=9qSQZo4ggrCNpybBExU8SypuUZV33igI11xw0P7rB3c%3D",
          "biological_replicate": 1,
          "technical_replicate": "1_1",
          "date_created": "2020-10-30T18:10:21.474736+00:00",
          "derived_from": [
            "/files/ENCFF463TZA/", 
            "/files/ENCFF285DRD/"
          ]
        },
        "quality_metrics": {
          "number_of_genes_detected": 13086
        },
        "pipeline": {
          "version": "1.2.1",
          "analysis_id": "ENCAN711WSR",
          "title": "Bulk RNA-seq",
          "accession": "ENCPL862USL",
          "software": "RSEM",
          "software_version": "1.2.31",
          "award_rfas": ["ENCODE4"]
        },
        "biosample_ontology": {
          "term_id": "EFO:0001187",
          "term_name": "HepG2",
          "classification": "cell line",
          "dbxrefs": ["Cellosaurus:CVCL_0027"],
          "organ_slims": ["liver", "endocrine gland", "epithelium", "exocrine gland"],
          "cell_slims": ["cancer cell", "epithelial cell"],
          "developmental_slims": ["endoderm"],
          "system_slims": ["endocrine system", "digestive system", "exocrine system"],
          "synonyms": ["HepG2 cell", "Hep-G2"]
        }
      },
      "GM23248": {
        "description": "Human skin fibroblast cell line",
        "tissue": "skin",
        "tissue_ontology": "UBERON:0002097",
        "disease": "normal",
        "disease_ontology": "PATO:0000461",
        "cell_type": "fibroblast",
        "cell_ontology": "CL:0000057",
        "sex": "male",
        "sex_ontology": "PATO:0000384",
        "organism": "human",
        "organism_ontology": "NCBITaxon:9606",
        "age": "53",
        "age_units": "year",
        "development_stage": "adult",
        "development_stage_ontology": "HsapDv:0000087",
        "ethnicity": "European",
        "subject_id": "ENCDO336AAA",
        "external_ids": [
          "Cellosaurus:CVCL_F183"
        ],
        "assay": "RNA-seq",
        "assay_ontology": "EFO:0009922",
        "assay_term_id": "OBI:0001271",
        "nucleic_acid_type": "RNA",
        "nucleic_acid_term_id": "SO:0000356",
        "assay_category": "total RNA-seq",
        "depletion": "rRNA",
        "depleted_in_term_id": ["SO:0000252"],
        "size_range": ">200",
        "experiment_accession": "ENCSR797BPP",
        "file": {
          "accession": "ENCFF640FPG",
          "format": "tsv",
          "type": "gene quantifications",
          "assembly": "hg19",
          "genome_annotation": "V19",
          "size": 9517666,
          "md5sum": "e695262d29249120665fd22334902382",
          "download_url": "https://www.encodeproject.org/files/ENCFF640FPG/@@download/ENCFF640FPG.tsv",
          "s3_uri": "s3://encode-public/2015/08/24/c878e6b3-a834-4b72-999f-88dd0fee5708/ENCFF640FPG.tsv",
          "azure_uri": "https://datasetencode.blob.core.windows.net/dataset/2015/08/24/c878e6b3-a834-4b72-999f-88dd0fee5708/ENCFF640FPG.tsv?sv=2019-10-10&si=prod&sr=c&sig=9qSQZo4ggrCNpybBExU8SypuUZV33igI11xw0P7rB3c%3D",
          "biological_replicate": 2,
          "technical_replicate": "2_1",
          "date_created": "2015-08-24T20:14:29.003512+00:00",
          "derived_from": [
            "/files/ENCFF791UML/",
            "/files/ENCFF826ONU/"
          ]
        },
        "quality_metrics": {
          "Pearson_correlation": 0.9349629,
          "Spearman_correlation": 0.9384167,
          "MAD_of_log_ratios": 0.455,
          "SD_of_log_ratios": 0.718
        },
        "pipeline": {
          "version": "1.2.19",
          "analysis_id": "ENCAN562HZW",
          "title": "RNA-seq of long RNAs (paired-end, stranded)",
          "accession": "ENCPL002LPE",
          "software": "RSEM",
          "software_version": "1.2.19",
          "award_rfas": ["ENCODE3"]
        },
        "biosample_ontology": {
          "term_id": "EFO:0005723",
          "term_name": "GM23248",
          "classification": "cell line",
          "dbxrefs": ["Cellosaurus:CVCL_F183"],
          "organ_slims": ["limb", "connective tissue", "skin of body"],
          "cell_slims": ["fibroblast", "connective tissue cell"],
          "developmental_slims": ["ectoderm"],
          "system_slims": ["integumental system"],
          "synonyms": ["GM23248 cell"]
        }
      },
      "Caki2": {
        "description": "Human kidney clear cell carcinoma cell line",
        "tissue": "kidney",
        "tissue_ontology": "UBERON:0002113",
        "disease": "clear cell renal cell carcinoma",
        "disease_ontology": "MONDO:0005086",
        "cell_type": "epithelial",
        "cell_ontology": "CL:0000066",
        "sex": "male",
        "sex_ontology": "PATO:0000384",
        "organism": "human",
        "organism_ontology": "NCBITaxon:9606",
        "age": "69",
        "age_units": "year",
        "development_stage": "adult",
        "development_stage_ontology": "HsapDv:0000087",
        "ethnicity": "European",
        "subject_id": "ENCDO140IFG",
        "external_ids": [
          "Cellosaurus:CVCL_0235"
        ],
        "assay": "RNA-seq",
        "assay_ontology": "EFO:0009922",
        "assay_term_id": "OBI:0001271",
        "nucleic_acid_type": "RNA",
        "nucleic_acid_term_id": "SO:0000356",
        "assay_category": "total RNA-seq",
        "depletion": "rRNA",
        "depleted_in_term_id": ["SO:0000252"],
        "size_range": ">200",
        "experiment_accession": "ENCSR584JXD",
        "file": {
          "accession": "ENCFF685WJV",
          "format": "tsv",
          "type": "gene quantifications",
          "assembly": "GRCh38",
          "genome_annotation": "V29",
          "size": 10847316,
          "md5sum": "5eb62a214db1abc18d83efcbecb94816",
          "download_url": "https://www.encodeproject.org/files/ENCFF685WJV/@@download/ENCFF685WJV.tsv",
          "s3_uri": "s3://encode-public/2021/04/17/2a74df90-2900-4a4f-9b5f-3e14300443b6/ENCFF685WJV.tsv",
          "azure_uri": "https://datasetencode.blob.core.windows.net/dataset/2021/04/17/2a74df90-2900-4a4f-9b5f-3e14300443b6/ENCFF685WJV.tsv?sv=2019-10-10&si=prod&sr=c&sig=9qSQZo4ggrCNpybBExU8SypuUZV33igI11xw0P7rB3c%3D",
          "biological_replicate": 2,
          "technical_replicate": "2_1",
          "date_created": "2021-04-17T09:07:21.528692+00:00",
          "derived_from": [
            "/files/ENCFF071QWV/",
            "/files/ENCFF285DRD/"
          ]
        },
        "quality_metrics": {
          "number_of_genes_detected": 11699,
          "Pearson_correlation": 0.8834402,
          "Spearman_correlation": 0.8942174,
          "MAD_of_log_ratios": 0.66,
          "SD_of_log_ratios": 1.046
        },
        "pipeline": {
          "version": "1.2.3",
          "analysis_id": "ENCAN995LNB",
          "title": "Bulk RNA-seq",
          "accession": "ENCPL862USL",
          "software": "RSEM",
          "software_version": "1.2.31",
          "award_rfas": ["ENCODE4"]
        },
        "biosample_ontology": {
          "term_id": "EFO:0002150",
          "term_name": "Caki2",
          "classification": "cell line",
          "dbxrefs": ["Cellosaurus:CVCL_0235"],
          "organ_slims": ["kidney", "epithelium"],
          "cell_slims": ["cancer cell", "epithelial cell"],
          "developmental_slims": ["mesoderm"],
          "system_slims": ["excretory system"],
          "synonyms": ["Caki-2"]
        }
      },
      "NCI-H460": {
        "description": "Human lung large cell carcinoma cell line",
        "tissue": "lung",
        "tissue_ontology": "UBERON:0002048",
        "disease": "large cell lung carcinoma",
        "disease_ontology": "MONDO:0005092",
        "cell_type": "epithelial",
        "cell_ontology": "CL:0000066",
        "sex": "male",
        "sex_ontology": "PATO:0000384",
        "organism": "human",
        "organism_ontology": "NCBITaxon:9606",
        "age": "",
        "age_units": "",
        "development_stage": "adult",
        "development_stage_ontology": "HsapDv:0000087",
        "ethnicity": "unknown",
        "subject_id": "ENCDO267PTG",
        "external_ids": [
          "Cellosaurus:CVCL_0459"
        ],
        "assay": "RNA-seq",
        "assay_ontology": "EFO:0009922",
        "assay_term_id": "OBI:0001271",
        "nucleic_acid_type": "RNA",
        "nucleic_acid_term_id": "SO:0000356",
        "assay_category": "total RNA-seq",
        "depletion": "rRNA",
        "depleted_in_term_id": ["SO:0000252"],
        "size_range": ">200",
        "experiment_accession": "ENCSR164OCT",
        "file": {
          "accession": "ENCFF876SRX",
          "format": "tsv",
          "type": "gene quantifications",
          "assembly": "GRCh38",
          "genome_annotation": "V29",
          "size": 10853793,
          "md5sum": "010ca18712a353e8eb5179dad8ede71a",
          "download_url": "https://www.encodeproject.org/files/ENCFF876SRX/@@download/ENCFF876SRX.tsv",
          "s3_uri": "s3://encode-public/2021/04/17/257eda93-2e5f-464c-a7df-1eecb2e93f70/ENCFF876SRX.tsv",
          "azure_uri": "https://datasetencode.blob.core.windows.net/dataset/2021/04/17/257eda93-2e5f-464c-a7df-1eecb2e93f70/ENCFF876SRX.tsv?sv=2019-10-10&si=prod&sr=c&sig=9qSQZo4ggrCNpybBExU8SypuUZV33igI11xw0P7rB3c%3D",
          "biological_replicate": 1,
          "technical_replicate": "1_1",
          "date_created": "2021-04-17T09:17:41.303339+00:00",
          "derived_from": [
            "/files/ENCFF169YJS/",
            "/files/ENCFF285DRD/"
          ]
        },
        "quality_metrics": {
          "number_of_genes_detected": 11941,
          "Pearson_correlation": 0.984119,
          "Spearman_correlation": 0.9843655,
          "MAD_of_log_ratios": 0.298,
          "SD_of_log_ratios": 0.377
        },
        "pipeline": {
          "version": "1.2.3",
          "analysis_id": "ENCAN187YBG",
          "title": "Bulk RNA-seq",
          "accession": "ENCPL862USL",
          "software": "RSEM",
          "software_version": "1.2.31",
          "award_rfas": ["ENCODE4"]
        },
        "biosample_ontology": {
          "term_id": "EFO:0003044",
          "term_name": "NCI-H460",
          "classification": "cell line",
          "dbxrefs": ["Cellosaurus:CVCL_0459"],
          "organ_slims": ["epithelium", "lung"],
          "cell_slims": ["cancer cell", "epithelial cell"],
          "developmental_slims": ["endoderm"],
          "system_slims": ["respiratory system"],
          "synonyms": []
        }
      },
      "Panc1": {
        "description": "Human pancreatic ductal carcinoma cell line",
        "tissue": "pancreas",
        "tissue_ontology": "UBERON:0001264",
        "disease": "pancreatic carcinoma",
        "disease_ontology": "MONDO:0007256",
        "cell_type": "ductal cell",
        "cell_ontology": "CL:0000068",
        "sex": "male",
        "sex_ontology": "PATO:0000384",
        "organism": "human",
        "organism_ontology": "NCBITaxon:9606",
        "age": "56",
        "age_units": "year",
        "development_stage": "adult",
        "development_stage_ontology": "HsapDv:0000087",
        "ethnicity": "European",
        "subject_id": "ENCDO000ABB",
        "external_ids": [
          "Cellosaurus:CVCL_0480"
        ],
        "assay": "RNA-seq",
        "assay_ontology": "EFO:0009922",
        "assay_term_id": "OBI:0001271",
        "nucleic_acid_type": "RNA",
        "nucleic_acid_term_id": "SO:0000356",
        "assay_category": "total RNA-seq",
        "depletion": "rRNA",
        "depleted_in_term_id": ["SO:0000252"],
        "size_range": ">200",
        "experiment_accession": "ENCSR128CYL",
        "file": {
          "accession": "ENCFF710IFD",
          "format": "tsv",
          "type": "gene quantifications",
          "assembly": "GRCh38",
          "genome_annotation": "V29",
          "size": 11064756,
          "md5sum": "72c97e2051b2ca81840d0f28b5ffbc9d",
          "download_url": "https://www.encodeproject.org/files/ENCFF710IFD/@@download/ENCFF710IFD.tsv",
          "s3_uri": "s3://encode-public/2020/10/30/6aa93aeb-7080-4f73-8e14-7bf4215c24b0/ENCFF710IFD.tsv",
          "azure_uri": "https://datasetencode.blob.core.windows.net/dataset/2020/10/30/6aa93aeb-7080-4f73-8e14-7bf4215c24b0/ENCFF710IFD.tsv?sv=2019-10-10&si=prod&sr=c&sig=9qSQZo4ggrCNpybBExU8SypuUZV33igI11xw0P7rB3c%3D",
          "biological_replicate": 1,
          "technical_replicate": "1_1",
          "date_created": "2020-10-30T17:48:29.648357+00:00",
          "derived_from": [
            "/files/ENCFF285DRD/",
            "/files/ENCFF227BTT/"
          ]
        },
        "quality_metrics": {
          "number_of_genes_detected": 13407
        },
        "pipeline": {
          "version": "1.2.1",
          "analysis_id": "ENCAN061KZF",
          "title": "Bulk RNA-seq",
          "accession": "ENCPL862USL",
          "software": "RSEM",
          "software_version": "1.2.31",
          "award_rfas": ["ENCODE4"]
        },
        "biosample_ontology": {
          "term_id": "EFO:0002713",
          "term_name": "Panc1",
          "classification": "cell line",
          "dbxrefs": ["Cellosaurus:CVCL_0480"],
          "organ_slims": ["pancreas", "epithelium"],
          "cell_slims": ["cancer cell", "epithelial cell"],
          "developmental_slims": ["endoderm"],
          "system_slims": ["endocrine system"],
          "synonyms": ["PANC-1"]
        }
      }
    },
    "data_repository": [
      {
        "name": "ENCODE Project",
        "url": "https://www.encodeproject.org/"
      }
    ]
  }
}
```

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gtex_metadata.json`

```json
{
  "label": "GTEx v10 RNA-seq Expression Data (Processed)",
  "dataset_version": "v10",
  "version": "1.0", 
  "gencode_version": "39",
  "harmonized_gencode_version": "v24",
  "reference_genome": "hg38",
  "harmonized_reference_genome": "hg38",
  "doi": [
    "10.1126/science.aaz1776",
    "10.1038/nature24277"
  ],
  "license_terms": {
    "license_url": "https://gtexportal.org/home/license",
    "download_date": "2025-03-10", 
    "dataset_owner": "GTEx Consortium"
  },
  "rna_seq_protocol": {
    "protocol_type": "RNA-seq polyA+",
    "protocol_confidence": "high",
    "selection_method": "polyA+ selection using Illumina TruSeq protocol",
    "sequencing_platform": "Illumina HiSeq 2000/2500",
    "read_configuration": "76-bp paired-end reads",
    "strand_specificity": "Non-strand specific"
  },
  "obs_columns": {
    "data_type": "RNA-seq",
    "assay_ontology": "EFO:0009922",
    "species": "human",
    "species_ontology": "NCBITaxon:9606",
    "expression_unit": "TPM" 
  },
  "metadata_source": {
    "citation": "GTEx Consortium (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science.",
    "url": "https://www.science.org/doi/10.1126/science.aaz1776",
    "notes": "GTEx v10 data, originally using GENCODE v39, harmonized to GENCODE v24 and hg38 for this project. RNA-seq libraries were prepared using the Illumina TruSeq protocol.",
    "v10_notes": "V10 analysis summary Gencode version: v39 RNA-seq alignment: STAR v2.7.10a Transcript quantification: RSEM v1.3.3 Genotype used for eQTL analysis: WGS",
    "v10_source": "https://gtexportal.org/home/tissueSummaryPage"
  },

  "race_source_column": "RACE", 
  "source_race_to_standard_label_map": {
    "1": "asian",
    "2": "black or african american",
    "3": "white",
    "4": "american indian or alaska native",
    "5": "native hawaiian or other pacific islander",
    "98": "unknown or not reported",
    "99": "unknown or not reported"
  },
  "ethnicity_source_column": "ETHNCTY", 
  "source_ethnicity_to_standard_label_map": {
    "1": "hispanic or latino",
    "2": "not hispanic or latino",
    "98": "unknown or not reported",
    "99": "unknown or not reported"
  },


  "extras": [
    {"key": "data_product_type", "value": "h5ad"},
    {"key": "dataset_name_long", "value": "GTEx v10 RNA-seq Expression Data (Processed by CZI DNA2Cell Pipeline)"},
    {"key": "project_description_short", "value": "Processed GTEx v10 RNA-seq TPM data, harmonized to GENCODE v24/hg38."},
    {"key": "assay_label", "value": "RNA-seq"},
    {"key": "organism_label", "value": "Homo sapiens"}
  ]
}
```

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/mage_metadata.json`

```json
{
  "label": "MAGE (1000 Genomes LCL) RNA-seq Expression Data (Processed)",
  "dataset_version": "PRJNA851328", 
  "version": "1.0",                
  "gencode_version": "v38",        
  "harmonized_gencode_version": "v24",
  "reference_genome": "GRCh38",   
  "harmonized_reference_genome": "hg38",
  "doi": [
    "10.1038/s41588-023-01346-7" 
  ],
  "license_terms": {
    "license_url": "Data derived from 1000 Genomes Project (open access) and available via SRA/Zenodo. Check specific terms if any apply to derived data.",
    "download_date": "2025-01-01", 
    "dataset_owner": "Taylor DJ, et al. / 1000 Genomes Project Consortium" 
  },
  "data_type": "RNA-seq", 
  "platform": "Illumina",   
  "genome_info": {
    "genome_version": "GRCh38",
    "confidence": "high",
    "notes": "Reference genome set from publication data."
  },
  "rna_seq_protocol": {
    "protocol_type": "RNA-seq total", 
    "protocol_confidence": "medium", 
    "selection_method": "Lymphoblastoid cell lines (LCLs)",
    "sequencing_platform": "Illumina",
    "read_configuration": "paired-end",
    "strand_specificity": "Not specified" 
  },
  "assay_ontology": "EFO:0009922", 
  "obs_columns": { 
    "data_type": "RNA-seq", 
    "assay_ontology": "EFO:0009922",
    "species": "human",
    "species_ontology": "NCBITaxon:9606",
    "expression_unit": "TPM",
    "tissue": "lymphoblast" 
  },
  "metadata_source": {
    "citation": "Taylor DJ, Chhetri SB, Tassia MG, et al. Sources of gene expression variation in a globally diverse human cohort. Nat Genet 55, 750–761 (2023).",
    "url": "https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA851328",
    "notes": "RNA-seq data from 1000 Genomes Project LCLs, harmonized to GENCODE v24/hg38."
  },

  "dataset_info": { 
    "accession": "PRJNA851328",
    "sample_size": 731,
    "library_count": 779,
    "populations": 26,
    "continental_groups": 5,
    "cell_type": "lymphoblastoid cell lines", 
    "data_repository": [
      {"name": "Sequence Read Archive (SRA)", "accession": "PRJNA851328"},
      {"name": "Zenodo", "doi": "10.5281/zenodo.10535719"}
    ]
  },
  "cell_type_info": {
    "anatomical_entity_id": "UBERON:0000178", 
    "cell_type_id": "CL:0000542" 
  },
    
  "population_source_column_in_ped": "Population",
  "pop_to_race_map": {
    "ACB": "black or african american", "ASW": "black or african american", "ESN": "black or african american", "GWD": "black or african american", "LWK": "black or african american", "MSL": "black or african american", "YRI": "black or african american",
    "CLM": "multiethnic", "MXL": "multiethnic", "PEL": "multiethnic", "PUR": "multiethnic",
    "CDX": "asian", "CHB": "asian", "CHS": "asian", "JPT": "asian", "KHV": "asian",
    "CEU": "white", "FIN": "white", "GBR": "white", "IBS": "white", "TSI": "white",
    "BEB": "asian", "GIH": "asian", "ITU": "asian", "PJL": "asian", "STU": "asian"
  },
  "pop_to_hispanic_map": {
    "CLM": "hispanic or latino", "MXL": "hispanic or latino", "PEL": "hispanic or latino", "PUR": "hispanic or latino"
  },
  "default_race_for_unmapped_pop": "unknown or not reported",
  
  "extras": [
    {"key": "data_product_type", "value": "h5ad"},
    {"key": "dataset_name_long", "value": "MAGE (1000 Genomes LCL) RNA-seq Expression Data (Processed by CZI DNA2Cell Pipeline)"},
    {"key": "project_description_short", "value": "Processed MAGE RNA-seq TPM data from LCLs, harmonized to GENCODE v24/hg38."},
    {"key": "assay_label", "value": "RNA-seq"},
    {"key": "organism_label", "value": "Homo sapiens"}
  ]
}
```

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json`

```json
{
  "Adipose - Subcutaneous": "UBERON:0002190",
  "Adipose - Visceral (Omentum)": "UBERON:0016529",
  "Adrenal": "UBERON:0002369",
  "Adrenal Gland": "UBERON:0002369",
  "Artery - Aorta": "UBERON:0000947",
  "Artery - Coronary": "UBERON:0001621",
  "Artery - Tibial": "UBERON:0007610",
  "Bladder": "UBERON:0001255",
  "Blood": "UBERON:0000178",
  "Bone Marrow": "UBERON:0002371",
  "Brain": "UBERON:0000955",
  "Brain - Amygdala": "UBERON:0001876",
  "Brain - Anterior cingulate cortex (BA24)": "UBERON:0009835",
  "Brain - Caudate (basal ganglia)": "UBERON:0001873",
  "Brain - Cerebellar Hemisphere": "UBERON:0002245",
  "Brain - Cerebellum": "UBERON:0002037",
  "Brain - Cortex": "UBERON:0000956",
  "Brain - Frontal Cortex (BA9)": "UBERON:0013529",
  "Brain - Hippocampus": "UBERON:0002310",
  "Brain - Hypothalamus": "UBERON:0001898",
  "Brain - Nucleus accumbens (basal ganglia)": "UBERON:0001882",
  "Brain - Putamen (basal ganglia)": "UBERON:0001874",
  "Brain - Spinal cord (cervical c-1)": "UBERON:0002726",
  "Brain - Substantia nigra": "UBERON:0002038",
  "Breast - Mammary Tissue": "UBERON:0001911",
  "Cells - Cultured fibroblasts": "CL:0000057",
  "Cells - EBV-transformed lymphocytes": "CL:0000542",
  "Cervix - Ectocervix": "UBERON:0012249",
  "Cervix - Endocervix": "UBERON:0000458",
  "Colon": "UBERON:0001155",
  "Colon - Sigmoid": "UBERON:0001159",
  "Colon - Transverse": "UBERON:0001157",
  "Esophagus - Gastroesophageal Junction": "UBERON:0007650",
  "Esophagus - Mucosa": "UBERON:0002469",
  "Esophagus - Muscularis": "UBERON:0004648",
  "Fallopian Tube": "UBERON:0003889",
  "Heart - Atrial Appendage": "UBERON:0006631",
  "Heart - Left Ventricle": "UBERON:0002084",
  "Kidney": "UBERON:0002113",
  "Kidney - Cortex": "UBERON:0001225",
  "Kidney - Medulla": "UBERON:0000362",
  "Liver": "UBERON:0002107",
  "Lung": "UBERON:0002048",
  "Lymphoblast": "CL:0000542",
  "Minor Salivary Gland": "UBERON:0001830",
  "Muscle": "UBERON:0002385",
  "Muscle - Skeletal": "UBERON:0001134",
  "Nerve - Tibial": "UBERON:0001323",
  "Ovary": "UBERON:0000992",
  "PBMC": "UBERON:0000178",
  "Pancreas": "UBERON:0001264",
  "Peripheral blood": "UBERON:0000178",
  "Pituitary": "UBERON:0000007",
  "Prostate": "UBERON:0002367",
  "Skin": "UBERON:0002097",
  "Skin - Not Sun Exposed (Suprapubic)": "UBERON:0036151",
  "Skin - Sun Exposed (Lower leg)": "UBERON:0036149",
  "Small Intestine - Terminal Ileum": "UBERON:0001211",
  "Spleen": "UBERON:0002106",
  "Stomach": "UBERON:0000945",
  "Testis": "UBERON:0000473",
  "Thyroid": "UBERON:0002046",
  "Uterus": "UBERON:0000995",
  "Vagina": "UBERON:0000996",
  "Whole Blood": "UBERON:0000178",
  "adipose - subcutaneous": "UBERON:0002190",
  "adipose - visceral (omentum)": "UBERON:0016529",
  "adrenal": "UBERON:0002369",
  "adrenal gland": "UBERON:0002369",
  "artery - aorta": "UBERON:0000947",
  "artery - coronary": "UBERON:0001621",
  "artery - tibial": "UBERON:0007610",
  "bladder": "UBERON:0001255",
  "blood": "UBERON:0000178",
  "bone marrow": "UBERON:0002371",
  "brain": "UBERON:0000955",
  "brain - amygdala": "UBERON:0001876",
  "brain - anterior cingulate cortex (ba24)": "UBERON:0009835",
  "brain - caudate (basal ganglia)": "UBERON:0001873",
  "brain - cerebellar hemisphere": "UBERON:0002245",
  "brain - cerebellum": "UBERON:0002037",
  "brain - cortex": "UBERON:0000956",
  "brain - frontal cortex (ba9)": "UBERON:0013529",
  "brain - hippocampus": "UBERON:0002310",
  "brain - hypothalamus": "UBERON:0001898",
  "brain - nucleus accumbens (basal ganglia)": "UBERON:0001882",
  "brain - putamen (basal ganglia)": "UBERON:0001874",
  "brain - spinal cord (cervical c-1)": "UBERON:0002726",
  "brain - substantia nigra": "UBERON:0002038",
  "breast - mammary tissue": "UBERON:0001911",
  "cells - cultured fibroblasts": "CL:0000057",
  "cells - ebv-transformed lymphocytes": "CL:0000542",
  "cervix - ectocervix": "UBERON:0012249",
  "cervix - endocervix": "UBERON:0000458",
  "colon": "UBERON:0001155",
  "colon - sigmoid": "UBERON:0001159",
  "colon - transverse": "UBERON:0001157",
  "esophagus - gastroesophageal junction": "UBERON:0007650",
  "esophagus - mucosa": "UBERON:0002469",
  "esophagus mucosa": "UBERON:0002469",
  "esophagus - muscularis": "UBERON:0004648",
  "esophagus muscularis": "UBERON:0004648",
  "fallopian tube": "UBERON:0003889",
  "heart - atrial appendage": "UBERON:0006631",
  "heart - left ventricle": "UBERON:0002084",
  "kidney": "UBERON:0002113",
  "kidney - cortex": "UBERON:0001225",
  "kidney - medulla": "UBERON:0000362",
  "liver": "UBERON:0002107",
  "lung": "UBERON:0002048",
  "lymphoblast": "CL:0000542",
  "minor salivary gland": "UBERON:0001830",
  "muscle": "UBERON:0002385",
  "muscle - skeletal": "UBERON:0001134",
  "nerve - tibial": "UBERON:0001323",
  "ovary": "UBERON:0000992",
  "pancreas": "UBERON:0001264",
  "pbmc": "UBERON:0000178",
  "peripheral blood": "UBERON:0000178",
  "pituitary": "UBERON:0000007",
  "prostate": "UBERON:0002367",
  "skin": "UBERON:0002097",
  "skin - not sun exposed (suprapubic)": "UBERON:0036151",
  "skin - sun exposed (lower leg)": "UBERON:0036149",
  "small intestine - terminal ileum": "UBERON:0001211",
  "spleen": "UBERON:0002106",
  "stomach": "UBERON:0000945",
  "testis": "UBERON:0000473",
  "thyroid": "UBERON:0002046",
  "uterus": "UBERON:0000995",
  "vagina": "UBERON:0000996",
  "whole blood": "UBERON:0000178"
}
```

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/assay_to_efo.json`

```json
{
    "RNA-seq": "EFO:0002768",
    "total RNA-seq": "EFO:0008896",
    "polyA RNA-seq": "EFO:0008771",
    "microarray": "EFO:0002772",
    "CAGE": "EFO:0007045",
    "ChIP-seq": "EFO:0002692",
    "ATAC-seq": "EFO:0007045",
    "DNase-seq": "EFO:0004428",
    "methylation array": "EFO:0002759",
    "single-cell RNA-seq": "EFO:0010550",
    "bulk RNA-seq": "EFO:0008896",
    "polyA_plus": "EFO:0008771",
    "polyA_minus": "EFO:0008896",
    "total": "EFO:0008896"
  }
```

## `/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/age_to_hsapdv.json`

```json

  {
    "infant": "HsapDv:0000083",
    "child": "HsapDv:0000084",
    "juvenile": "HsapDv:0000085",
    "adolescent": "HsapDv:0000086",
    "adult": "HsapDv:0000087",
    "elderly": "HsapDv:0000088",
    "0-2": "HsapDv:0000083",
    "3-11": "HsapDv:0000084",
    "12-14": "HsapDv:0000085",
    "15-19": "HsapDv:0000086",
    "20-65": "HsapDv:0000087",
    "66+": "HsapDv:0000088",
    "fetal": "HsapDv:0000081",
    "newborn": "HsapDv:0000082",
    "prenatal": "HsapDv:0000081",
    "15": "HsapDv:0000086",
    "53": "HsapDv:0000087",
    "56": "HsapDv:0000087",
    "58": "HsapDv:0000087",
    "69": "HsapDv:0000087"
  }
```

