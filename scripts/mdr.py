#!/usr/bin/env python3
#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/mdr.py
"""
Enhanced RNA-seq Metadata Collection Script

This script collects, standardizes, and stores metadata from multiple RNA-seq datasets,
integrating information from standardization logs and metadata JSON files for MDR registration.
Includes improved handling of ontology mappings and cross-dataset harmonization.
"""

import os
import re
import json
import logging
import argparse
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path
import scanpy as sc
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('metadata_collection')

# Define dataset names
DATASETS = ['gtex', 'encode', 'entex', 'mage', 'adni', 'combined']

class MetadataCollector:
    def __init__(self, data_dir, log_file, output_dir, mapping_dir=None):
        """
        Initialize the metadata collector
        
        Args:
            data_dir: Directory containing h5ad files
            log_file: Path to standardization log file
            output_dir: Directory to save standardized metadata
            mapping_dir: Directory containing ontology mapping files
        """
        self.data_dir = Path(data_dir)
        self.log_file = Path(log_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Create mapping directory if specified
        self.mapping_dir = None
        if mapping_dir:
            self.mapping_dir = Path(mapping_dir)
            self.mapping_dir.mkdir(exist_ok=True, parents=True)
        
        # Template for standardized metadata
        self.metadata_template = {
            "dataset_info": {
                "dataset_id": "",
                "dataset_name": "",
                "version": "",
                "publication_date": "",
                "update_date": datetime.now().strftime("%Y-%m-%d"),
                "doi": [],
                "contributors": [],
                "description": "",
                "steward": "",
                "license": {
                    "license_url": "",
                    "download_date": "",
                    "dataset_owner": ""
                }
            },
            "dataset_statistics": {
                "sample_count": 0,
                "donor_count": 0,
                "gene_count": 0,
                "tissue_count": 0
            },
            "demographics": {
                "age_distribution": {},
                "sex_distribution": {},
                "ethnicity_distribution": []
            },
            "sample_attributes": {
                "subject_id_field": "subject_id",
                "sex_field": "sex",
                "age_field": "age",
                "tissue_field": "tissue",
                "disease_field": ""
            },
            "ontology_mappings": {
                "assay": {
                    "ontology": "",
                    "label": ""
                },
                "organism": {
                    "ontology": "NCBITaxon:9606",
                    "label": "Homo sapiens"
                },
                "development_stage": {
                    "ontology": "",
                    "label": ""
                },
                "tissues": []
            },
            "genomic_metadata": {
                "reference_genome": "",
                "genome_source": "",
                "original_gencode_version": "",
                "harmonized_gencode_version": "",
                "gencode_mapping_notes": "",
                "expression_unit": ""
            },
            "rna_seq_metadata": {
                "protocol_type": "",
                "protocol_confidence": "",
                "selection_method": "",
                "sequencing_platform": "",
                "read_configuration": "",
                "strand_specificity": ""
            },
            "data_completeness": {
                "subject_id": 0,
                "sex": 0,
                "age": 0,
                "tissue": 0,
                "species": 0,
                "tissue_ontology": 0,
                "assay_ontology": 0,
                "developmental_stage_ontology": 0
            },
            "file_locations": {
                "internal": [],
                "external": []
            },
            "related_models_benchmarks": {
                "models": [],
                "benchmarks": []
            },
            "dataset_specific": {}
        }
        
        # Parse log file once
        self.log_data = self._parse_log_file()
        
        # Load or create tissue-to-ontology mapping table
        self.tissue_to_ontology = self._initialize_tissue_mapping()
        
    def _initialize_tissue_mapping(self):
        """Initialize or load tissue-to-ontology mapping table"""
        mapping_file = None
        if self.mapping_dir:
            mapping_file = self.mapping_dir / "tissue_to_uberon.csv"
            
        if mapping_file and mapping_file.exists():
            # Load existing mapping table
            try:
                df = pd.read_csv(mapping_file)
                mapping = dict(zip(df['tissue_name'], df['ontology_id']))
                logger.info(f"Loaded {len(mapping)} tissue-to-ontology mappings from {mapping_file}")
                return mapping
            except Exception as e:
                logger.error(f"Error loading tissue mapping: {e}")
        
        # Create default mapping for GTEx tissues
        # This is a starter set of common GTEx tissues mapped to UBERON
        mapping = {
            'Whole Blood': 'UBERON:0000178',
            'Muscle - Skeletal': 'UBERON:0001134',
            'Thyroid': 'UBERON:0002046',
            'Skin - Sun Exposed (Lower leg)': 'UBERON:0036149',
            'Lung': 'UBERON:0002048',
            'Brain - Cortex': 'UBERON:0000956',
            'Heart - Left Ventricle': 'UBERON:0001103',
            'Liver': 'UBERON:0002107',
            'Adipose - Subcutaneous': 'UBERON:0002190',
            'Breast - Mammary Tissue': 'UBERON:0001911'
            # Add more mappings as they become available
        }
        
        # Save default mapping if mapping directory exists
        if mapping_file:
            try:
                df = pd.DataFrame({
                    'tissue_name': list(mapping.keys()),
                    'ontology_id': list(mapping.values()),
                    'confidence': ['high'] * len(mapping),
                    'source': ['default'] * len(mapping)
                })
                df.to_csv(mapping_file, index=False)
                logger.info(f"Saved default tissue mapping to {mapping_file}")
            except Exception as e:
                logger.error(f"Error saving tissue mapping: {e}")
        
        return mapping
    
    def _parse_log_file(self):
        """Parse the standardization log file to extract metadata information"""
        logger.info(f"Parsing log file: {self.log_file}")
        
        if not self.log_file.exists():
            logger.error(f"Log file not found: {self.log_file}")
            return {}
        
        log_text = self.log_file.read_text()
        datasets_info = {}
        
        # Process each dataset section in the log
        for dataset_name in DATASETS:
            # Create case-insensitive pattern for dataset name
            pattern = re.compile(f"Metadata standardization complete for {dataset_name}", re.IGNORECASE)
            match = pattern.search(log_text)
            
            if not match:
                # Try alternative formatting
                upper_dataset = dataset_name.upper()
                pattern = re.compile(f"Metadata standardization complete for {upper_dataset}", re.IGNORECASE)
                match = pattern.search(log_text)
                
            if match:
                # Get the section of the log for this dataset
                start_idx = match.start()
                next_section = re.search(r"Metadata standardization complete for", log_text[start_idx + 1:])
                
                if next_section:
                    end_idx = start_idx + 1 + next_section.start()
                    section_text = log_text[start_idx:end_idx]
                else:
                    section_text = log_text[start_idx:]
                
                # Extract key information
                dataset_info = {
                    "sample_count": self._extract_value(section_text, r"Samples: (\d+)"),
                    "gene_count": self._extract_value(section_text, r"Genes: (\d+)"),
                    "reference_genome": self._extract_value(section_text, r"Reference Genome: (\S+)"),
                    "genome_source": self._extract_value(section_text, r"Reference Genome: \S+ \(source: (\S+)\)"),
                    "gencode_version": self._extract_value(section_text, r"GENCODE Version: (\d+)"),
                    "gencode_source": self._extract_value(section_text, r"GENCODE Version: \d+ \(source: (\S+)\)"),
                    "rna_seq_type": self._extract_value(section_text, r"RNA-seq Type: (.+?) \(confidence:"),
                    "protocol_confidence": self._extract_value(section_text, r"confidence: (\w+)\)"),
                    "data_completeness": self._extract_completeness(section_text)
                }
                
                datasets_info[dataset_name.lower()] = dataset_info
        
        return datasets_info
    
    def _extract_value(self, text, pattern):
        """Extract value using regex pattern"""
        match = re.search(pattern, text)
        return match.group(1) if match else ""
    
    def _extract_completeness(self, text):
        """Extract completeness metrics from log section"""
        completeness = {}
        
        # Find the Quality control checks section
        qc_section_match = re.search(r"Quality control checks:(.*?)(?=\n\d{4}-\d{2}-\d{2}|\Z)", 
                                     text, re.DOTALL)
        
        if qc_section_match:
            qc_section = qc_section_match.group(1)
            
            # Extract percentages for each field
            fields = ["subject_id", "sex", "age", "tissue", "species", 
                      "tissue_ontology", "assay_ontology", "developmental_stage_ontology"]
            
            for field in fields:
                pattern = rf"{field}: (\d+)/(\d+) \((\d+\.\d+)%\)"
                match = re.search(pattern, qc_section)
                
                if match:
                    present, total, percentage = match.groups()
                    completeness[field] = {
                        "present": int(present),
                        "total": int(total),
                        "percentage": float(percentage)
                    }
        
        return completeness
    
    def load_json_metadata(self, json_file):
        """Load metadata from JSON file"""
        try:
            with open(json_file, 'r') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.error(f"Error loading JSON file {json_file}: {e}")
            return {}
    
    def load_anndata_metadata(self, h5ad_file):
        """Load metadata from AnnData file"""
        try:
            adata = sc.read_h5ad(h5ad_file)
            # Extract relevant metadata
            metadata = {
                "uns": dict(adata.uns),
                "obs_keys": list(adata.obs.columns),
                "var_keys": list(adata.var.columns),
                "shape": adata.shape
            }
            
            # Extract sample count and gene count
            metadata["sample_count"] = adata.n_obs
            metadata["gene_count"] = adata.n_vars
            
            # Extract subject and tissue information if available
            if "subject_id" in adata.obs:
                metadata["donor_count"] = len(adata.obs["subject_id"].unique())
            
            if "tissue" in adata.obs:
                metadata["tissue_count"] = len(adata.obs["tissue"].unique())
                metadata["tissue_list"] = adata.obs["tissue"].unique().tolist()
            
            # Extract tissue ontology mapping stats if available
            if "tissue_ontology" in adata.obs and "tissue" in adata.obs:
                tissue_mapping = {}
                for tissue in adata.obs["tissue"].unique():
                    mask = adata.obs["tissue"] == tissue
                    ontology_ids = adata.obs.loc[mask, "tissue_ontology"].unique()
                    if len(ontology_ids) > 0 and not pd.isna(ontology_ids[0]):
                        tissue_mapping[tissue] = ontology_ids[0]
                
                metadata["tissue_mapping"] = tissue_mapping
                metadata["tissue_mapping_count"] = len(tissue_mapping)
                metadata["tissue_mapping_percentage"] = (len(tissue_mapping) / len(adata.obs["tissue"].unique())) * 100
            
            # Extract RNA-seq protocol info if available
            if "rna_seq_protocol" in adata.uns:
                metadata["rna_seq_protocol"] = adata.uns["rna_seq_protocol"]
            
            return metadata
        except Exception as e:
            logger.error(f"Error loading AnnData file {h5ad_file}: {e}")
            return {}
    
    def _extract_tissue_list(self, adata_metadata):
        """Extract unique tissue list from AnnData metadata with ontology mappings"""
        tissues = []
        
        # Get tissue list and mappings if available
        tissue_list = adata_metadata.get("tissue_list", [])
        tissue_mapping = adata_metadata.get("tissue_mapping", {})
        
        for tissue in tissue_list:
            tissue_entry = {
                "label": tissue,
                "ontology": tissue_mapping.get(tissue, "")
            }
            
            # Check if we have a mapping in our default table
            if not tissue_entry["ontology"] and tissue in self.tissue_to_ontology:
                tissue_entry["ontology"] = self.tissue_to_ontology[tissue]
            
            # Count samples for this tissue if possible
            if "obs" in adata_metadata and "tissue" in adata_metadata["obs"]:
                tissue_samples = sum(adata_metadata["obs"]["tissue"] == tissue)
                tissue_entry["sample_count"] = tissue_samples
            
            tissues.append(tissue_entry)
        
        return tissues
    
    def standardize_metadata(self, dataset_name, json_metadata=None, anndata_metadata=None):
        """
        Standardize metadata for a dataset
        
        Args:
            dataset_name: Name of the dataset
            json_metadata: Optional JSON metadata
            anndata_metadata: Optional AnnData metadata
        
        Returns:
            Standardized metadata dictionary
        """
        # Start with template copy
        metadata = self.metadata_template.copy()
        
        # Set dataset info
        metadata["dataset_info"]["dataset_id"] = dataset_name.lower()
        metadata["dataset_info"]["dataset_name"] = f"{dataset_name.upper()} RNA-seq Expression Data"
        
        # Add log-derived data if available
        if dataset_name.lower() in self.log_data:
            log_info = self.log_data[dataset_name.lower()]
            
            # Dataset statistics
            metadata["dataset_statistics"]["sample_count"] = int(log_info.get("sample_count", 0))
            metadata["dataset_statistics"]["gene_count"] = int(log_info.get("gene_count", 0))
            
            # Genomic metadata
            metadata["genomic_metadata"]["reference_genome"] = log_info.get("reference_genome", "")
            metadata["genomic_metadata"]["genome_source"] = log_info.get("genome_source", "")
            metadata["genomic_metadata"]["original_gencode_version"] = log_info.get("gencode_version", "")
            metadata["genomic_metadata"]["harmonized_gencode_version"] = "24"  # Our target version
            
            # Generate mapping notes if versions differ
            if log_info.get("gencode_version", "") != "24":
                metadata["genomic_metadata"]["gencode_mapping_notes"] = (
                    f"Original {dataset_name.upper()} data was aligned using GENCODE v{log_info.get('gencode_version', '')} "
                    f"annotation. For harmonization across datasets, gene IDs were mapped to GENCODE v24 to ensure "
                    f"consistent comparison with other RNA-seq datasets in the collection."
                )
            
            # RNA-seq metadata
            metadata["rna_seq_metadata"]["protocol_type"] = log_info.get("rna_seq_type", "")
            metadata["rna_seq_metadata"]["protocol_confidence"] = log_info.get("protocol_confidence", "")
            
            # Enhanced protocol information for GTEx
            if dataset_name.lower() == "gtex":
                metadata["rna_seq_metadata"]["protocol_type"] = "RNA-seq polyA+"
                metadata["rna_seq_metadata"]["protocol_confidence"] = "high"
                metadata["rna_seq_metadata"]["selection_method"] = "polyA+ selection using Illumina TruSeq protocol"
                metadata["rna_seq_metadata"]["sequencing_platform"] = "Illumina HiSeq 2000/2500"
                metadata["rna_seq_metadata"]["read_configuration"] = "76-bp paired-end reads"
                metadata["rna_seq_metadata"]["strand_specificity"] = "Non-strand specific"
            
            # Data completeness
            if "data_completeness" in log_info:
                for field, values in log_info["data_completeness"].items():
                    metadata["data_completeness"][field] = values["percentage"]
        
        # Add JSON-derived data if available
        if json_metadata:
            # License information
            if "license_terms" in json_metadata:
                metadata["dataset_info"]["license"] = json_metadata["license_terms"]
            
            # DOIs
            if "doi" in json_metadata:
                metadata["dataset_info"]["doi"] = json_metadata["doi"]
            
            # Dataset steward
            if "steward" in json_metadata:
                metadata["dataset_info"]["steward"] = json_metadata["steward"]
            
            # File locations
            if "internal_sources" in json_metadata:
                metadata["file_locations"]["internal"] = json_metadata["internal_sources"]
            
            if "external_sources" in json_metadata:
                metadata["file_locations"]["external"] = json_metadata["external_sources"]
            
            # Related models and benchmarks
            if "models" in json_metadata:
                metadata["related_models_benchmarks"]["models"] = json_metadata["models"]
            
            if "benchmarks" in json_metadata:
                metadata["related_models_benchmarks"]["benchmarks"] = json_metadata["benchmarks"]
            
            # Process extras
            if "extras" in json_metadata:
                for extra in json_metadata["extras"]:
                    key = extra["key"]
                    value = extra["value"]
                    
                    if key == "dataset_name":
                        metadata["dataset_info"]["dataset_name"] = value
                    elif key == "project_description":
                        metadata["dataset_info"]["description"] = value
                    elif key == "published_at":
                        metadata["dataset_info"]["publication_date"] = value
                    elif key == "updated_at":
                        metadata["dataset_info"]["update_date"] = value
                    elif key == "contributors":
                        metadata["dataset_info"]["contributors"] = [value] if isinstance(value, str) else value
                    elif key == "sample_count":
                        metadata["dataset_statistics"]["sample_count"] = int(value)
                    elif key == "donor_count":
                        metadata["dataset_statistics"]["donor_count"] = int(value)
                    elif key == "gene_count":
                        metadata["dataset_statistics"]["gene_count"] = int(value)
                    elif key == "tissue_count":
                        metadata["dataset_statistics"]["tissue_count"] = int(value)
                    elif key == "age_distribution":
                        metadata["demographics"]["age_distribution"] = value
                    elif key == "sex_distribution":
                        metadata["demographics"]["sex_distribution"] = value
                    elif key == "ethnicity_from_publication":
                        metadata["demographics"]["ethnicity_distribution"] = value
                    elif key == "assay":
                        metadata["ontology_mappings"]["assay"] = value
                    elif key == "organism":
                        metadata["ontology_mappings"]["organism"] = value
                    elif key == "development_stage":
                        metadata["ontology_mappings"]["development_stage"] = value
                    elif key == "tissues":
                        metadata["ontology_mappings"]["tissues"] = value
                    elif key == "age":
                        metadata["sample_attributes"]["age_field"] = value
                    elif key == "sex":
                        metadata["sample_attributes"]["sex_field"] = value
                    elif key == "tissue":
                        metadata["sample_attributes"]["tissue_field"] = value
                    elif key == "disease":
                        metadata["sample_attributes"]["disease_field"] = value
                    elif key == "universal_ids":
                        metadata["sample_attributes"]["subject_id_field"] = value
                    else:
                        # Store other fields in dataset_specific
                        metadata["dataset_specific"][key] = value
        
        # Add AnnData-derived data if available
        if anndata_metadata:
            # Override with actual counts from AnnData if available
            if "sample_count" in anndata_metadata:
                metadata["dataset_statistics"]["sample_count"] = anndata_metadata["sample_count"]
            
            if "gene_count" in anndata_metadata:
                metadata["dataset_statistics"]["gene_count"] = anndata_metadata["gene_count"]
            
            if "donor_count" in anndata_metadata:
                metadata["dataset_statistics"]["donor_count"] = anndata_metadata["donor_count"]
            
            if "tissue_count" in anndata_metadata:
                metadata["dataset_statistics"]["tissue_count"] = anndata_metadata["tissue_count"]
            
            # Add tissue ontology mappings
            if "tissue_list" in anndata_metadata:
                tissue_entries = self._extract_tissue_list(anndata_metadata)
                metadata["ontology_mappings"]["tissues"] = tissue_entries
            
            # Add RNA-seq protocol info if available
            if "rna_seq_protocol" in anndata_metadata:
                protocol_info = anndata_metadata["rna_seq_protocol"]
                metadata["rna_seq_metadata"]["protocol_type"] = protocol_info.get("protocol_type", "")
                metadata["rna_seq_metadata"]["protocol_confidence"] = protocol_info.get("protocol_confidence", "")
                metadata["rna_seq_metadata"]["selection_method"] = protocol_info.get("selection_method", "")
                metadata["rna_seq_metadata"]["sequencing_platform"] = protocol_info.get("sequencing_platform", "")
                metadata["rna_seq_metadata"]["read_configuration"] = protocol_info.get("read_configuration", "")
                metadata["rna_seq_metadata"]["strand_specificity"] = protocol_info.get("strand_specificity", "")
            
            # Store anndata-specific keys
            metadata["dataset_specific"]["anndata_obs_keys"] = anndata_metadata.get("obs_keys", [])
            metadata["dataset_specific"]["anndata_var_keys"] = anndata_metadata.get("var_keys", [])
        
        return metadata
    
    def process_dataset(self, dataset_name, json_path=None):
        """
        Process a dataset to collect and standardize metadata
        
        Args:
            dataset_name: Name of the dataset
            json_path: Optional path to JSON metadata file
        """
        logger.info(f"Processing {dataset_name} dataset")
        
        # Load JSON metadata if available
        json_metadata = None
        if json_path and os.path.exists(json_path):
            json_metadata = self.load_json_metadata(json_path)
        
        # Find and load AnnData file
        h5ad_path = self.data_dir / f"{dataset_name.lower()}_standardized_v2.h5ad"
        if not h5ad_path.exists():
            h5ad_path = self.data_dir / f"{dataset_name.lower()}_standardized.h5ad"
        
        anndata_metadata = None
        if h5ad_path.exists():
            anndata_metadata = self.load_anndata_metadata(h5ad_path)
        
        # Standardize metadata
        metadata = self.standardize_metadata(dataset_name, json_metadata, anndata_metadata)
        
        # Save standardized metadata
        output_path = self.output_dir / f"{dataset_name.lower()}_metadata.json"
        with open(output_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Saved standardized metadata to {output_path}")
        
        return metadata
    
    def process_all_datasets(self, json_dir=None):
        """
        Process all datasets
        
        Args:
            json_dir: Optional directory containing JSON metadata files
        """
        all_metadata = {}
        
        for dataset_name in DATASETS:
            json_path = None
            if json_dir:
                # Try to find matching JSON file
                potential_paths = [
                    os.path.join(json_dir, f"{dataset_name.lower()}_metadata.json"),
                    os.path.join(json_dir, f"{dataset_name}_metadata.json"),
                    os.path.join(json_dir, f"{dataset_name.lower()}.json"),
                    os.path.join(json_dir, f"{dataset_name}.json")
                ]
                
                for path in potential_paths:
                    if os.path.exists(path):
                        json_path = path
                        break
            
            metadata = self.process_dataset(dataset_name, json_path)
            all_metadata[dataset_name] = metadata
        
        # Create a summary file
        summary = {
            "collection_date": datetime.now().strftime("%Y-%m-%d"),
            "datasets": list(all_metadata.keys()),
            "summary": {
                dataset: {
                    "sample_count": meta["dataset_statistics"]["sample_count"],
                    "gene_count": meta["dataset_statistics"]["gene_count"],
                    "reference_genome": meta["genomic_metadata"]["reference_genome"],
                    "original_gencode_version": meta["genomic_metadata"]["original_gencode_version"],
                    "harmonized_gencode_version": meta["genomic_metadata"]["harmonized_gencode_version"],
                    "tissue_ontology_completion": meta["data_completeness"]["tissue_ontology"]
                }
                for dataset, meta in all_metadata.items()
            }
        }
        
        # Calculate overall tissue ontology completion
        total_samples = sum(meta["dataset_statistics"]["sample_count"] for meta in all_metadata.values())
        mapped_samples = sum(meta["dataset_statistics"]["sample_count"] * meta["data_completeness"]["tissue_ontology"] / 100 
                             for meta in all_metadata.values())
        completion_percentage = (mapped_samples / total_samples * 100) if total_samples > 0 else 0
        
        summary["overall_tissue_ontology_completion"] = round(completion_percentage, 1)
        
        # Add action items for improving mapping
        summary["action_items"] = {
            "tissue_ontology_mapping": [
                "Improve GTEx tissue-to-UBERON mapping",
                "Create missing mappings for brain subregions",
                "Handle cell lines using Cell Ontology (CL)",
                "Validate mappings with domain experts"
            ],
            "protocol_information": [
                "Extract more detailed protocol information for non-GTEx datasets",
                "Add sequencing platform information where missing",
                "Clarify strand specificity for all datasets"
            ],
            "gencode_harmonization": [
                "Document gene ID mapping process between GENCODE versions",
                "Validate gene ID consistency across datasets",
                "Create gene ID concordance table"
            ]
        }
        
        summary_path = self.output_dir / "metadata_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Saved metadata summary to {summary_path}")
        
        # Create tissue mapping summary
        self._generate_tissue_mapping_report(all_metadata)
    
    def _generate_tissue_mapping_report(self, all_metadata):
        """Generate a report on tissue-to-ontology mapping status"""
        if not self.mapping_dir:
            return
            
        # Collect all unique tissues across datasets
        all_tissues = set()
        tissue_datasets = {}
        tissue_mappings = {}
        
        for dataset, metadata in all_metadata.items():
            if "ontology_mappings" in metadata and "tissues" in metadata["ontology_mappings"]:
                for tissue_entry in metadata["ontology_mappings"]["tissues"]:
                    tissue_name = tissue_entry.get("label", "")
                    if not tissue_name:
                        continue
                        
                    all_tissues.add(tissue_name)
                    
                    # Track which datasets contain this tissue
                    if tissue_name not in tissue_datasets:
                        tissue_datasets[tissue_name] = []
                    tissue_datasets[tissue_name].append(dataset)
                    
                    # Track ontology mapping if available
                    if "ontology" in tissue_entry and tissue_entry["ontology"]:
                        tissue_mappings[tissue_name] = tissue_entry["ontology"]
        
        # Create mapping status report
        mapping_rows = []
        for tissue in sorted(all_tissues):
            mapping_status = "mapped" if tissue in tissue_mappings else "unmapped"
            datasets = ", ".join(tissue_datasets.get(tissue, []))
            ontology_id = tissue_mappings.get(tissue, "")
            
            mapping_rows.append({
                "tissue_name": tissue,
                "ontology_id": ontology_id,
                "mapping_status": mapping_status,
                "datasets": datasets
            })
        
        # Save mapping status report
        report_path = self.mapping_dir / "tissue_mapping_status.csv"
        df = pd.DataFrame(mapping_rows)
        df.to_csv(report_path, index=False)
        
        # Calculate mapping statistics
        total_tissues = len(all_tissues)
        mapped_tissues = len(tissue_mappings)
        mapping_percentage = (mapped_tissues / total_tissues * 100) if total_tissues > 0 else 0
        
        # Create summary statistics
        summary = {
            "total_tissues": total_tissues,
            "mapped_tissues": mapped_tissues,
            "mapping_percentage": round(mapping_percentage, 1),
            "unmapped_tissues": sorted(list(all_tissues - set(tissue_mappings.keys())))
        }
        
        summary_path = self.mapping_dir / "tissue_mapping_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Generated tissue mapping report: {report_path}")
        logger.info(f"Tissue mapping status: {mapped_tissues}/{total_tissues} tissues mapped ({mapping_percentage:.1f}%)")

def main():
    parser = argparse.ArgumentParser(description="Collect and standardize RNA-seq metadata")
    parser.add_argument("--data-dir", required=True, help="Directory containing h5ad files")
    parser.add_argument("--log-file", required=True, help="Path to standardization log file")
    parser.add_argument("--output-dir", required=True, help="Directory to save standardized metadata")
    parser.add_argument("--json-dir", help="Directory containing JSON metadata files")
    parser.add_argument("--mapping-dir", help="Directory for ontology mapping files")
    
    args = parser.parse_args()
    
    collector = MetadataCollector(
        args.data_dir, 
        args.log_file, 
        args.output_dir,
        args.mapping_dir
    )
    collector.process_all_datasets(args.json_dir)

if __name__ == "__main__":
    main()