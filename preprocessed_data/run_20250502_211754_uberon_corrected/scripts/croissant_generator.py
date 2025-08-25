#!/usr/bin/env python3
"""
Croissant JSON-LD Generator for h5ad Genomics Datasets

This script generates Croissant metadata files for genomics datasets stored in h5ad format.
It extracts metadata from AnnData objects and creates standardized JSON-LD files following
the MLCommons Croissant specification.

Based on:
- Croissant Format Specification: https://docs.mlcommons.org/croissant/docs/croissant-spec.html
- Cell by Gene Schema: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md
"""

import scanpy as sc
import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional
import uuid


class DatasetInfo:
    """Container for dataset-specific information from literature."""
    
    GTEx = {
        'name': 'GTEx (Genotype-Tissue Expression)',
        'description': 'The Genotype-Tissue Expression (GTEx) project characterizes human transcriptomes within and across individuals for a wide variety of primary tissues and cell types. This dataset contains RNA-seq data from 838 individuals across 49 tissues.',
        'url': 'https://www.gtexportal.org/',
        'citations': [
            'Aguet, F. et al. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318-1330 (2020).',
            'The GTEx Consortium. Genetic effects on gene expression across human tissues. Nature 550, 204-213 (2017).'
        ],
        'license': 'dbGaP Authorized Access',
        'creator': 'The GTEx Consortium',
        'keywords': ['gene expression', 'eQTL', 'tissue', 'RNA-seq', 'human', 'multi-tissue'],
        'assay_type': 'RNA sequencing',
        'organism': 'Homo sapiens',
        'tissue_types': 'Multi-tissue (49 tissues)',
        'subject_demographics': 'Adults aged 20-79, predominantly male (67%)'
    }
    
    ENCODE = {
        'name': 'ENCODE (Encyclopedia of DNA Elements)',
        'description': 'The Encyclopedia of DNA Elements (ENCODE) project systematically maps regions of transcription, transcription factor association, chromatin structure and histone modification. This phase III dataset includes RNA-seq data from various cell lines and tissues.',
        'url': 'https://www.encodeproject.org/',
        'citations': [
            'ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 583, 699-710 (2020).',
            'ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57-74 (2012).'
        ],
        'license': 'CC0',
        'creator': 'ENCODE Project Consortium',
        'keywords': ['gene expression', 'functional genomics', 'transcription', 'cell lines', 'RNA-seq'],
        'assay_type': 'RNA sequencing',
        'organism': 'Homo sapiens',
        'tissue_types': 'Various cell lines and primary tissues',
        'subject_demographics': 'Various cell lines and donor tissues'
    }
    
    ADNI = {
        'name': 'ADNI (Alzheimer\'s Disease Neuroimaging Initiative)',
        'description': 'The Alzheimer\'s Disease Neuroimaging Initiative is a longitudinal multicenter study designed to develop clinical, imaging, genetic, and biochemical biomarkers for the early detection and tracking of Alzheimer\'s disease. This dataset contains microarray gene expression data from brain tissue.',
        'url': 'https://adni.loni.usc.edu/',
        'citations': [
            'Petersen RC, et al. Alzheimer\'s Disease Neuroimaging Initiative (ADNI): Clinical characterization. Neurology 74, 201-209 (2010).'
        ],
        'license': 'ADNI Data Use Agreement',
        'creator': 'Alzheimer\'s Disease Neuroimaging Initiative',
        'keywords': ['Alzheimer disease', 'neuroimaging', 'microarray', 'brain', 'neurodegenerative'],
        'assay_type': 'Microarray',
        'organism': 'Homo sapiens',
        'tissue_types': 'Brain tissue',
        'subject_demographics': 'Adults aged 55-90, cognitively normal, MCI, and AD patients'
    }
    
    MAGE = {
        'name': 'MAGE (Multi-ancestry Analysis of Gene Expression)',
        'description': 'RNA-seq data from lymphoblastoid cell lines derived from 731 individuals from the 1000 Genomes Project, representing 26 globally-distributed populations across five continental groups. This dataset provides a large, geographically diverse resource for studying human transcriptome variation.',
        'url': 'https://github.com/mccoy-lab/MAGE',
        'citations': [
            'Taylor DJ, et al. Sources of gene expression variation in a globally diverse human cohort. Nature (2024). DOI: 10.1038/s41586-024-07708-2'
        ],
        'license': 'CC0',
        'creator': 'McCoy Lab, Johns Hopkins University',
        'keywords': ['gene expression', 'population genetics', 'diversity', 'lymphoblastoid cell lines', 'RNA-seq'],
        'assay_type': 'RNA sequencing',
        'organism': 'Homo sapiens',
        'tissue_types': 'Lymphoblastoid cell lines',
        'subject_demographics': '731 individuals from 26 populations across 5 continental groups'
    }


class CroissantGenerator:
    """Generator for Croissant JSON-LD metadata files from h5ad datasets."""
    
    def __init__(self, base_url: str = "https://example.com/datasets"):
        self.base_url = base_url
        self.dataset_info = DatasetInfo()
        
    def extract_h5ad_metadata(self, file_path: Path) -> Dict[str, Any]:
        """Extract metadata from an h5ad file."""
        
        adata = sc.read_h5ad(file_path)
        
        # Basic structure information
        metadata = {
            'shape': adata.shape,
            'n_genes': adata.n_vars,
            'n_samples': adata.n_obs,
            'obs_columns': list(adata.obs.columns),
            'var_columns': list(adata.var.columns),
            'uns_keys': list(adata.uns.keys()) if adata.uns else [],
            'layers': list(adata.layers.keys()) if adata.layers else [],
        }
        
        # Extract key metadata from observations
        obs_metadata = {}
        important_cols = ['dataset', 'tissue', 'cell_type', 'sex', 'age', 'ethnicity', 
                         'species', 'assay', 'data_type', 'expression_unit']
        
        for col in important_cols:
            if col in adata.obs.columns:
                unique_values = adata.obs[col].unique()
                # Convert numpy types to Python native types for JSON serialization
                values_list = []
                for val in unique_values[:10]:  # Limit to first 10 values
                    if pd.isna(val):
                        values_list.append(None)
                    elif hasattr(val, 'item'):  # numpy scalar
                        values_list.append(val.item())
                    else:
                        values_list.append(str(val))
                
                obs_metadata[col] = {
                    'unique_count': int(len(unique_values)),
                    'values': values_list
                }
        
        # Extract gene information
        var_metadata = {}
        if 'gene_type' in adata.var.columns:
            gene_types = adata.var['gene_type'].value_counts()
            # Convert to regular Python dict with int values
            var_metadata['gene_types'] = {str(k): int(v) for k, v in gene_types.head(10).items()}
        
        if 'chromosome' in adata.var.columns:
            chromosomes = adata.var['chromosome'].unique()
            # Convert to regular Python list with string values
            var_metadata['chromosomes'] = [str(c) for c in chromosomes if pd.notna(c)]
        
        # Extract processing information from .uns
        processing_info = {}
        if adata.uns:
            for key in ['processing_date', 'gencode_version', 'reference_genome']:
                if key in adata.uns:
                    value = adata.uns[key]
                    # Convert numpy types to Python native types
                    if hasattr(value, 'item'):
                        processing_info[key] = value.item()
                    else:
                        processing_info[key] = str(value)
        
        metadata.update({
            'observations': obs_metadata,
            'variables': var_metadata,
            'processing': processing_info
        })
        
        return metadata
    
    def get_dataset_info(self, dataset_name: str) -> Dict[str, Any]:
        """Get dataset-specific information from literature."""
        dataset_name = dataset_name.upper()
        
        # Map common variations to standard names
        name_mapping = {
            'GTEX': 'GTEx',
            'ENCODE': 'ENCODE', 
            'ADNI': 'ADNI',
            'MAGE': 'MAGE'
        }
        
        mapped_name = name_mapping.get(dataset_name, dataset_name)
        if hasattr(self.dataset_info, mapped_name):
            return getattr(self.dataset_info, mapped_name)
        else:
            # Default information for unknown datasets
            return {
                'name': dataset_name,
                'description': f'Gene expression dataset: {dataset_name}',
                'url': self.base_url,
                'citations': [],
                'license': 'Unknown',
                'creator': 'Unknown',
                'keywords': ['gene expression', 'genomics'],
                'assay_type': 'Unknown',
                'organism': 'Homo sapiens',
                'tissue_types': 'Various',
                'subject_demographics': 'Unknown'
            }
    
    def generate_croissant_metadata(self, file_path: Path) -> Dict[str, Any]:
        """Generate Croissant JSON-LD metadata for an h5ad file."""
        
        # Extract dataset name from filename
        dataset_name = file_path.stem.replace('_standardized_preprocessed', '').upper()
        
        # Get h5ad metadata
        h5ad_meta = self.extract_h5ad_metadata(file_path)
        
        # Get dataset-specific information
        dataset_info = self.get_dataset_info(dataset_name)
        
        # Generate unique identifier
        dataset_id = f"{dataset_name.lower()}_h5ad_{datetime.now().strftime('%Y%m%d')}"
        
        # Build Croissant metadata
        croissant_metadata = {
            "@context": {
                "@language": "en",
                "@vocab": "https://schema.org/",
                "citeAs": "cr:citeAs",
                "column": "cr:column",
                "conformsTo": "dct:conformsTo",
                "cr": "http://mlcommons.org/croissant/",
                "rai": "http://mlcommons.org/croissant/RAI/",
                "data": {
                    "@id": "cr:data",
                    "@type": "@json"
                },
                "dataType": {
                    "@id": "cr:dataType",
                    "@type": "@vocab"
                },
                "dct": "http://purl.org/dc/terms/",
                "examples": {
                    "@id": "cr:examples",
                    "@type": "@json"
                },
                "extract": "cr:extract",
                "field": "cr:field",
                "fileProperty": "cr:fileProperty",
                "fileObject": "cr:fileObject",
                "fileSet": "cr:fileSet",
                "format": "cr:format",
                "includes": "cr:includes",
                "isLiveDataset": "cr:isLiveDataset",
                "jsonPath": "cr:jsonPath",
                "key": "cr:key",
                "md5": "cr:md5",
                "parentField": "cr:parentField",
                "path": "cr:path",
                "recordSet": "cr:recordSet",
                "references": "cr:references",
                "regex": "cr:regex",
                "repeated": "cr:repeated",
                "replace": "cr:replace",
                "sc": "https://schema.org/",
                "separator": "cr:separator",
                "source": "cr:source",
                "subField": "cr:subField",
                "transform": "cr:transform",
                "wd": "https://www.wikidata.org/wiki/"
            },
            "@type": "sc:Dataset",
            "@id": f"{self.base_url}/{dataset_id}",
            "dct:conformsTo": "http://mlcommons.org/croissant/1.0",
            "name": f"{dataset_info['name']} - Standardized H5AD",
            "description": f"{dataset_info['description']} This dataset has been standardized and preprocessed for analysis, containing {h5ad_meta['n_samples']:,} samples and {h5ad_meta['n_genes']:,} genes.",
            "url": dataset_info['url'],
            "license": dataset_info['license'],
            "creator": {
                "@type": "Organization",
                "name": dataset_info['creator']
            },
            "citation": dataset_info['citations'],
            "keywords": dataset_info['keywords'],
            "datePublished": datetime.now().isoformat(),
            "version": "1.0",
            "identifier": dataset_id,
            
            # Genomics-specific metadata
            "sc:organism": dataset_info['organism'],
            "sc:assayType": dataset_info['assay_type'],
            "sc:tissueType": dataset_info['tissue_types'],
            
            # File distribution
            "distribution": [
                {
                    "@type": "cr:FileObject",
                    "name": file_path.name,
                    "description": f"H5AD file containing gene expression data for {dataset_name}",
                    "contentUrl": f"{self.base_url}/{file_path.name}",
                    "encodingFormat": "application/x-hdf5",
                    "contentSize": file_path.stat().st_size if file_path.exists() else None,
                    "sha256": self._calculate_file_hash(file_path) if file_path.exists() else None
                }
            ],
            
            # Data structure description
            "recordSet": [
                {
                    "@type": "cr:RecordSet",
                    "name": "gene_expression_matrix",
                    "description": "Gene expression matrix with samples as observations and genes as variables",
                    "field": [
                        {
                            "@type": "cr:Field",
                            "name": "sample_id",
                            "description": "Unique identifier for each sample",
                            "dataType": "sc:Text"
                        },
                        {
                            "@type": "cr:Field", 
                            "name": "gene_expression",
                            "description": f"Gene expression values ({h5ad_meta.get('observations', {}).get('expression_unit', {}).get('values', ['Unknown'])[0]})",
                            "dataType": "sc:Float"
                        }
                    ]
                },
                {
                    "@type": "cr:RecordSet",
                    "name": "sample_metadata",
                    "description": "Sample-level metadata including tissue type, demographics, and experimental conditions",
                    "field": self._generate_metadata_fields(h5ad_meta['obs_columns'])
                },
                {
                    "@type": "cr:RecordSet",
                    "name": "gene_metadata",
                    "description": "Gene-level metadata including gene names, types, and chromosomal locations",
                    "field": self._generate_metadata_fields(h5ad_meta['var_columns'])
                }
            ],
            
            # Provenance and processing information
            "sc:processingInformation": {
                "standardization": "Cell by Gene schema v3.0.0",
                "preprocessing": "Standardized and corrected for tissue ontology mappings",
                "processing_date": h5ad_meta.get('processing', {}).get('processing_date', 'Unknown'),
                "reference_genome": h5ad_meta.get('processing', {}).get('reference_genome', 'Unknown'),
                "gencode_version": h5ad_meta.get('processing', {}).get('gencode_version', 'Unknown')
            },
            
            # Dataset statistics
            "sc:datasetStatistics": {
                "numberOfSamples": h5ad_meta['n_samples'],
                "numberOfGenes": h5ad_meta['n_genes'],
                "dataShape": h5ad_meta['shape'],
                "tissueTypes": list(h5ad_meta.get('observations', {}).get('tissue', {}).get('values', [])),
                "subjectDemographics": dataset_info['subject_demographics']
            }
        }
        
        return croissant_metadata
    
    def _generate_metadata_fields(self, column_names: List[str]) -> List[Dict[str, Any]]:
        """Generate field descriptions for metadata columns."""
        
        field_descriptions = {
            'sample_id': 'Unique sample identifier',
            'subject_id': 'Unique subject identifier',
            'tissue': 'Tissue or cell type',
            'tissue_ontology': 'Tissue ontology term ID',
            'cell_type': 'Cell type annotation',
            'sex': 'Biological sex',
            'age': 'Age or age range',
            'ethnicity': 'Ethnicity or ancestry',
            'disease': 'Disease status or condition',
            'dataset': 'Source dataset name',
            'assay': 'Experimental assay type',
            'data_type': 'Data type (e.g., RNA-seq, microarray)',
            'expression_unit': 'Units of expression measurement',
            'gene_id': 'Unique gene identifier',
            'gene_name': 'Gene symbol or name',
            'gene_type': 'Gene biotype classification',
            'chromosome': 'Chromosomal location'
        }
        
        fields = []
        for col in column_names:
            description = field_descriptions.get(col, f"Metadata field: {col}")
            
            # Infer data type from column name
            if col in ['age', 'rna_integrity_number', 'ischemic_time']:
                data_type = "sc:Float"
            elif col.endswith('_id') or col in ['sample_id', 'subject_id', 'gene_id']:
                data_type = "sc:Text" 
            else:
                data_type = "sc:Text"  # Default to text
            
            fields.append({
                "@type": "cr:Field",
                "name": col,
                "description": description,
                "dataType": data_type
            })
        
        return fields
    
    def _calculate_file_hash(self, file_path: Path) -> Optional[str]:
        """Calculate SHA256 hash of a file (placeholder implementation)."""
        # In a production environment, implement actual hash calculation
        return None
    
    def save_croissant_metadata(self, metadata: Dict[str, Any], output_path: Path) -> None:
        """Save Croissant metadata to a JSON-LD file."""
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        
        print(f"Croissant metadata saved to: {output_path}")


def main():
    """Main function to generate Croissant metadata for all h5ad files."""
    
    # Initialize generator
    generator = CroissantGenerator()
    
    # Find all h5ad files in current directory
    current_dir = Path('.')
    h5ad_files = list(current_dir.glob('*.h5ad'))
    
    if not h5ad_files:
        print("No h5ad files found in current directory")
        return
    
    print(f"Found {len(h5ad_files)} h5ad files")
    
    # Generate Croissant metadata for each file
    for h5ad_file in sorted(h5ad_files):
        print(f"\nProcessing: {h5ad_file.name}")
        
        try:
            # Generate metadata
            metadata = generator.generate_croissant_metadata(h5ad_file)
            
            # Create output filename
            output_name = h5ad_file.stem + '_croissant.jsonld'
            output_path = current_dir / output_name
            
            # Save metadata
            generator.save_croissant_metadata(metadata, output_path)
            
        except Exception as e:
            print(f"Error processing {h5ad_file.name}: {str(e)}")
    
    print(f"\nCroissant metadata generation complete!")


if __name__ == "__main__":
    main()