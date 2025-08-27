#!/usr/bin/env python3
"""
Dataset Manager for RNA-seq Gene Quantification Data
---------------------------------------------------

This module handles dataset metadata management and file access operations.
It provides functions to:
- Fetch and cache experiment metadata from ENCODE
- List available datasets
- Get TSV files for experiments
- Download TSV files
"""

import json
import requests
from pathlib import Path
import os

class DatasetManager:
    def __init__(self, config_file=None):
        """Initialize the dataset manager with configuration"""
        # Load configuration
        if config_file is None:
            script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
            config_file = script_dir.parent / "config" / "settings.json"
            
        with open(config_file, 'r') as f:
            self.config = json.load(f)
            
        # Set up paths
        self.base_dir = Path(self.config['base_dir'])
        self.raw_data_dir = Path(self.config['raw_data_dir'])
        self.analysis_dir = Path(self.config['analysis_dir'])
        self.metadata_dir = Path(self.config['metadata_dir'])
        
        # Create directories if they don't exist
        self.raw_data_dir.mkdir(parents=True, exist_ok=True)
        self.analysis_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_dir.mkdir(parents=True, exist_ok=True)
        
        # Load cell line information
        self.cell_lines = self.config['cell_lines']
        
    def list_cell_lines(self):
        """List all cell lines/datasets in the configuration"""
        return list(self.cell_lines.keys())
    
    def get_cell_line_info(self, cell_line):
        """Get information about a specific cell line"""
        if cell_line in self.cell_lines:
            return self.cell_lines[cell_line]
        return None
    
    def add_cell_line(self, cell_line, experiment_id, rna_type, is_tissue=False, **kwargs):
        """Add a new cell line or dataset to the configuration"""
        info = {
            'experiment_id': experiment_id,
            'rna_type': rna_type,
            'is_tissue': is_tissue
        }
        info.update(kwargs)
        self.cell_lines[cell_line] = info
        
        # Update the configuration file
        self.config['cell_lines'] = self.cell_lines
        with open(self.base_dir / "config" / "settings.json", 'w') as f:
            json.dump(self.config, f, indent=2)
            
        return True
    
    def fetch_experiment_metadata(self, experiment_id):
        """Fetch experiment metadata from ENCODE API"""
        metadata_file = self.metadata_dir / f"{experiment_id}.json"
        
        # Try to load from cache first
        if metadata_file.exists():
            try:
                with open(metadata_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Error loading cached metadata: {e}")
        
        # Fetch from ENCODE API
        try:
            print(f"Fetching metadata for experiment {experiment_id}...")
            url = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
            response = requests.get(url)
            if response.status_code == 200:
                metadata = response.json()
                # Save to cache
                with open(metadata_file, 'w') as f:
                    json.dump(metadata, f, indent=2)
                return metadata
            else:
                print(f"Error fetching metadata: HTTP {response.status_code}")
        except Exception as e:
            print(f"Error fetching metadata: {e}")
        
        return None
    
    def get_tsv_files(self, experiment_id):
        """Find gene quantification TSV files for a given experiment ID"""
        # Fetch metadata if not already cached
        metadata = self.fetch_experiment_metadata(experiment_id)
        
        if not metadata:
            return []
        
        # Look for gene quantification files
        tsv_files = []
        for file_info in metadata.get('files', []):
            if (file_info.get('output_type') == 'gene quantifications' and 
                file_info.get('file_format') == 'tsv'):
                
                rep_info = file_info.get('biological_replicates', [])
                rep_num = rep_info[0] if rep_info else 0
                
                tsv_files.append({
                    'accession': file_info.get('accession'),
                    'url': f"https://www.encodeproject.org/files/{file_info.get('accession')}/@@download/{file_info.get('accession')}.tsv",
                    'rep': rep_num
                })
        
        return tsv_files
    
    def download_tsv(self, accession, url, cell_line):
        """Download a TSV file"""
        output_dir = self.raw_data_dir / cell_line
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{accession}.tsv"
        
        if output_file.exists():
            print(f"File {accession}.tsv already exists for {cell_line}")
            return output_file
        
        try:
            print(f"Downloading {accession}.tsv for {cell_line}...")
            response = requests.get(url)
            if response.status_code == 200:
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded {accession}.tsv")
                return output_file
            else:
                print(f"Error downloading file: HTTP {response.status_code}")
        except Exception as e:
            print(f"Error downloading {accession}: {str(e)}")
        
        return None
    
    def add_entex_datasets(self):
        """Add ENTEx tissue datasets to the configuration"""
        entex_metadata_file = self.metadata_dir / "entex_experiments.json"
        
        # Try to load from cache first
        if entex_metadata_file.exists():
            with open(entex_metadata_file, 'r') as f:
                entex_data = json.load(f)
        else:
            # Query ENCODE API for ENTEx datasets
            entex_url = "https://www.encodeproject.org/search/?type=Experiment&internal_tags=ENTEx&assay_title=total+RNA-seq&status=released&format=json"
            response = requests.get(entex_url)
            if response.status_code == 200:
                entex_data = response.json()
                # Save to cache
                with open(entex_metadata_file, 'w') as f:
                    json.dump(entex_data, f, indent=2)
            else:
                print(f"Error fetching ENTEx metadata: HTTP {response.status_code}")
                return 0
        
        # Process and add datasets
        added = 0
        for exp in entex_data.get('@graph', []):
            exp_id = exp.get('@id').split('/')[2]
            biosample = exp.get('biosample_ontology', {}).get('term_name', 'unknown')
            donor = exp.get('biosample', {}).get('donor', {}).get('accession', 'unknown')
            
            dataset_name = f"ENTEx_{biosample}_{donor}"
            if dataset_name not in self.cell_lines:
                self.add_cell_line(
                    dataset_name,
                    experiment_id=exp_id,
                    rna_type='total',
                    is_tissue=True,
                    biosample=biosample,
                    donor=donor
                )
                added += 1
        
        print(f"Added {added} ENTEx tissue datasets")
        return added