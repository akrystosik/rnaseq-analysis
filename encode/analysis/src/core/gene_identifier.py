#!/usr/bin/env python3
"""
Gene Identifier for RNA-seq Gene Quantification Data
---------------------------------------------------

This module handles gene identification and mapping, including:
- Mapping between gene symbols and Ensembl IDs
- Extracting gene expression data from TSV files
- Building comprehensive gene mapping databases
"""

import json
import pandas as pd
import os
from pathlib import Path
import re
import requests

class GeneIdentifier:
    def __init__(self, config_file=None):
        """Initialize the gene identifier with configuration"""
        # Load configuration
        if config_file is None:
            script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
            config_file = script_dir.parent / "config" / "settings.json"
            
        with open(config_file, 'r') as f:
            self.config = json.load(f)
            
        # Set up paths
        self.base_dir = Path(self.config['base_dir'])
        self.metadata_dir = Path(self.config['metadata_dir'])
        
        # Ensure metadata directory exists
        self.metadata_dir.mkdir(parents=True, exist_ok=True)
        
        # Load housekeeping genes
        self.housekeeping_genes = self.config['housekeeping_genes']
        
        # Load known gene mappings from config
        self.known_gene_mappings = self.config.get('known_gene_mappings', {})
        
        # Initialize gene ID mapping
        self.gene_mapping = {}
        self.load_gene_mapping()
    
    def load_gene_mapping(self):
        """Load gene mapping from file or initialize from known mappings"""
        mapping_file = self.metadata_dir / "gene_mapping.json"
        
        if mapping_file.exists():
            try:
                with open(mapping_file, 'r') as f:
                    self.gene_mapping = json.load(f)
                print(f"Loaded mapping for {len(self.gene_mapping)} genes")
            except Exception as e:
                print(f"Error loading gene mapping: {e}")
                self.gene_mapping = self.known_gene_mappings.copy()
        else:
            # Initialize with known mappings from config
            self.gene_mapping = self.known_gene_mappings.copy()
    
    def save_gene_mapping(self):
        """Save gene mapping to file"""
        mapping_file = self.metadata_dir / "gene_mapping.json"
        
        with open(mapping_file, 'w') as f:
            json.dump(self.gene_mapping, f, indent=2)
        
        print(f"Saved mapping for {len(self.gene_mapping)} genes")
    
    def add_mapping(self, gene_symbol, ensembl_id):
        """Add a new gene mapping"""
        self.gene_mapping[gene_symbol] = ensembl_id
    
    def get_ensembl_id(self, gene_symbol):
        """Get Ensembl ID for a gene symbol"""
        return self.gene_mapping.get(gene_symbol)
    
    def extract_gene_symbol(self, transcript_id):
        """Extract gene symbol from transcript_id format like 'ENST00000123456|GENE_SYMBOL'"""
        if isinstance(transcript_id, str) and '|' in transcript_id:
            return transcript_id.split('|')[1]
        return None
    
    def get_base_ensembl_id(self, gene_id):
        """Strip version number from an Ensembl gene ID if present"""
        if isinstance(gene_id, str) and gene_id.startswith("ENSG"):
            return gene_id.split('.')[0]
        return gene_id
    
    def parse_tsv(self, file_path):
        """Parse a gene quantification TSV file with improved string handling"""
        try:
            df = pd.read_csv(file_path, sep='\t')
            
            # Immediately convert gene_id to string to avoid issues with string operations
            df['gene_id'] = df['gene_id'].astype(str)
            
            # Check columns
            required_cols = ['gene_id', 'TPM']
            if not all(col in df.columns for col in required_cols):
                # Try to find alternative columns
                if 'gene' in df.columns and 'TPM' in df.columns:
                    df = df.rename(columns={'gene': 'gene_id'})
                elif 'gene_id' in df.columns and 'FPKM' in df.columns and 'TPM' not in df.columns:
                    # Some older ENCODE files have FPKM but not TPM
                    df = df.rename(columns={'FPKM': 'TPM'})
                    print(f"Using FPKM as TPM for {file_path}")
                else:
                    print(f"Required columns not found in {file_path}")
                    print(f"Available columns: {df.columns.tolist()}")
                    return None
            
            # If transcript_id(s) column is present, extract gene symbols from it
            if 'transcript_id(s)' in df.columns:
                df['gene_symbol'] = df['transcript_id(s)'].apply(self.extract_gene_symbol)
            
            return df
            
        except Exception as e:
            print(f"Error parsing {file_path}: {str(e)}")
        
        return None
    
    def extract_gene_data(self, df, gene_symbol):
        """Extract data for a specific gene from the dataframe using improved identification methods"""
        if df is None:
            return 0
        
        # Try using the gene_symbol column if available (from transcript_id extraction)
        if 'gene_symbol' in df.columns:
            symbol_match = df[df['gene_symbol'] == gene_symbol]
            if not symbol_match.empty:
                return symbol_match['TPM'].values[0]
        
        # Try using the Ensembl ID if we have it (with version number handling)
        ensembl_id = self.get_ensembl_id(gene_symbol)
        if ensembl_id:
            # Remove version numbers from gene_id for comparison
            base_gene_ids = df['gene_id'].apply(self.get_base_ensembl_id)
            mask = base_gene_ids == ensembl_id
            matched = df[mask]
            
            if not matched.empty:
                return matched['TPM'].values[0]
        
        # Fallbacks if the ensembl ID approach doesn't work
        
        # Try direct string match with gene symbol
        mask = df['gene_id'].str.contains(gene_symbol, case=False, na=False)
        matched = df[mask]
        
        if not matched.empty:
            return matched['TPM'].values[0]
        
        # Try with transcript ID column
        if 'transcript_id(s)' in df.columns:
            mask = df['transcript_id(s)'].str.contains(gene_symbol, case=False, na=False)
            matched = df[mask]
            
            if not matched.empty:
                return matched['TPM'].values[0]
        
        # If we still couldn't find it, return 0
        return 0
    
    def build_mapping_from_file(self, file_path):
        """Build gene mapping from a single TSV file"""
        df = self.parse_tsv(file_path)
        
        if df is None:
            return 0
        
        added = 0
        if 'gene_id' in df.columns and 'transcript_id(s)' in df.columns:
            for _, row in df.iterrows():
                gene_id = row['gene_id']
                if isinstance(gene_id, str) and gene_id.startswith('ENSG'):
                    # Extract base ID (without version)
                    base_id = self.get_base_ensembl_id(gene_id)
                    
                    # Extract gene symbol from transcript info
                    transcript_info = row['transcript_id(s)']
                    if isinstance(transcript_info, str) and '|' in transcript_info:
                        symbol = transcript_info.split('|')[1]
                        
                        # Only add if not already present
                        if symbol not in self.gene_mapping:
                            self.gene_mapping[symbol] = base_id
                            added += 1
        
        return added
    
    def build_comprehensive_mapping(self, tsv_files):
        """Build a comprehensive mapping from multiple TSV files"""
        total_added = 0
        
        for file_path in tsv_files:
            added = self.build_mapping_from_file(file_path)
            total_added += added
            
            # If we've added enough genes, stop
            if total_added > 20000:
                break
        
        # Save the mapping
        self.save_gene_mapping()
        
        return total_added
    
    def fetch_mapping_from_biomart(self, gene_symbols):
        """Fetch gene mapping from BioMart for a list of gene symbols"""
        # Only try this if we have requests library
        try:
            # BioMart URL
            url = "https://ensembl.org/biomart/martservice"
            
            # Prepare XML query
            xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
                <Dataset name="hsapiens_gene_ensembl" interface="default">
                    <Filter name="external_gene_name" value="{','.join(gene_symbols)}"/>
                    <Attribute name="ensembl_gene_id"/>
                    <Attribute name="external_gene_name"/>
                </Dataset>
            </Query>
            """
            
            # Make request
            response = requests.get(url, params={'query': xml_query})
            
            if response.status_code == 200:
                # Parse response
                lines = response.text.strip().split('\n')
                if len(lines) > 1:  # At least header + 1 data row
                    header = lines[0].split('\t')
                    
                    # Extract gene-to-Ensembl mapping
                    added = 0
                    for i in range(1, len(lines)):
                        parts = lines[i].split('\t')
                        if len(parts) >= 2:
                            ensembl_id = parts[0]
                            gene_symbol = parts[1]
                            
                            if gene_symbol and ensembl_id and gene_symbol not in self.gene_mapping:
                                self.gene_mapping[gene_symbol] = ensembl_id
                                added += 1
                    
                    return added
            
            # If something went wrong
            print(f"Error fetching from BioMart: {response.status_code}")
            return 0
            
        except Exception as e:
            print(f"Error fetching from BioMart: {str(e)}")
            return 0