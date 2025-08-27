#!/usr/bin/env python3
"""
Gene Mapping Builder Script
--------------------------

This script builds a comprehensive gene mapping database from:
1. Existing TSV files
2. BioMart queries
3. Additional online databases

Having a good gene mapping is critical for accurate gene identification.
"""

import argparse
import json
import os
from pathlib import Path
import sys
import glob
import requests
import pandas as pd

# Add the project root to the path to import modules
script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
project_root = script_dir.parent
sys.path.append(str(project_root))

from core.dataset_manager import DatasetManager
from core.gene_identifier import GeneIdentifier

def scan_tsv_directory(directory, gene_identifier):
    """Scan directory for TSV files to build mapping"""
    print(f"Scanning directory: {directory}")
    
    # Find all TSV files recursively
    tsv_files = glob.glob(str(directory) + "/**/*.tsv", recursive=True)
    
    print(f"Found {len(tsv_files)} TSV files")
    
    added_total = 0
    processed_files = 0
    
    # Process each file
    for file_path in tsv_files:
        try:
            added = gene_identifier.build_mapping_from_file(file_path)
            added_total += added
            processed_files += 1
            
            if processed_files % 10 == 0:
                print(f"Processed {processed_files}/{len(tsv_files)} files, found {added_total} mappings")
                
            # If we've added enough genes, stop
            if added_total > 30000:
                break
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    print(f"Completed processing {processed_files} files")
    print(f"Added {added_total} gene mappings")
    
    return added_total

def fetch_from_biomart(gene_identifier, additional_genes=None):
    """Fetch gene mappings from BioMart"""
    print("Fetching gene mappings from BioMart...")
    
    # Start with housekeeping genes that don't already have mappings
    missing_genes = [
        gene for gene in gene_identifier.housekeeping_genes 
        if gene not in gene_identifier.gene_mapping
    ]
    
    # Add additional genes if provided
    if additional_genes:
        missing_genes.extend([g for g in additional_genes if g not in gene_identifier.gene_mapping])
    
    if not missing_genes:
        print("No missing genes to fetch")
        return 0
    
    print(f"Fetching mappings for {len(missing_genes)} genes")
    
    # Fetch in batches to avoid timeout
    batch_size = 100
    total_added = 0
    
    for i in range(0, len(missing_genes), batch_size):
        batch = missing_genes[i:i+batch_size]
        added = gene_identifier.fetch_mapping_from_biomart(batch)
        total_added += added
        print(f"Batch {i//batch_size + 1}: Added {added} mappings")
    
    print(f"Added {total_added} mappings from BioMart")
    return total_added

def fetch_from_mygene(gene_identifier, additional_genes=None):
    """Fetch gene mappings from MyGene.info"""
    try:
        import mygene
    except ImportError:
        print("mygene package not installed. Skipping MyGene.info query.")
        print("Install with: pip install mygene")
        return 0
    
    print("Fetching gene mappings from MyGene.info...")
    
    # Get list of genes that don't have mappings yet
    missing_genes = [
        gene for gene in gene_identifier.housekeeping_genes 
        if gene not in gene_identifier.gene_mapping
    ]
    
    # Add additional genes if provided
    if additional_genes:
        missing_genes.extend([g for g in additional_genes if g not in gene_identifier.gene_mapping])
    
    if not missing_genes:
        print("No missing genes to fetch")
        return 0
    
    print(f"Fetching mappings for {len(missing_genes)} genes")
    
    # Initialize MyGene client
    mg = mygene.MyGeneInfo()
    
    # Query MyGene.info
    try:
        results = mg.querymany(missing_genes, scopes='symbol', fields='ensembl.gene')
        
        # Process results
        added = 0
        for result in results:
            if 'ensembl' in result and 'gene' in result['ensembl']:
                gene_symbol = result['query']
                
                # Handle multiple ensembl IDs
                if isinstance(result['ensembl']['gene'], list):
                    ensembl_id = result['ensembl']['gene'][0]
                else:
                    ensembl_id = result['ensembl']['gene']
                
                # Add to mapping
                if gene_symbol not in gene_identifier.gene_mapping:
                    gene_identifier.gene_mapping[gene_symbol] = ensembl_id
                    added += 1
        
        print(f"Added {added} mappings from MyGene.info")
        return added
        
    except Exception as e:
        print(f"Error querying MyGene.info: {e}")
        return 0

def main():
    """Build comprehensive gene mapping"""
    parser = argparse.ArgumentParser(description='Build gene mapping database')
    
    # Input options
    parser.add_argument('--scan-dir', help='Directory to scan for TSV files')
    parser.add_argument('--additional-genes', help='File with additional genes to map (one per line)')
    
    # Source options
    parser.add_argument('--use-biomart', action='store_true', help='Query BioMart for missing genes')
    parser.add_argument('--use-mygene', action='store_true', help='Query MyGene.info for missing genes')
    
    # Output options
    parser.add_argument('--output-file', help='Custom output file for gene mapping')
    parser.add_argument('--config', help='Path to custom config file')
    
    args = parser.parse_args()
    
    print("\n===== Building Gene Mapping Database =====")
    
    # Set up components
    gene_identifier = GeneIdentifier(args.config)
    dataset_manager = DatasetManager(args.config)
    
    # Load initial mapping
    initial_count = len(gene_identifier.gene_mapping)
    print(f"Starting with {initial_count} gene mappings")
    
    # Load additional genes if provided
    additional_genes = []
    if args.additional_genes and os.path.exists(args.additional_genes):
        with open(args.additional_genes, 'r') as f:
            additional_genes = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(additional_genes)} additional genes to map")
    
    # Scan directory for TSV files
    if args.scan_dir:
        scan_path = Path(args.scan_dir)
        if scan_path.exists() and scan_path.is_dir():
            scan_tsv_directory(scan_path, gene_identifier)
        else:
            print(f"Error: {args.scan_dir} is not a valid directory")
    else:
        # Default to scanning the raw data directory
        scan_tsv_directory(dataset_manager.raw_data_dir, gene_identifier)
    
    # Query BioMart if requested
    if args.use_biomart:
        fetch_from_biomart(gene_identifier, additional_genes)
    
    # Query MyGene.info if requested
    if args.use_mygene:
        fetch_from_mygene(gene_identifier, additional_genes)
    
    # Save mapping
    if args.output_file:
        output_file = Path(args.output_file)
    else:
        output_file = Path(dataset_manager.metadata_dir) / "gene_mapping.json"
    
    with open(output_file, 'w') as f:
        json.dump(gene_identifier.gene_mapping, f, indent=2)
    
    final_count = len(gene_identifier.gene_mapping)
    new_count = final_count - initial_count
    
    print(f"\nGene mapping database built successfully:")
    print(f"  Initial mappings: {initial_count}")
    print(f"  New mappings added: {new_count}")
    print(f"  Total mappings: {final_count}")
    print(f"  Saved to: {output_file}")
    
    # Check coverage of housekeeping genes
    hk_mapped = sum(1 for gene in gene_identifier.housekeeping_genes if gene in gene_identifier.gene_mapping)
    print(f"\nHousekeeping gene coverage: {hk_mapped}/{len(gene_identifier.housekeeping_genes)} ({hk_mapped/len(gene_identifier.housekeeping_genes)*100:.1f}%)")
    
    # Check coverage of additional genes if provided
    if additional_genes:
        add_mapped = sum(1 for gene in additional_genes if gene in gene_identifier.gene_mapping)
        print(f"Additional gene coverage: {add_mapped}/{len(additional_genes)} ({add_mapped/len(additional_genes)*100:.1f}%)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())