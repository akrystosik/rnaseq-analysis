#!/usr/bin/env python3
"""
Create Affymetrix U219 to Ensembl ID Mapping Using GEO SOFT Format File
With auto-installation of dependencies

Usage:
  python u219_mapping_self_install.py

Output:
  Creates a CSV file with probe ID to Ensembl ID mappings
"""

import sys
import subprocess
import os

# Check and install required packages
required_packages = ['mygene', 'pandas', 'requests']
for package in required_packages:
    try:
        __import__(package)
    except ImportError:
        print(f"Installing required package: {package}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Now import after installation
import requests
import pandas as pd
import re
import mygene

def download_soft_file(geo_id="GPL13667"):
    """
    Download the SOFT format file for a GEO platform.
    
    Parameters:
    -----------
    geo_id : str
        GEO platform ID (e.g., GPL13667 for Affymetrix U219)
    
    Returns:
    --------
    str
        Content of the SOFT file
    """
    print(f"Downloading SOFT file for {geo_id}...")
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}&targ=self&form=text&view=full"
    
    try:
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Failed to download SOFT file: HTTP {response.status_code}")
            return None
            
        print(f"Downloaded {len(response.text)} bytes")
        return response.text
    except Exception as e:
        print(f"Error downloading SOFT file: {e}")
        return None

def parse_soft_file(soft_content):
    """
    Parse a GEO SOFT format file to extract probe mappings.
    
    Parameters:
    -----------
    soft_content : str
        Content of the SOFT file
    
    Returns:
    --------
    dict, dict, dict
        Dictionaries mapping probe IDs to Ensembl IDs, Entrez IDs, and gene symbols
    """
    print("Parsing SOFT file...")
    
    # Initialize dictionaries
    probe_to_ensembl = {}
    probe_to_entrez = {}
    probe_to_symbol = {}
    
    # Find the table section
    table_start = False
    column_indices = {}
    
    lines = soft_content.split('\n')
    for i, line in enumerate(lines):
        # Find table header
        if line.startswith('!platform_table_begin'):
            table_start = True
            # The next line should be the header
            if i+1 < len(lines):
                header = lines[i+1].split('\t')
                # Find column indices
                for j, col in enumerate(header):
                    if col == 'ID':
                        column_indices['probe_id'] = j
                    elif col == 'Ensembl':
                        column_indices['ensembl'] = j
                    elif col == 'Entrez Gene':
                        column_indices['entrez'] = j
                    elif col == 'Gene Symbol':
                        column_indices['symbol'] = j
            continue
            
        # End of table
        if line.startswith('!platform_table_end'):
            break
            
        # Process data lines
        if table_start and not line.startswith('!') and not line.startswith('#') and line.strip():
            fields = line.split('\t')
            
            if len(fields) <= max(column_indices.values()):
                continue
                
            # Get probe ID
            probe_id = fields[column_indices['probe_id']]
            
            # Get Ensembl ID if available
            if 'ensembl' in column_indices:
                ensembl_column = fields[column_indices['ensembl']]
                # Extract Ensembl gene ID if present
                if ensembl_column and ensembl_column != '---':
                    # Look for ENSG pattern
                    ensembl_match = re.search(r'(ENSG\d+)', ensembl_column)
                    if ensembl_match:
                        ensembl_id = ensembl_match.group(1)
                        probe_to_ensembl[probe_id] = ensembl_id
            
            # Get Entrez ID if available
            if 'entrez' in column_indices:
                entrez_id = fields[column_indices['entrez']]
                if entrez_id and entrez_id != '---' and entrez_id.strip():
                    probe_to_entrez[probe_id] = entrez_id
            
            # Get gene symbol if available
            if 'symbol' in column_indices:
                symbol = fields[column_indices['symbol']]
                if symbol and symbol != '---' and symbol.strip():
                    probe_to_symbol[probe_id] = symbol
    
    print(f"Parsed {len(probe_to_ensembl)} direct Ensembl mappings")
    print(f"Found {len(probe_to_entrez)} probes with Entrez IDs")
    print(f"Found {len(probe_to_symbol)} probes with gene symbols")
    
    return probe_to_ensembl, probe_to_entrez, probe_to_symbol

def map_entrez_to_ensembl_mygene(entrez_ids):
    """
    Map Entrez Gene IDs to Ensembl IDs using MyGene.
    
    Parameters:
    -----------
    entrez_ids : list
        List of Entrez Gene IDs
    
    Returns:
    --------
    dict
        Dictionary mapping Entrez IDs to Ensembl IDs
    """
    print(f"Mapping {len(entrez_ids)} Entrez IDs to Ensembl IDs using MyGene...")
    
    entrez_to_ensembl = {}
    
    try:
        # Initialize MyGene client
        mg = mygene.MyGeneInfo()
        
        # Process in batches to avoid timeout
        batch_size = 1000
        for i in range(0, len(entrez_ids), batch_size):
            batch = entrez_ids[i:i+batch_size]
            print(f"Processing batch {i//batch_size + 1}/{(len(entrez_ids) + batch_size - 1)//batch_size}")
            
            # Query for batch of Entrez IDs
            result = mg.querymany(batch, scopes='entrezgene', fields='ensembl.gene', species='human')
            
            # Extract mappings
            for item in result:
                if 'ensembl' in item and not item.get('notfound', False):
                    entrez_id = str(item['query'])
                    if isinstance(item['ensembl'], dict):
                        ensembl_id = item['ensembl'].get('gene')
                    elif isinstance(item['ensembl'], list):
                        ensembl_id = item['ensembl'][0].get('gene')
                    else:
                        continue
                        
                    if ensembl_id:
                        entrez_to_ensembl[entrez_id] = ensembl_id
        
        print(f"Mapped {len(entrez_to_ensembl)} Entrez IDs to Ensembl IDs")
    
    except Exception as e:
        print(f"Error using MyGene for Entrez-to-Ensembl mapping: {e}")
    
    return entrez_to_ensembl

def map_symbols_to_ensembl_mygene(gene_symbols):
    """
    Map gene symbols to Ensembl IDs using MyGene.
    
    Parameters:
    -----------
    gene_symbols : list
        List of gene symbols
    
    Returns:
    --------
    dict
        Dictionary mapping gene symbols to Ensembl IDs
    """
    print(f"Mapping {len(gene_symbols)} gene symbols to Ensembl IDs using MyGene...")
    
    symbol_to_ensembl = {}
    
    try:
        # Initialize MyGene client
        mg = mygene.MyGeneInfo()
        
        # Process in batches to avoid timeout
        batch_size = 1000
        for i in range(0, len(gene_symbols), batch_size):
            batch = gene_symbols[i:i+batch_size]
            print(f"Processing batch {i//batch_size + 1}/{(len(gene_symbols) + batch_size - 1)//batch_size}")
            
            # Query for batch of gene symbols
            result = mg.querymany(batch, scopes='symbol', fields='ensembl.gene', species='human')
            
            # Extract mappings
            for item in result:
                if 'ensembl' in item and not item.get('notfound', False):
                    symbol = str(item['query'])
                    if isinstance(item['ensembl'], dict):
                        ensembl_id = item['ensembl'].get('gene')
                    elif isinstance(item['ensembl'], list):
                        ensembl_id = item['ensembl'][0].get('gene')
                    else:
                        continue
                        
                    if ensembl_id:
                        symbol_to_ensembl[symbol] = ensembl_id
        
        print(f"Mapped {len(symbol_to_ensembl)} gene symbols to Ensembl IDs")
    
    except Exception as e:
        print(f"Error using MyGene for Symbol-to-Ensembl mapping: {e}")
    
    return symbol_to_ensembl

def create_mapping_file():
    """
    Create a comprehensive mapping file from Affymetrix U219 probe IDs to Ensembl gene IDs.
    """
    # 1. Download and parse SOFT file
    soft_content = download_soft_file()
    
    if not soft_content:
        print("Failed to download SOFT file. Exiting.")
        return
    
    # 2. Parse the SOFT file
    probe_to_ensembl, probe_to_entrez, probe_to_symbol = parse_soft_file(soft_content)
    
    # 3. For probes with Entrez IDs but no Ensembl IDs, map via Entrez using MyGene
    entrez_ids_to_map = set()
    probes_needing_entrez_mapping = {}
    
    for probe_id, entrez_id in probe_to_entrez.items():
        if probe_id not in probe_to_ensembl:
            entrez_ids_to_map.add(entrez_id)
            probes_needing_entrez_mapping[probe_id] = entrez_id
    
    if entrez_ids_to_map:
        print(f"Found {len(entrez_ids_to_map)} unique Entrez IDs for probes without direct Ensembl mapping")
        entrez_to_ensembl = map_entrez_to_ensembl_mygene(list(entrez_ids_to_map))
        
        # Add the mappings
        for probe_id, entrez_id in probes_needing_entrez_mapping.items():
            if entrez_id in entrez_to_ensembl:
                probe_to_ensembl[probe_id] = entrez_to_ensembl[entrez_id]
    
    # 4. For probes with gene symbols but no Ensembl IDs, map via symbols using MyGene
    symbols_to_map = set()
    probes_needing_symbol_mapping = {}
    
    for probe_id, symbol in probe_to_symbol.items():
        if probe_id not in probe_to_ensembl:
            symbols_to_map.add(symbol)
            probes_needing_symbol_mapping[probe_id] = symbol
    
    if symbols_to_map:
        print(f"Found {len(symbols_to_map)} unique gene symbols for probes without direct Ensembl mapping")
        symbol_to_ensembl = map_symbols_to_ensembl_mygene(list(symbols_to_map))
        
        # Add the mappings
        for probe_id, symbol in probes_needing_symbol_mapping.items():
            if symbol in symbol_to_ensembl:
                probe_to_ensembl[probe_id] = symbol_to_ensembl[symbol]
    
    # 5. Create the final mapping file
    mapping_df = pd.DataFrame({
        'probe_id': list(probe_to_ensembl.keys()),
        'ensembl_id': list(probe_to_ensembl.values())
    })
    
    # Save to CSV
    output_file = "u219_probe_to_ensembl.csv"
    mapping_df.to_csv(output_file, index=False)
    
    print(f"Created mapping file with {len(mapping_df)} probe-to-Ensembl mappings")
    print(f"Saved to {output_file}")
    
    return output_file

if __name__ == "__main__":
    mapping_file = create_mapping_file()
    print(f"Mapping file created: {mapping_file}")