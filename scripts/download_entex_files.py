#!/usr/bin/env python3
"""
Script to download ENTEx gene quantification files from ENCODE
and organize them in the appropriate directory structure.
"""

import os
import requests
import concurrent.futures
import json
import re
from pathlib import Path
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Download ENTEx gene quantification files from ENCODE')
    parser.add_argument('--input', '-i', required=True, help='Input file containing URLs of files to download')
    parser.add_argument('--output-dir', '-o', required=True, help='Output directory for downloaded files')
    parser.add_argument('--threads', '-t', type=int, default=4, help='Number of download threads')
    parser.add_argument('--metadata-file', '-m', help='Path to output metadata JSON file')
    return parser.parse_args()

def download_file(url, output_dir):
    """Download a single file from a URL and save it to the output directory."""
    # Extract filename from URL
    filename = url.split('/')[-1]
    filepath = os.path.join(output_dir, filename)
    
    # Skip if file already exists
    if os.path.exists(filepath):
        print(f"File {filename} already exists, skipping")
        return filepath, True
    
    try:
        print(f"Downloading {filename} from {url}")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        # Create directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Save file
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print(f"Successfully downloaded {filename}")
        return filepath, True
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
        return url, False

def fetch_experiment_metadata(file_ids):
    """Fetch metadata for the given file IDs from the ENCODE API."""
    metadata = {}
    
    for file_id in file_ids:
        # Extract file ID from URL or filename
        if isinstance(file_id, str) and 'ENCFF' in file_id:
            # Extract ENCFF ID from string
            match = re.search(r'(ENCFF\w+)', file_id)
            if match:
                file_id = match.group(1)
        
        url = f"https://www.encodeproject.org/files/{file_id}/?format=json"
        try:
            print(f"Fetching metadata for {file_id}")
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            # Extract relevant metadata
            experiment_id = data.get('dataset', '').split('/')[-2] if data.get('dataset') else None
            
            if experiment_id:
                exp_url = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
                exp_response = requests.get(exp_url)
                exp_response.raise_for_status()
                exp_data = exp_response.json()
                
                # Extract tissue, donor, and other metadata
                biosample = exp_data.get('biosample_ontology', {})
                tissue = biosample.get('term_name', 'unknown_tissue')
                donor = None
                for rep in exp_data.get('replicates', []):
                    if rep.get('library', {}).get('biosample', {}).get('donor', {}).get('accession'):
                        donor = rep.get('library', {}).get('biosample', {}).get('donor', {}).get('accession')
                        break
                
                metadata[file_id] = {
                    'file_id': file_id,
                    'experiment_id': experiment_id,
                    'tissue': tissue,
                    'donor': donor,
                    'assay_type': exp_data.get('assay_title', 'unknown_assay'),
                    'file_format': data.get('file_format', 'unknown_format'),
                    'output_type': data.get('output_type', 'unknown_output_type'),
                    'biological_replicates': data.get('biological_replicates', []),
                    'technical_replicates': data.get('technical_replicates', []),
                    'genome_assembly': data.get('assembly', 'unknown_assembly')
                }
                print(f"Metadata fetched for {file_id}: {metadata[file_id]}")
            else:
                print(f"Could not find experiment for file {file_id}")
        
        except requests.exceptions.RequestException as e:
            print(f"Error fetching metadata for {file_id}: {e}")
    
    return metadata

def main():
    args = parse_args()
    
    # Read URLs from input file
    with open(args.input, 'r') as f:
        urls = [line.strip() for line in f if line.strip() and line.strip().startswith('http')]
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Download files using thread pool
    print(f"Downloading {len(urls)} files with {args.threads} threads")
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_url = {executor.submit(download_file, url, args.output_dir): url for url in urls}
        for future in concurrent.futures.as_completed(future_to_url):
            url = future_to_url[future]
            try:
                filepath, success = future.result()
                results.append((filepath, success))
            except Exception as e:
                print(f"Error processing {url}: {e}")
                results.append((url, False))
    
    # Count successes and failures
    successes = [r for r in results if r[1]]
    failures = [r for r in results if not r[1]]
    
    print(f"Download summary: {len(successes)} successes, {len(failures)} failures")
    
    # Generate list of file IDs from successful downloads
    file_ids = []
    for filepath, success in successes:
        if success:
            if isinstance(filepath, str) and os.path.exists(filepath):
                filename = os.path.basename(filepath)
                if 'ENCFF' in filename:
                    match = re.search(r'(ENCFF\w+)', filename)
                    if match:
                        file_ids.append(match.group(1))
    
    # Fetch metadata for downloaded files
    if file_ids:
        print(f"Fetching metadata for {len(file_ids)} files")
        metadata = fetch_experiment_metadata(file_ids)
        
        # Save metadata to file if specified
        if args.metadata_file:
            print(f"Saving metadata to {args.metadata_file}")
            with open(args.metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            print(f"Metadata saved to {args.metadata_file}")
    
    print("Processing complete!")

if __name__ == "__main__":
    main()