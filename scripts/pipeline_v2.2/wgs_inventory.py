#!/usr/bin/env python3
"""
WGS Data Inventory Script
Systematically catalog all WGS sample and subject IDs across datasets
"""

import os
import re
import csv
from pathlib import Path
import gzip

def get_cell_line_samples():
    """Extract ENCODE cell line samples from BAM and VCF files"""
    samples = []
    
    # BAM files
    bam_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/bam")
    if bam_dir.exists():
        for cell_line_dir in bam_dir.iterdir():
            if cell_line_dir.is_dir():
                cell_line = cell_line_dir.name
                for bam_file in cell_line_dir.glob("*.bam"):
                    sample_id = bam_file.stem.replace(".marked_duplicates", "")
                    samples.append({
                        'sample_id': sample_id,
                        'subject_id': cell_line,
                        'dataset': 'ENCODE',
                        'data_type': 'germline',
                        'file_types_available': 'BAM'
                    })
    
    # DeepVariant VCF files
    deepvariant_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/deepvariant")
    if deepvariant_dir.exists():
        for vcf_file in deepvariant_dir.glob("*.deepvariant.vcf.gz"):
            sample_id = vcf_file.name.split('.')[0]
            cell_line = sample_id.split('_')[0]
            # Find existing sample or create new
            existing = next((s for s in samples if s['sample_id'] == sample_id), None)
            if existing:
                existing['file_types_available'] += ',VCF_germline'
            else:
                samples.append({
                    'sample_id': sample_id,
                    'subject_id': cell_line,
                    'dataset': 'ENCODE',
                    'data_type': 'germline',
                    'file_types_available': 'VCF_germline'
                })
    
    # DeepSomatic VCF files
    deepsomatic_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/deepsomatic")
    if deepsomatic_dir.exists():
        for vcf_file in deepsomatic_dir.glob("*.deepsomatic.vcf.gz"):
            sample_id = vcf_file.name.split('.')[0]
            cell_line = sample_id.split('_')[0]
            # Find existing sample or create new
            existing = next((s for s in samples if s['sample_id'] == sample_id), None)
            if existing:
                existing['file_types_available'] += ',VCF_somatic'
                if existing['data_type'] == 'germline':
                    existing['data_type'] = 'germline+somatic'
            else:
                samples.append({
                    'sample_id': sample_id,
                    'subject_id': cell_line,
                    'dataset': 'ENCODE',
                    'data_type': 'somatic',
                    'file_types_available': 'VCF_somatic'
                })
    
    # ENCODE lifted VCF files (legacy/reference data)
    lifted_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/encode_lifted")
    if lifted_dir.exists():
        for vcf_file in lifted_dir.glob("ENCFF*.hg38.sorted.vcf.gz"):
            sample_id = vcf_file.name.split('.')[0]
            # These are ENCODE file IDs, treat as separate reference samples
            samples.append({
                'sample_id': sample_id,
                'subject_id': f"ENCODE_ref_{sample_id}",
                'dataset': 'ENCODE_legacy',
                'data_type': 'germline',
                'file_types_available': 'VCF_lifted_hg38'
            })
    
    return samples

def get_gtex_samples():
    """Extract GTEx samples from summary files and VCF"""
    samples = []
    
    # Read GTEx samples from summary files
    gtex_summary_file = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/gtex_counts/gtex_summary_complete.tsv")
    if gtex_summary_file.exists():
        with open(gtex_summary_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['Sample'].startswith('GTEX-') and int(row['Total']) > 0:
                    samples.append({
                        'sample_id': row['Sample'],
                        'subject_id': row['Sample'],  # GTEx uses same ID for sample and subject
                        'dataset': 'GTEx',
                        'data_type': 'germline',
                        'file_types_available': 'VCF_germline'
                    })
    
    # Also check all_samples.txt for additional samples
    all_samples_file = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/gtex_counts/all_samples.txt")
    if all_samples_file.exists():
        with open(all_samples_file, 'r') as f:
            for line in f:
                sample_id = line.strip()
                if sample_id.startswith('GTEX-'):
                    # Check if already in samples
                    if not any(s['sample_id'] == sample_id for s in samples):
                        samples.append({
                            'sample_id': sample_id,
                            'subject_id': sample_id,
                            'dataset': 'GTEx',
                            'data_type': 'germline',
                            'file_types_available': 'VCF_germline'
                        })
    
    return samples

def analyze_identifier_patterns(samples):
    """Analyze naming patterns and conventions"""
    patterns = {
        'ENCODE': {'cell_lines': set(), 'pair_formats': set()},
        'GTEx': {'id_format': set(), 'prefix_counts': {}}
    }
    
    for sample in samples:
        if sample['dataset'] == 'ENCODE':
            # Extract cell line name
            cell_line = sample['subject_id']
            patterns['ENCODE']['cell_lines'].add(cell_line)
            
            # Check for pair format
            if '_pair' in sample['sample_id']:
                pair_match = re.search(r'_pair(\d+)', sample['sample_id'])
                if pair_match:
                    patterns['ENCODE']['pair_formats'].add(f"pair{pair_match.group(1)}")
        
        elif sample['dataset'] == 'GTEx':
            # Analyze GTEx ID format
            if sample['sample_id'].startswith('GTEX-'):
                suffix = sample['sample_id'][5:]  # Remove 'GTEX-' prefix
                patterns['GTEx']['id_format'].add(len(suffix))
                
                # Count prefix patterns
                prefix = suffix[:2] if len(suffix) >= 2 else suffix
                patterns['GTEx']['prefix_counts'][prefix] = patterns['GTEx']['prefix_counts'].get(prefix, 0) + 1
    
    return patterns

def main():
    print("Starting WGS data inventory...")
    
    # Collect samples from all datasets
    all_samples = []
    
    # ENCODE cell lines
    print("Collecting ENCODE cell line samples...")
    encode_samples = get_cell_line_samples()
    all_samples.extend(encode_samples)
    print(f"Found {len(encode_samples)} ENCODE samples")
    
    # GTEx
    print("Collecting GTEx samples...")
    gtex_samples = get_gtex_samples()
    all_samples.extend(gtex_samples)
    print(f"Found {len(gtex_samples)} GTEx samples")
    
    # Analyze patterns
    print("Analyzing identifier patterns...")
    patterns = analyze_identifier_patterns(all_samples)
    
    # Save results to CSV
    output_file = "wgs_sample_inventory.csv"
    with open(output_file, 'w', newline='') as f:
        fieldnames = ['sample_id', 'subject_id', 'dataset', 'data_type', 'file_types_available']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_samples)
    
    print(f"\nInventory saved to {output_file}")
    
    # Print summary
    print("\n" + "="*60)
    print("WGS DATA INVENTORY SUMMARY")
    print("="*60)
    
    # Dataset counts
    dataset_counts = {}
    data_type_counts = {}
    for sample in all_samples:
        dataset = sample['dataset']
        data_type = sample['data_type']
        dataset_counts[dataset] = dataset_counts.get(dataset, 0) + 1
        data_type_counts[data_type] = data_type_counts.get(data_type, 0) + 1
    
    print(f"\nTotal WGS samples: {len(all_samples)}")
    print("\nBy Dataset:")
    for dataset, count in sorted(dataset_counts.items()):
        print(f"  {dataset}: {count} samples")
    
    print("\nBy Data Type:")
    for data_type, count in sorted(data_type_counts.items()):
        print(f"  {data_type}: {count} samples")
    
    # ENCODE patterns
    print(f"\nENCODE Cell Lines ({len(patterns['ENCODE']['cell_lines'])}):")
    for cell_line in sorted(patterns['ENCODE']['cell_lines']):
        count = sum(1 for s in all_samples if s['dataset'] == 'ENCODE' and s['subject_id'] == cell_line)
        print(f"  {cell_line}: {count} samples")
    
    print(f"\nENCODE Pair Formats: {sorted(patterns['ENCODE']['pair_formats'])}")
    
    # GTEx patterns
    print(f"\nGTEx ID Format:")
    print(f"  ID lengths: {sorted(patterns['GTEx']['id_format'])}")
    print(f"  Common prefixes (top 10):")
    top_prefixes = sorted(patterns['GTEx']['prefix_counts'].items(), key=lambda x: x[1], reverse=True)[:10]
    for prefix, count in top_prefixes:
        print(f"    {prefix}*: {count} samples")
    
    # File type availability
    print(f"\nFile Types Available:")
    file_types = {}
    for sample in all_samples:
        types = sample['file_types_available']
        file_types[types] = file_types.get(types, 0) + 1
    
    for file_type, count in sorted(file_types.items()):
        print(f"  {file_type}: {count} samples")

if __name__ == "__main__":
    main()