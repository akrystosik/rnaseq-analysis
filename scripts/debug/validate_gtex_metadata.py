#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/validate_gtex_metadata.py
#!/usr/bin/env python3
# Fixed GTEx Metadata Validation Script
# Usage: python3 validate_gtex_metadata.py

import os
import pandas as pd
import numpy as np
import gzip
import io
import re

# Set the metadata directory
METADATA_DIR = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/phenotype"

print("Starting GTEx metadata validation...")
print(f"Looking for files in: {METADATA_DIR}")

# List of compressed data files to analyze
subject_file = os.path.join(METADATA_DIR, "phs000424.v10.pht002740.v9.p2.GTEx_Subject.MULTI.txt.gz")
sample_file = os.path.join(METADATA_DIR, "phs000424.v10.pht002741.v10.p2.GTEx_Sample.MULTI.txt.gz")
subject_phenotypes_file = os.path.join(METADATA_DIR, "phs000424.v10.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz")
sample_attributes_file = os.path.join(METADATA_DIR, "phs000424.v10.pht002743.v10.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz")

# Check file existence
print("\nChecking data files:")
for file_path in [subject_file, sample_file, subject_phenotypes_file, sample_attributes_file]:
    if os.path.exists(file_path):
        file_size = os.path.getsize(file_path) / (1024 * 1024)  # Size in MB
        print(f"  ✓ {os.path.basename(file_path)} ({file_size:.2f} MB)")
    else:
        print(f"  ✗ {os.path.basename(file_path)} (Not found)")

# Function to load gzipped text file with explicit delimiter detection
def load_gzipped_file(file_path):
    try:
        # First read a few lines to analyze the file structure
        with gzip.open(file_path, 'rt') as f:
            # Skip potential header comments
            line = f.readline()
            while line.startswith('#'):
                line = f.readline()
                
            # Try to detect the delimiter
            if '\t' in line:
                delimiter = '\t'
            elif ',' in line:
                delimiter = ','
            else:
                # If can't detect, default to tab
                delimiter = '\t'
                
            # Add the first non-comment line back to the buffer
            buffer = io.StringIO(line)
            # Add the rest of the file
            buffer.write(f.read())
            buffer.seek(0)
            
            # Now read with the detected delimiter
            df = pd.read_csv(buffer, sep=delimiter, on_bad_lines='skip', comment='#')
            return df
    except Exception as e:
        # Try a more brute-force approach if the above fails
        try:
            print(f"  Trying alternative parsing method for {os.path.basename(file_path)}...")
            
            # Read file contents as text
            with gzip.open(file_path, 'rt') as f:
                content = f.readlines()
            
            # Skip comments and find header line
            header_line = None
            data_lines = []
            
            for i, line in enumerate(content):
                if line.startswith('#'):
                    continue
                if header_line is None:
                    header_line = line.strip()
                else:
                    data_lines.append(line.strip())
            
            if header_line:
                # Split by tab and create dataframe
                headers = re.split(r'\t|,', header_line)
                
                # Process data lines
                data = []
                for line in data_lines:
                    fields = re.split(r'\t|,', line)
                    # Ensure the row has the same number of fields as headers
                    if len(fields) >= len(headers):
                        data.append(fields[:len(headers)])
                    else:
                        # Pad with None if needed
                        data.append(fields + [None] * (len(headers) - len(fields)))
                
                # Create dataframe
                df = pd.DataFrame(data, columns=headers)
                return df
            
            return None
        except Exception as detailed_e:
            print(f"  Error reading {os.path.basename(file_path)}: {detailed_e}")
            return None

# Function to analyze and print dataframe information
def analyze_dataframe(df, name):
    print(f"\n=== {name} Analysis ===")
    if df is None or df.empty:
        print("  No data available")
        return None
    
    print(f"  Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    print(f"  Columns: {', '.join(df.columns[:10])}{'...' if len(df.columns) > 10 else ''}")
    
    return df

# Load subject data
print("\nLoading subject data...")
subject_df = load_gzipped_file(subject_file)
subject_df = analyze_dataframe(subject_df, "Subject Data")

# Load sample data
print("\nLoading sample data...")
sample_df = load_gzipped_file(sample_file)
sample_df = analyze_dataframe(sample_df, "Sample Data")

# Load subject phenotypes data
print("\nLoading subject phenotypes data...")
subject_phenotypes_df = load_gzipped_file(subject_phenotypes_file)
subject_phenotypes_df = analyze_dataframe(subject_phenotypes_df, "Subject Phenotypes Data")

# Load sample attributes data
print("\nLoading sample attributes data...")
sample_attributes_df = load_gzipped_file(sample_attributes_file)
sample_attributes_df = analyze_dataframe(sample_attributes_df, "Sample Attributes Data")

# Analyze donor demographics
print("\n=== Donor Demographics Analysis ===")
if subject_df is not None and not subject_df.empty:
    # Count unique donors
    if 'SUBJID' in subject_df.columns:
        unique_donors = subject_df['SUBJID'].nunique()
        print(f"  Total unique donors: {unique_donors}")
    
    # Sex distribution
    if 'SEX' in subject_df.columns:
        sex_mapping = {1: 'Male', 2: 'Female', '1': 'Male', '2': 'Female'}
        sex_counts = subject_df['SEX'].map(sex_mapping).value_counts(normalize=True) * 100
        
        print(f"  Sex distribution:")
        for sex, percentage in sex_counts.items():
            print(f"    - {sex}: {percentage:.1f}%")
    else:
        print("  SEX column not found")
else:
    print("  No subject data available")

# Try to get ethnicity data from subject phenotypes
if subject_phenotypes_df is not None and not subject_phenotypes_df.empty:
    race_cols = [col for col in subject_phenotypes_df.columns if 'RACE' in col.upper()]
    if race_cols:
        print("\n  Ethnic distribution:")
        for col in race_cols:
            race_counts = subject_phenotypes_df[col].value_counts(normalize=True) * 100
            for race, percentage in race_counts.items():
                print(f"    - {race}: {percentage:.1f}%")
    else:
        print("\n  No race/ethnicity columns found in subject phenotypes")
else:
    print("\n  No subject phenotypes data available")

# Analyze sample data
print("\n=== Sample Analysis ===")
if sample_df is not None and not sample_df.empty:
    # Count samples
    total_samples = sample_df.shape[0]
    print(f"  Total samples: {total_samples}")
else:
    print("  No sample data available")

# Get tissue information from sample attributes
if sample_attributes_df is not None and not sample_attributes_df.empty:
    # Look for tissue columns
    tissue_cols = [col for col in sample_attributes_df.columns if 'SMTS' in col or 'SMTSD' in col]
    
    if tissue_cols:
        print("\n  Tissue information:")
        for col in tissue_cols:
            unique_tissues = sample_attributes_df[col].nunique()
            print(f"    - {col}: {unique_tissues} unique values")
            
            # Top tissues by sample count
            top_tissues = sample_attributes_df[col].value_counts().head(10)
            print(f"    Top 10 tissues by sample count ({col}):")
            for tissue, count in top_tissues.items():
                print(f"      - {tissue}: {count} samples")
    else:
        print("  No tissue columns found in sample attributes")
else:
    print("  No sample attributes data available")

# Try to manually count columns in files
print("\n=== Manual File Analysis ===")
for file_path, name in [
    (subject_file, "Subject File"),
    (sample_file, "Sample File"),
    (subject_phenotypes_file, "Subject Phenotypes File"),
    (sample_attributes_file, "Sample Attributes File")
]:
    try:
        with gzip.open(file_path, 'rt') as f:
            lines = []
            for i, line in enumerate(f):
                if i < 10:  # Only read first 10 lines
                    lines.append(line.strip())
                else:
                    break
            
            print(f"\n  {name} - First few lines:")
            for i, line in enumerate(lines):
                if i < 5:  # Only show first 5 lines
                    if len(line) > 100:
                        # Truncate long lines
                        print(f"    {line[:100]}...")
                    else:
                        print(f"    {line}")
            
            # Count rows
            with gzip.open(file_path, 'rt') as f:
                row_count = sum(1 for line in f if not line.startswith('#'))
            print(f"  Total rows in {name}: {row_count}")
    except Exception as e:
        print(f"  Error processing {name}: {e}")

# Generate summary
print("\n=== Summary of Findings vs. Document Claims ===")
print("  Document Claim | Validation Result")
print("  -------------- | -----------------")

# Compare donor count
document_donors = 961
if 'unique_donors' in locals() and unique_donors is not None:
    validated_donors = unique_donors
else:
    # Try to get donor count from manual file analysis
    validated_donors = "See manual file analysis"
print(f"  Donors: {document_donors} | {validated_donors}")

# Compare tissue types
document_tissues = 54
validated_tissues = "See tissue breakdown above"
print(f"  Tissue types: {document_tissues} | {validated_tissues}")

# Compare sample count
document_samples = 19788
if 'total_samples' in locals() and total_samples is not None:
    validated_samples = total_samples
else:
    # Try to get sample count from manual file analysis
    validated_samples = "See manual file analysis"
print(f"  Samples: {document_samples} | {validated_samples}")

print("\nValidation complete!")