#!/usr/bin/env python3
"""
ADNI Metadata Processor - Revised with Visit Date

Extracts demographic information from ADNI dataset and creates a standardized 
JSON metadata file for integration with the RNA-seq pipeline.
"""

import pandas as pd
import json
import os
from pathlib import Path
import datetime
from dateutil.relativedelta import relativedelta

# Paths
adni_metadata_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadataADNI/subject_demographics/PTDEMOG_25Apr2025.csv"
output_json_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/adni_metadata.json"

# Load the ADNI demographics data
print(f"Loading ADNI metadata from {adni_metadata_path}")
adni_demo = pd.read_csv(adni_metadata_path)

# Print basic info
print(f"Loaded {len(adni_demo)} records")
print(f"Columns: {adni_demo.columns.tolist()}")

# Define demographic field mappings based on data dictionary
sex_mapping = {
    "M": "male",
    "F": "female",
    "m": "male",
    "f": "female",
    "Male": "male",
    "Female": "female",
    "1": "male",
    "2": "female",
    1: "male",
    2: "female"
}

ethnicity_mapping = {
    1: "Hispanic or Latino",
    2: "Not Hispanic or Latino",
    3: "Unknown",
    "1": "Hispanic or Latino",
    "2": "Not Hispanic or Latino",
    "3": "Unknown"
}

race_mapping = {
    1: "American Indian or Alaskan Native",
    2: "Asian",
    3: "Native Hawaiian or Other Pacific Islander",
    4: "Black or African American",
    5: "White",
    6: "More than one race",
    7: "Unknown",
    "1": "American Indian or Alaskan Native",
    "2": "Asian",
    "3": "Native Hawaiian or Other Pacific Islander",
    "4": "Black or African American",
    "5": "White",
    "6": "More than one race",
    "7": "Unknown"
}

# Create a mapping from RID to demographic information
subject_demographics = {}

# Track which columns we find
found_columns = {
    'sex': False,
    'ethnicity': False,
    'race': False,
    'age': False
}

# First, get the first visit date for each subject
subject_first_visits = {}
for _, row in adni_demo.iterrows():
    rid = row.get('RID')
    if pd.isna(rid):
        continue
    
    rid = int(rid)
    
    # Check for visit date
    if 'VISDATE' in row and not pd.isna(row['VISDATE']):
        try:
            visit_date = pd.to_datetime(row['VISDATE'])
            
            # Update first visit date if this is the first or earlier than current
            if rid not in subject_first_visits or visit_date < subject_first_visits[rid]:
                subject_first_visits[rid] = visit_date
        except:
            pass  # Skip if date conversion fails

# Now process the demographic data with visit dates
for _, row in adni_demo.iterrows():
    rid = row.get('RID')
    if pd.isna(rid):
        continue
    
    rid = int(rid)
    adni_id = f"ADNI_{rid}"
    
    # Skip if we already processed this subject
    if adni_id in subject_demographics:
        continue
    
    # Extract demographic information
    demographics = {}
    
    # Sex - Check for PTGENDER
    if 'PTGENDER' in row and not pd.isna(row['PTGENDER']):
        sex_code = row['PTGENDER']
        demographics['sex'] = sex_mapping.get(sex_code, "unknown")
        found_columns['sex'] = True
    
    # Age - Calculate from birth date and visit date if available
    if 'PTDOB' in row and not pd.isna(row['PTDOB']) and rid in subject_first_visits:
        try:
            birth_date = pd.to_datetime(row['PTDOB'])
            visit_date = subject_first_visits[rid]
            
            # Calculate age in years
            age_delta = relativedelta(visit_date, birth_date)
            age = age_delta.years
            
            demographics['age'] = str(age)
            found_columns['age'] = True
        except:
            # Fallback to year calculation if exact date fails
            if 'PTDOBYY' in row and not pd.isna(row['PTDOBYY']):
                try:
                    birth_year = int(row['PTDOBYY'])
                    visit_year = subject_first_visits[rid].year
                    age = visit_year - birth_year
                    demographics['age'] = str(age)
                    found_columns['age'] = True
                except:
                    pass
    elif 'PTDOBYY' in row and not pd.isna(row['PTDOBYY']) and rid in subject_first_visits:
        # Fallback to year only calculation
        try:
            birth_year = int(row['PTDOBYY'])
            visit_year = subject_first_visits[rid].year
            age = visit_year - birth_year
            demographics['age'] = str(age)
            found_columns['age'] = True
        except:
            pass
    
    # Race - Check for PTRACCAT
    if 'PTRACCAT' in row and not pd.isna(row['PTRACCAT']):
        race_code = row['PTRACCAT']
        demographics['race'] = race_mapping.get(race_code, "Unknown")
        found_columns['race'] = True
    
    # Ethnicity - Check for PTETHCAT
    if 'PTETHCAT' in row and not pd.isna(row['PTETHCAT']):
        ethnicity_code = row['PTETHCAT']
        demographics['ethnicity'] = ethnicity_mapping.get(ethnicity_code, "Unknown")
        found_columns['ethnicity'] = True
    
    # Add to dictionary if we have at least some demographic info
    if demographics:
        subject_demographics[adni_id] = demographics

# Report on columns found
print("\nDemographic fields found:")
for field, found in found_columns.items():
    print(f"  {field}: {'Yes' if found else 'No'}")

# Save as JSON
output_dir = os.path.dirname(output_json_path)
os.makedirs(output_dir, exist_ok=True)

# Create the metadata JSON structure
adni_metadata = {
    "dataset_info": {
        "source": "ADNI",
        "description": "Alzheimer's Disease Neuroimaging Initiative microarray data with demographic information",
        "version": "April 2025",
        "data_type": "microarray",
        "reference_genome": "hg38",
        "gencode_version": "24"
    },
    "subject_demographics": subject_demographics
}

# Save to JSON file
with open(output_json_path, 'w') as f:
    json.dump(adni_metadata, f, indent=2)

print(f"Saved ADNI metadata for {len(subject_demographics)} subjects to {output_json_path}")

# Create a simple report of the demographics distribution
print("\nDemographic Distribution:")
sex_count = {"male": 0, "female": 0, "unknown": 0}
race_count = {}
ethnicity_count = {}
age_ranges = {"<65": 0, "65-75": 0, "75-85": 0, ">85": 0, "unknown": 0}

for subject, demo in subject_demographics.items():
    # Sex
    sex = demo.get('sex', 'unknown')
    sex_count[sex] = sex_count.get(sex, 0) + 1
    
    # Race
    race = demo.get('race', 'Unknown')
    race_count[race] = race_count.get(race, 0) + 1
    
    # Ethnicity
    ethnicity = demo.get('ethnicity', 'Unknown')
    ethnicity_count[ethnicity] = ethnicity_count.get(ethnicity, 0) + 1
    
    # Age
    if 'age' in demo:
        try:
            age = int(float(demo['age']))
            if age < 65:
                age_ranges["<65"] += 1
            elif age < 75:
                age_ranges["65-75"] += 1
            elif age < 85:
                age_ranges["75-85"] += 1
            else:
                age_ranges[">85"] += 1
        except ValueError:
            age_ranges["unknown"] += 1
    else:
        age_ranges["unknown"] += 1

print("\nSex distribution:")
total = sum(sex_count.values())
for sex, count in sex_count.items():
    print(f"  {sex}: {count} ({count/total*100:.1f}%)")

print("\nRace distribution:")
total = sum(race_count.values())
for race, count in sorted(race_count.items(), key=lambda x: x[1], reverse=True):
    print(f"  {race}: {count} ({count/total*100:.1f}%)")

print("\nEthnicity distribution:")
total = sum(ethnicity_count.values())
for ethnicity, count in sorted(ethnicity_count.items(), key=lambda x: x[1], reverse=True):
    print(f"  {ethnicity}: {count} ({count/total*100:.1f}%)")

print("\nAge distribution:")
total = sum(age_ranges.values())
for age_range, count in age_ranges.items():
    print(f"  {age_range}: {count} ({count/total*100:.1f}%)")