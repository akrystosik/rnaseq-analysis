#!/usr/bin/env python3
"""
ADNI Metadata Integration Script

Updates the standardize_datasets.py script to incorporate ADNI demographic metadata
during the standardization process.
"""

import re
import os

# Path to the script to modify
script_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py'

# Create a backup of the original script
backup_path = f"{script_path}.bak_{os.path.basename(os.path.dirname(script_path))}"
os.system(f"cp {script_path} {backup_path}")
print(f"Created backup at {backup_path}")

# Load the contents of the script
with open(script_path, 'r') as f:
    content = f.read()

# Find the process_adni_data function
adni_function_pattern = re.compile(r'def process_adni_data\([^)]*\):.*?return adata', re.DOTALL)
adni_function_match = adni_function_pattern.search(content)

if not adni_function_match:
    print("Could not find the process_adni_data function")
    exit(1)

original_function = adni_function_match.group(0)

# Create the modified function with ADNI metadata integration
modified_function = original_function.replace(
    "# Create metadata dictionary",
    """# Load ADNI subject demographics if available
    adni_metadata_file = os.path.join(metadata_dir, 'adni_metadata.json') if metadata_dir else None
    subject_demographics = {}
    
    if adni_metadata_file and os.path.exists(adni_metadata_file):
        try:
            with open(adni_metadata_file, 'r') as f:
                adni_metadata = json.load(f)
                if 'subject_demographics' in adni_metadata:
                    subject_demographics = adni_metadata['subject_demographics']
                    logger.info(f"Loaded demographics for {len(subject_demographics)} ADNI subjects")
        except Exception as e:
            logger.warning(f"Error loading ADNI metadata: {e}")
    
    # Create metadata dictionary"""
)

# Add code to incorporate subject demographics
modified_function = modified_function.replace(
    "metadata = {",
    """# Extract subject ID from file name
    subject_id_match = re.search(r'(\\d+)_S_(\\d+)', subject_id)
    if subject_id_match:
        adni_rid = f"ADNI_{subject_id_match.group(2)}"
        # Look up demographics if available
        if adni_rid in subject_demographics:
            demo = subject_demographics[adni_rid]
            # Add demographics to metadata
            if 'sex' in demo:
                metadata['sex'] = demo['sex']
            if 'age' in demo:
                metadata['age'] = demo['age']
            if 'race' in demo:
                metadata['race'] = demo['race']
            if 'ethnicity' in demo:
                metadata['ethnicity'] = demo['ethnicity']
            logger.debug(f"Added demographics for {adni_rid}")
    
    metadata = {"""
)

# Replace the original function with the modified version
updated_content = content.replace(original_function, modified_function)

# Make sure the imports exist
if "import json" not in updated_content:
    updated_content = updated_content.replace(
        "import logging",
        "import logging\nimport json"
    )

if "import re" not in updated_content:
    updated_content = updated_content.replace(
        "import logging",
        "import logging\nimport re"
    )

# Write the updated content back to the file
with open(script_path, 'w') as f:
    f.write(updated_content)

print("Successfully added ADNI metadata integration to standardize_datasets.py")
print("You can now run the pipeline with enhanced ADNI demographic information")