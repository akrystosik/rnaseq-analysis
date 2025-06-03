#!/usr/bin/env python3
import json
import os
import glob

# Load the original metadata
metadata_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/metadata/entex_metadata.json"
with open(metadata_path, 'r') as f:
    metadata = json.load(f)

# Find all TSV files in the ENCODE directory
base_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/encode/"
tsv_files = glob.glob(f"{base_dir}**/*.tsv", recursive=True)

# Create a mapping from file ID to file path
file_id_to_path = {}
for file_path in tsv_files:
    file_name = os.path.basename(file_path)
    file_id = file_name.split('.')[0]  # e.g., ENCFF008RKC
    file_id_to_path[file_id] = file_path

# Update file_path values in metadata
updated_count = 0
not_found_count = 0

# First, try the entex_metadata array if it exists
if "entex_metadata" in metadata:
    for sample_info in metadata["entex_metadata"]:
        sample_id = sample_info.get("sample_id")
        if sample_id in file_id_to_path:
            sample_info["file_path"] = file_id_to_path[sample_id]
            updated_count += 1
        else:
            not_found_count += 1

# Then update the file-specific entries
for file_id in list(metadata.keys()):
    if file_id.startswith("ENCFF"):
        if file_id in file_id_to_path:
            metadata[file_id]["file_path"] = file_id_to_path[file_id]
            updated_count += 1
        else:
            not_found_count += 1

# Save updated metadata
backup_path = metadata_path + ".bak"
os.rename(metadata_path, backup_path)
with open(metadata_path, 'w') as f:
    json.dump(metadata, f, indent=2)

print(f"Updated {updated_count} file paths, couldn't find {not_found_count} files.")
print(f"Original metadata backed up to {backup_path}")