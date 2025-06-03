# generate_metadata_list.py
import os
import sys

# --- Define the EXACT metadata files used by the pipeline ---
# Adjust these paths if they differ in your setup
METADATA_FILES_TO_INCLUDE = [
    # Dataset-specific JSONs
#    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/adni_metadata.json",
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/encode_metadata.json",
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gtex_metadata.json",
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/mage_metadata.json",

    # Ontology Mapping JSONs
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/tissue_to_uberon.json",
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/assay_to_efo.json",
    "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/age_to_hsapdv.json",

    # # Gene ID / Mapping Files (CSVs)
    # "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/entrez_to_ensembl_mapping.csv",
    # "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_id_reference_mapping.csv",
    # "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gene_mapping/encode_id_to_ensembl_mapping.csv",
    # "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/gencode_v24_complete_mapping.csv", # Used by rnaseq_utils

    # Add any other specific files if needed, e.g., GTEx attributes if directly used
    # "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
    # "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
]

def generate_markdown(file_list, output_file="metadata_file_list.md"):
    """
    Generate a Markdown file listing specific files and their contents.

    Args:
        file_list (list): A list of absolute file paths to include.
        output_file (str): The path to the output Markdown file.
    """
    # Write the file paths and their contents to the Markdown file
    with open(output_file, "w") as md_file:
        md_file.write("# Metadata File Contents\n\n")
        md_file.write("The following metadata files were used or referenced by the pipeline, along with their contents:\n\n")

        for path in file_list:
            # Determine the file extension for syntax highlighting hint
            _, ext = os.path.splitext(path)
            lang_hint = ext.lstrip('.').lower()
            if lang_hint == 'py': lang_hint = 'python'
            if lang_hint == 'sh': lang_hint = 'bash'
            if lang_hint == 'txt': lang_hint = '' # No hint for txt

            md_file.write(f"## `{path}`\n\n")
            if os.path.exists(path):
                md_file.write(f"```{lang_hint}\n")
                try:
                    # For large CSVs, maybe just show head? - Let's include full for now.
                    # Limit read size?
                    with open(path, "r", encoding='utf-8', errors='ignore') as file:
                        # Decide whether to limit large files - e.g., mapping files
                        # file_size = os.path.getsize(path)
                        # if file_size > 5 * 1024 * 1024: # Limit files > 5MB
                        #      md_file.write(f"[File larger than 5MB - Showing first 100 lines]\n\n")
                        #      for i, line in enumerate(file):
                        #          if i >= 100:
                        #              break
                        #          md_file.write(line)
                        # else:
                        md_file.write(file.read())

                except Exception as e:
                    md_file.write(f"Error reading file: {e}")
                md_file.write("\n```\n\n")
            else:
                md_file.write("*File not found.*\n\n")

if __name__ == "__main__":
    # Define the output Markdown file
    output_markdown_file = "metadata_file_list.md"

    # Generate the Markdown file
    generate_markdown(METADATA_FILES_TO_INCLUDE, output_markdown_file)
    print(f"Metadata Markdown file generated: {output_markdown_file}")