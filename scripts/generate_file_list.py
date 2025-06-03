import os
import sys

def generate_markdown(directory, output_file="file_list.md"):
    """
    Generate a Markdown file listing all .py and .sh files in the given directory,
    including their contents.

    Args:
        directory (str): The directory to scan for files.
        output_file (str): The path to the output Markdown file.
    """
    # Initialize a list to store file paths
    file_paths = []

    # Walk through the directory
    for root, _, files in os.walk(directory):
        for file in files:
            # Check if the file has a .py or .sh extension and exclude .md files and this script
            if (file.endswith(".py") or file.endswith(".sh")) and file != os.path.basename(__file__):
                # Append the relative file path
                file_paths.append(os.path.join(root, file))

    # Write the file paths and their contents to the Markdown file
    with open(output_file, "w") as md_file:
        md_file.write("# File List with Contents\n\n")
        md_file.write("The following files were found, along with their contents:\n\n")
        for path in file_paths:
            md_file.write(f"## `{path}`\n\n")
            md_file.write("````\n")  # Use quadruple backticks for Markdown code blocks
            try:
                with open(path, "r") as file:
                    md_file.write(file.read())
            except Exception as e:
                md_file.write(f"Error reading file: {e}")
            md_file.write("\n````\n\n")

if __name__ == "__main__":
    # Check if a directory argument is provided
    if len(sys.argv) > 1:
        directory_to_scan = sys.argv[1]
    else:
        directory_to_scan = os.getcwd()  # Default to the current working directory

    # Define the output Markdown file
    output_markdown_file = "file_list.md"

    # Generate the Markdown file
    generate_markdown(directory_to_scan, output_markdown_file)
    print(f"Markdown file generated: {output_markdown_file}")