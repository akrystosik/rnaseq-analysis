#!/usr/bin/env python3
"""
Script to locate and fix indentation issues in standardize_datasets.py
"""

import sys
import re

# Path to the input file
input_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets.py"
# Path to the output file
output_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/standardize_datasets_fixed.py"

# Read the file
with open(input_file, 'r') as f:
    lines = f.readlines()

# Print the problematic area around line 1885
problem_area_start = max(1875, 0)
problem_area_end = min(1895, len(lines))
print("Problematic code area:")
for i in range(problem_area_start, problem_area_end):
    print(f"{i+1}: {lines[i]}", end='')
    
# Look for indentation issues around line 1885
# Specifically, looking for lines with 'metadata = {'
metadata_pattern = re.compile(r'^\s+metadata = \{')
fixed_lines = []
fixed_something = False

for i, line in enumerate(lines):
    if i+1 >= 1880 and i+1 <= 1890 and metadata_pattern.match(line):
        # This might be the problematic line
        # Check previous line for indentation level
        prev_line = lines[i-1]
        prev_indent = len(prev_line) - len(prev_line.lstrip())
        curr_indent = len(line) - len(line.lstrip())
        
        if curr_indent > prev_indent:
            # This is likely the problem - fix the indentation
            fixed_line = ' ' * prev_indent + line.lstrip()
            fixed_lines.append(fixed_line)
            fixed_something = True
            print(f"\nFound and fixed problematic line {i+1}:")
            print(f"Original: {line}")
            print(f"Fixed: {fixed_line}")
        else:
            fixed_lines.append(line)
    else:
        fixed_lines.append(line)

if fixed_something:
    # Write the fixed content to a new file
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    print(f"\nFixed file saved to {output_file}")
    print("You can replace the original file with the fixed version using:")
    print(f"cp {output_file} {input_file}")
else:
    print("\nCouldn't automatically identify the indentation issue.")
    print("Manual inspection is needed. Try checking if there are tabs mixed with spaces.")