import sys
sys.path.append('/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts')
from gene_expression_tools import load_expression_data, get_gene_expression, ExpressionDataLoader

# Define the file path for the dataset
combined_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/combined_all_genes_standardized.h5ad"
adni_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/standardized_data/adni_standardized.h5ad"

print(f"Loading files directly...")

# Method 1: Load using specific_files parameter
loader = load_expression_data(specific_files={
    "combined_direct": combined_file,
    "adni": adni_file
})

# Test combined dataset first
print("\nTesting combined dataset:")
expr = get_gene_expression("combined_direct", "ENSG00000000003", data_loader=loader)
if expr:
    print(f"  Gene: ENSG00000000003")
    print(f"  Mean expression: {expr.get('mean_expression')}")
    print(f"  Sample count: {expr.get('sample_count')}")

# Test adni dataset
print("\nTesting adni dataset:")
expr = get_gene_expression("adni", "ENSG00000000003", data_loader=loader)
if expr:
    print(f"  Gene: ENSG00000000003")
    print(f"  Mean expression: {expr.get('mean_expression')}")
    print(f"  Sample count: {expr.get('sample_count')}")
else:
    print("  No expression data found for ENSG00000000003 in adni dataset")

# Method 2: Create a loader and use load_file_directly
print("\nTesting with ExpressionDataLoader.load_file_directly:")
loader2 = ExpressionDataLoader()
loader2.load_file_directly("adni_test", adni_file)

expr = get_gene_expression("adni_test", "ENSG00000000003", data_loader=loader2)
if expr:
    print(f"  Gene: ENSG00000000003")
    print(f"  Mean expression: {expr.get('mean_expression')}")
    print(f"  Sample count: {expr.get('sample_count')}")
else:
    print("  No expression data found for ENSG00000000003 in adni_test dataset")