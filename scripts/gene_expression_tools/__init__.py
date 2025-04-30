"""
Gene Expression Tools

Package for high-performance gene expression data analysis.
"""

from .optimized_gene_expression import (
    load_expression_data,
    get_gene_expression,
    get_gene_expression_matrix,
    get_available_tissues,
    get_available_donors,
    get_all_tissues_gene_expression,
    get_tissue_gene_expression,
    get_donor_gene_expression,
    ExpressionDataLoader
)

__all__ = [
    'load_expression_data',
    'get_gene_expression',
    'get_gene_expression_matrix',
    'get_available_tissues',
    'get_available_donors',
    'get_all_tissues_gene_expression',
    'get_tissue_gene_expression',
    'get_donor_gene_expression',
    'ExpressionDataLoader'
]
