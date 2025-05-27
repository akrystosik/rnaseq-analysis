#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def analyze_gene_overlap():
    """Analyze gene overlap between datasets"""
    
    data_dir = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137'
    dataset_files = {
        'ADNI': 'adni_standardized_preprocessed.h5ad',
        'ENCODE': 'encode_standardized_preprocessed.h5ad', 
        'GTEx': 'gtex_standardized_preprocessed.h5ad',
        'MAGE': 'mage_standardized_preprocessed.h5ad'
    }
    
    print('🧬 GENE OVERLAP ANALYSIS')
    print('=' * 50)
    
    gene_sets = {}
    dataset_info = {}
    
    # Load all datasets and collect gene information
    for name, filename in dataset_files.items():
        file_path = f'{data_dir}/{filename}'
        
        try:
            adata = sc.read_h5ad(file_path)
            
            # Get gene information
            gene_ids = set(adata.var_names)
            gene_sets[name] = gene_ids
            
            dataset_info[name] = {
                'total_genes': len(gene_ids),
                'samples': adata.n_obs,
                'sample_gene_ids': list(gene_ids)[:5]
            }
            
            print(f'\n{name} Dataset:')
            print(f'  Total genes: {len(gene_ids):,}')
            print(f'  Sample gene IDs: {list(gene_ids)[:5]}')
            
            # Check gene ID formats
            ensembl_count = sum(1 for g in gene_ids if g.startswith('ENSG'))
            entrez_count = sum(1 for g in gene_ids if g.startswith('ENTREZ:'))
            other_count = len(gene_ids) - ensembl_count - entrez_count
            
            print(f'  Ensembl IDs: {ensembl_count:,} ({ensembl_count/len(gene_ids)*100:.1f}%)')
            print(f'  Entrez IDs: {entrez_count:,}')
            print(f'  Other IDs: {other_count:,}')
            
        except Exception as e:
            print(f'❌ Error loading {name}: {e}')
    
    # Analyze overlaps
    print(f'\n🔍 GENE OVERLAP ANALYSIS:')
    print('=' * 40)
    
    if len(gene_sets) >= 2:
        # Calculate pairwise overlaps
        dataset_names = list(gene_sets.keys())
        
        print('Pairwise overlaps:')
        for i, name1 in enumerate(dataset_names):
            for j, name2 in enumerate(dataset_names):
                if i < j:
                    overlap = gene_sets[name1] & gene_sets[name2]
                    union = gene_sets[name1] | gene_sets[name2]
                    jaccard = len(overlap) / len(union) if len(union) > 0 else 0
                    
                    print(f'  {name1} ∩ {name2}: {len(overlap):,} genes')
                    print(f'    ({len(overlap)/len(gene_sets[name1])*100:.1f}% of {name1}, '
                          f'{len(overlap)/len(gene_sets[name2])*100:.1f}% of {name2})')
                    print(f'    Jaccard similarity: {jaccard:.3f}')
        
        # Calculate total unique genes
        all_genes = set()
        for gene_set in gene_sets.values():
            all_genes.update(gene_set)
        
        print(f'\n📊 SUMMARY STATISTICS:')
        print(f'  Total unique genes across all datasets: {len(all_genes):,}')
        print(f'  Sum of individual dataset genes: {sum(len(gs) for gs in gene_sets.values()):,}')
        print(f'  Redundancy factor: {sum(len(gs) for gs in gene_sets.values()) / len(all_genes):.2f}x')
        
        # Find core genes (present in all datasets)
        core_genes = gene_sets[dataset_names[0]]
        for name in dataset_names[1:]:
            core_genes = core_genes & gene_sets[name]
        
        print(f'  Core genes (in all datasets): {len(core_genes):,}')
        print(f'  Core gene examples: {list(core_genes)[:5] if core_genes else "None"}')
        
        # Find dataset-specific genes
        print(f'\n🎯 DATASET-SPECIFIC GENES:')
        for name in dataset_names:
            other_genes = set()
            for other_name in dataset_names:
                if other_name != name:
                    other_genes.update(gene_sets[other_name])
            
            unique_genes = gene_sets[name] - other_genes
            print(f'  {name} unique genes: {len(unique_genes):,}')
            if unique_genes:
                print(f'    Examples: {list(unique_genes)[:3]}')
        
        # Analyze the 161k total issue
        print(f'\n🚨 VALIDATION REPORT ANALYSIS:')
        print('=' * 35)
        print(f'The validation report shows 161,993 "total genes"')
        print(f'This appears to be the SUM of genes across datasets: {sum(len(gs) for gs in gene_sets.values()):,}')
        print(f'But the actual UNIQUE genes across all datasets: {len(all_genes):,}')
        print(f'')
        print(f'📊 Breakdown:')
        for name, gene_set in gene_sets.items():
            print(f'  {name}: {len(gene_set):,} genes')
        print(f'  ──────────────────')
        print(f'  SUM: {sum(len(gs) for gs in gene_sets.values()):,} genes (what validation reports)')
        print(f'  UNIQUE: {len(all_genes):,} genes (actual unique genes)')
        
        return {
            'individual_counts': {name: len(gs) for name, gs in gene_sets.items()},
            'total_sum': sum(len(gs) for gs in gene_sets.values()),
            'unique_total': len(all_genes),
            'core_genes': len(core_genes),
            'all_genes': all_genes,
            'gene_sets': gene_sets
        }

if __name__ == '__main__':
    results = analyze_gene_overlap()