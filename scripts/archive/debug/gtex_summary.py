# GTEx Single Cell RNA-seq Validation Script
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

print("Starting GTEx scRNA-seq validation...")

# File path
file_path = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/raw_data/snRNAseq_atlas/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad'

try:
    # Load the h5ad file
    print(f"Loading {file_path}...")
    adata = sc.read_h5ad(file_path)
    
    # Basic dataset information
    print("\n----- BASIC DATASET INFORMATION -----")
    print(f"Dataset dimensions: {adata.shape[0]} cells × {adata.shape[1]} genes")
    print(f"Available observation (cell) metadata: {list(adata.obs.columns)}")
    
    # Tissue analysis
    print("\n----- TISSUE ANALYSIS -----")
    if 'tissue' in adata.obs.columns:
        tissues = adata.obs['tissue'].unique()
        print(f"Number of tissues: {len(tissues)}")
        print("Tissues in dataset:")
        
        # Get cell count per tissue
        tissue_counts = adata.obs['tissue'].value_counts().sort_values(ascending=False)
        for tissue, count in tissue_counts.items():
            print(f"- {tissue}: {count:,} cells")
    else:
        # Try to find alternate tissue column names
        potential_tissue_cols = [col for col in adata.obs.columns if 'tissue' in col.lower()]
        if potential_tissue_cols:
            print(f"Tissue column not found, but found potential alternatives: {potential_tissue_cols}")
        else:
            print("No tissue column found in dataset")
    
    # Donor analysis
    print("\n----- DONOR ANALYSIS -----")
    donor_cols = [col for col in adata.obs.columns if 'donor' in col.lower()]
    if donor_cols:
        donor_col = donor_cols[0]
        donors = adata.obs[donor_col].unique()
        print(f"Number of donors: {len(donors)}")
        
        # Count cells per donor
        donor_counts = adata.obs[donor_col].value_counts().sort_values(ascending=False)
        print("Top 10 donors by cell count:")
        for i, (donor, count) in enumerate(donor_counts.head(10).items()):
            print(f"- {donor}: {count:,} cells")
            
        # Check donors per tissue
        if 'tissue' in adata.obs.columns:
            print("\nDonors per tissue:")
            for tissue in sorted(tissues):
                tissue_donors = adata.obs[adata.obs['tissue'] == tissue][donor_col].unique()
                print(f"- {tissue}: {len(tissue_donors)} donors")
    else:
        print("No donor column found in dataset")
    
    # Cell type analysis
    print("\n----- CELL TYPE ANALYSIS -----")
    celltype_cols = [col for col in adata.obs.columns if 'cell' in col.lower() and 'type' in col.lower()]
    if celltype_cols:
        celltype_col = celltype_cols[0]
        celltypes = adata.obs[celltype_col].unique()
        print(f"Number of cell types: {len(celltypes)}")
        
        # Count cells per cell type
        celltype_counts = adata.obs[celltype_col].value_counts().sort_values(ascending=False)
        print("Top 15 cell types by cell count:")
        for celltype, count in celltype_counts.head(15).items():
            print(f"- {celltype}: {count:,} cells")
            
        # Count cell types per tissue
        if 'tissue' in adata.obs.columns:
            print("\nNumber of cell types per tissue:")
            for tissue in sorted(tissues):
                tissue_celltypes = adata.obs[adata.obs['tissue'] == tissue][celltype_col].unique()
                print(f"- {tissue}: {len(tissue_celltypes)} cell types")
    else:
        print("No cell type column found in dataset")
    
    # Check for breast tissue specifically
    print("\n----- BREAST TISSUE CHECK -----")
    if 'tissue' in adata.obs.columns:
        breast_cells = adata.obs[adata.obs['tissue'].str.lower().str.contains('breast')].shape[0]
        if breast_cells > 0:
            print(f"FOUND breast tissue: {breast_cells:,} cells")
        else:
            print("NO breast tissue found in the dataset")
    
    # Check for skeletal muscle specifically
    print("\n----- SKELETAL MUSCLE CHECK -----")
    if 'tissue' in adata.obs.columns:
        muscle_cells = adata.obs[adata.obs['tissue'].str.lower().str.contains('muscle')].shape[0]
        if muscle_cells > 0:
            muscle_tissues = adata.obs[adata.obs['tissue'].str.lower().str.contains('muscle')]['tissue'].unique()
            print(f"FOUND muscle tissue: {muscle_tissues}")
            for mt in muscle_tissues:
                count = adata.obs[adata.obs['tissue'] == mt].shape[0]
                print(f"- {mt}: {count:,} cells")
        else:
            print("NO muscle tissue found in the dataset")
    
    # Total cell count verification
    print("\n----- TOTAL CELL COUNT VERIFICATION -----")
    total_cells = adata.shape[0]
    print(f"Total cells in dataset: {total_cells:,}")
    print(f"Matches document claim of 209,126 cells: {'Yes' if total_cells == 209126 else 'No'}")
    
    print("\nValidation complete!")

except Exception as e:
    print(f"Error while processing file: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)