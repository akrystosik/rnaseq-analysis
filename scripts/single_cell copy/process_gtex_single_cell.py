#!/usr/bin/env python3
"""
Process GTEx Single-Cell RNA-seq Data:
1. Load h5ad file.
2. Standardize cell-level metadata using rnaseq_utils.
3. Perform pseudobulking by averaging gene expression per donor per cell type.
4. Create a new AnnData object for pseudobulk data.
5. Calculate gene expression variance across pseudobulk samples.
6. Save the pseudobulk AnnData.
"""

import os
import argparse
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
import scipy.sparse as sp
import sys 

# Import shared utilities
try:
    from rnaseq_utils import (
        standardize_metadata,
        load_mappings,
        load_dataset_specific_metadata,
        apply_dataset_specific_metadata,
        ensure_serializable,
        standardize_ensembl_id
    )
    from standardize_datasets import save_anndata as robust_save_anndata
    SAVE_METHOD = "robust"
except ImportError as e:
    logging.error(f"Could not import from rnaseq_utils or standardize_datasets: {e}")
    logging.error("Ensure 'scripts/pipeline' directory is in PYTHONPATH or accessible.")
    logging.warning("Using basic AnnData write as a fallback for saving.")
    
    def robust_save_anndata(adata, file_path):
        try:
            adata.write_h5ad(file_path)
            logging.info(f"Fallback basic AnnData write successful to {file_path}")
            return True
        except Exception as write_err:
            logging.error(f"Fallback basic AnnData write FAILED for {file_path}: {write_err}")
            return False
    SAVE_METHOD = "basic_fallback"
    
    # Minimal fallbacks for other functions if needed early, though current structure avoids this
    def load_mappings(metadata_json_dir=None): return {}
    def standardize_metadata(obs_df, dataset_name, mappings, **kwargs): return obs_df # Passthrough
    def load_dataset_specific_metadata(metadata_dir, dataset_name): return None
    def apply_dataset_specific_metadata(adata, metadata): return adata

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gtex_sc_processor')


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process GTEx single-cell RNA-seq data for pseudobulking.')
    parser.add_argument('--input-h5ad', type=str, required=True,
                        help='Path to the input GTEx single-cell h5ad file.')
    parser.add_argument('--output-h5ad', type=str, required=True,
                        help='Path to save the output pseudobulk h5ad file.')
    parser.add_argument('--metadata-json-dir', type=str, required=True,
                        help='Directory containing metadata JSON files.')
    parser.add_argument('--gene-reference-mapping-csv', type=str, required=True,
                        help='Path to the gene_id_mapping_reference.csv file.') # <<<< NEW ARGUMENT
    return parser.parse_args()


def load_gene_reference_mapping(gene_ref_csv_path: str) -> pd.DataFrame:
    logger.info(f"Loading gene ID reference mapping from: {gene_ref_csv_path}")
    try:
        ref_df = pd.read_csv(gene_ref_csv_path, low_memory=False)
        # Ensure key columns are string type for mapping
        for col in ['gene_id', 'original_gene_id', 'gene_name', 'numeric_id']:
            if col in ref_df.columns:
                ref_df[col] = ref_df[col].astype(str).fillna('')
        logger.info(f"Loaded gene reference mapping with {len(ref_df)} entries.")
        return ref_df
    except Exception as e:
        logger.error(f"Failed to load gene reference mapping from {gene_ref_csv_path}: {e}")
        raise

def standardize_adata_gene_ids(adata_sc: ad.AnnData, gene_ref_df: pd.DataFrame) -> ad.AnnData:
    logger.info("Standardizing gene IDs in input AnnData...")
    original_var = adata_sc.var.copy()
    original_var_names = list(adata_sc.var_names) # Keep a copy of original var_names

    # Create mappings from reference for quick lookup
    # We want to map original IDs (symbols, ENSG with version, Entrez) to a clean ENSG (no version)
    
    # 1. Map by 'gene_id' in ref_df (this is already supposed to be the clean ENSG)
    #    to itself, useful for validating if input is already clean ENSG.
    map_ensg_to_ensg = dict(zip(gene_ref_df['gene_id'], gene_ref_df['gene_id']))
    
    # 2. Map by 'original_gene_id' in ref_df (e.g., ENSG.version) to clean 'gene_id'
    map_orig_ensg_to_ensg = {}
    if 'original_gene_id' in gene_ref_df.columns: # original_gene_id from GENCODE GTF (ENSG.version)
        for _, row in gene_ref_df.iterrows():
            if row['original_gene_id'] and row['gene_id'] and row['original_gene_id'] != row['gene_id']:
                 # original_gene_id might have multiple entries if it's a general symbol, so be careful
                 # Here, we assume original_gene_id from GTF parsing (parse_gencode_gtf) is quite specific (ENSG.version)
                 map_orig_ensg_to_ensg[row['original_gene_id']] = row['gene_id']

    # 3. Map by 'gene_name' (symbol) in ref_df to clean 'gene_id'
    # Handle with care: a symbol can map to multiple ENSG IDs. Prioritize protein_coding.
    map_symbol_to_ensg = {}
    if 'gene_name' in gene_ref_df.columns and 'gene_type' in gene_ref_df.columns:
        # Prioritize protein_coding, then others
        for gene_type_priority in ['protein_coding', 'lncRNA', 'miRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'other']:
            subset_df = gene_ref_df[gene_ref_df['gene_type'] == gene_type_priority]
            for _, row in subset_df.iterrows():
                if row['gene_name'] and row['gene_id'] and row['gene_name'] not in map_symbol_to_ensg: # Add only if not already mapped by higher priority
                    map_symbol_to_ensg[row['gene_name']] = row['gene_id']
        # Add remaining symbols that weren't in priority types
        for _, row in gene_ref_df.iterrows():
             if row['gene_name'] and row['gene_id'] and row['gene_name'] not in map_symbol_to_ensg:
                  map_symbol_to_ensg[row['gene_name']] = row['gene_id']

    # 4. Map by 'numeric_id' (Entrez) in ref_df to clean 'gene_id'
    map_entrez_to_ensg = {}
    if 'numeric_id' in gene_ref_df.columns:
         for _, row in gene_ref_df.iterrows():
            if row['numeric_id'] and row['gene_id'] and not row['numeric_id'].upper().startswith("ENTREZ:"):
                map_entrez_to_ensg[f"ENTREZ:{row['numeric_id']}"] = row['gene_id'] # Add prefix if not present
                map_entrez_to_ensg[row['numeric_id']] = row['gene_id'] # Also map raw numeric string
            elif row['numeric_id'] and row['gene_id']: # Already has ENTREZ: prefix
                map_entrez_to_ensg[row['numeric_id']] = row['gene_id']


    mapped_ids = {}
    unmapped_count = 0
    for original_id in original_var_names:
        std_id = str(original_id) # Start with original
        target_ensg = None
        
        # Attempt mapping (order of preference)
        # a. Is it already a clean ENSG?
        if std_id.startswith("ENSG") and "." not in std_id and std_id in map_ensg_to_ensg:
            target_ensg = std_id
        # b. Try symbol map
        elif std_id in map_symbol_to_ensg:
            target_ensg = map_symbol_to_ensg[std_id]
        # c. Try original Ensembl (with version) map
        elif std_id in map_orig_ensg_to_ensg:
            target_ensg = map_orig_ensg_to_ensg[std_id]
        # d. Try stripping version manually and check if that ENSG base is in ref
        elif std_id.startswith("ENSG") and "." in std_id:
            base_ensg = standardize_ensembl_id(std_id) # Utility from rnaseq_utils
            if base_ensg in map_ensg_to_ensg:
                target_ensg = base_ensg
        # e. Try Entrez map (assuming original_id could be numeric string or ENTREZ:numeric)
        elif std_id in map_entrez_to_ensg:
             target_ensg = map_entrez_to_ensg[std_id]
        elif f"ENTREZ:{std_id}" in map_entrez_to_ensg: # If original_id is just numeric string
             target_ensg = map_entrez_to_ensg[f"ENTREZ:{std_id}"]

        if target_ensg and target_ensg.startswith("ENSG"):
            mapped_ids[original_id] = target_ensg
        else:
            # If still no mapping, keep original for now, or decide on placeholder/drop
            # For GTEx SC, the var_names are often symbols or ENSG with version.
            # If it's already ENSG with version, standardize_ensembl_id can clean it.
            cleaned_original = standardize_ensembl_id(std_id)
            if cleaned_original.startswith("ENSG"):
                mapped_ids[original_id] = cleaned_original # Use cleaned version if it's ENSG
                if not target_ensg: logger.debug(f"Gene ID '{original_id}' cleaned to '{cleaned_original}' (not found in ref map directly).")
            else:
                unmapped_count += 1
                logger.warning(f"Gene ID '{original_id}' could not be mapped to a standard Ensembl ID and is not ENSG. Kept as original for now.")
                mapped_ids[original_id] = original_id # Keep original if unmappable and not ENSG-like

    logger.info(f"Gene ID mapping: {len(mapped_ids) - unmapped_count} mapped, {unmapped_count} unmapped or kept as original non-ENSG.")

    # --- Aggregate expression for genes mapping to the same standard ID ---
    # Create a new var_names list and an aggregation mapping
    # new_var_meta_list will store rows for the new var_df
    new_var_meta_list = []
    # temp_X_data will store columns of the new expression matrix
    temp_X_data_cols = []
    
    # Group original IDs by their new standard target ID
    target_to_originals_map = {}
    for original_id, target_id in mapped_ids.items():
        target_to_originals_map.setdefault(target_id, []).append(original_id)

    new_standard_var_names = sorted(list(target_to_originals_map.keys()))
    
    # Create new expression matrix by summing counts for genes that map to the same standard ID
    # Assuming adata_sc.X is cells x genes
    if sp.issparse(adata_sc.X):
        new_X = sp.lil_matrix((adata_sc.n_obs, len(new_standard_var_names)), dtype=adata_sc.X.dtype)
    else:
        new_X = np.zeros((adata_sc.n_obs, len(new_standard_var_names)), dtype=adata_sc.X.dtype)

    for i, target_id in enumerate(new_standard_var_names):
        original_ids_for_target = target_to_originals_map[target_id]
        original_indices = [original_var_names.index(orig_id) for orig_id in original_ids_for_target if orig_id in original_var_names]
        
        if original_indices:
            # Sum expression across the original gene columns for this target_id
            current_sum = adata_sc.X[:, original_indices[0]].copy() # Start with the first one
            for k in range(1, len(original_indices)):
                current_sum += adata_sc.X[:, original_indices[k]]
            new_X[:, i] = current_sum
            
            # For var metadata, take from the first original ID, or merge if needed
            # Here, we just take from the first original ID that mapped to this target
            first_original_id = original_ids_for_target[0]
            var_row_data = original_var.loc[first_original_id].copy()
            var_row_data['original_ids_mapped_here'] = ";".join(original_ids_for_target) # Store which original IDs were summed
            var_row_data.name = target_id # Set the name for the new Series/row
            new_var_meta_list.append(var_row_data)
        else:
            logger.warning(f"No original indices found for target_id {target_id}. This shouldn't happen.")


    new_var_df = pd.DataFrame(new_var_meta_list)
    if not new_var_df.empty:
         new_var_df.index = new_standard_var_names # Index is now the standard gene ID
    else: # Handle case where no genes were processed successfully
        logger.error("new_var_df is empty after gene ID standardization. Aborting further processing.")
        # Create an empty DataFrame with the expected index to avoid downstream errors if processing must continue partially
        new_var_df = pd.DataFrame(index=pd.Index(new_standard_var_names, name='feature_id'))


    new_var_df['gene_id'] = new_var_df.index.astype(str) # Ensure 'gene_id' column matches new index

    # Create new AnnData with standardized genes
    adata_standardized_genes = ad.AnnData(X=new_X, obs=adata_sc.obs.copy(), var=new_var_df)
    # Copy .uns if necessary (though pseudobulking happens next, then .uns is set)
    adata_standardized_genes.uns = adata_sc.uns.copy() 
    
    logger.info(f"Gene ID standardization complete. New .var shape: {adata_standardized_genes.var.shape}, .X shape: {adata_standardized_genes.X.shape}")
    return adata_standardized_genes


def prepare_metadata_for_standardization(obs_df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Preparing raw metadata for standardization...")
    obs_prep = obs_df.copy()

    # --- ADJUSTED RENAME MAP based on actual logs ---
    # Prioritize specific source columns for target standard names
    source_to_target_map = {
        'subject_id': ['individual', 'Participant ID'], # Try 'individual' first, then 'Participant ID'
        'original_cell_id': ['Sample ID', 'barcode'],   # Try 'Sample ID' first, then 'barcode' (if 'barcode' is a column)
        'tissue': ['Tissue', 'tissue'],                 # Try 'Tissue' (uppercase), then 'tissue' (lowercase)
        'sex_original': ['Sex'],
        'age_original': ['Age_bin']
        # cell_type is handled separately
    }

    # Process renames carefully
    for target_col, source_candidates in source_to_target_map.items():
        found_source = False
        for source_col in source_candidates:
            if source_col in obs_prep.columns:
                if target_col not in obs_prep.columns:
                    obs_prep.rename(columns={source_col: target_col}, inplace=True)
                    logger.info(f"Renamed column '{source_col}' to '{target_col}'")
                elif source_col != target_col: # Target exists, source is different, overwrite
                    logger.warning(f"Target column '{target_col}' exists. Overwriting with data from '{source_col}'.")
                    obs_prep[target_col] = obs_prep[source_col]
                    if source_col not in source_to_target_map.keys(): # Avoid deleting if source_col is also a target_col for another field
                        obs_prep.drop(columns=[source_col], inplace=True)
                # If source_col == target_col and target_col already exists, do nothing.
                found_source = True
                break # Found a source for this target, move to next target
        if not found_source and target_col not in obs_prep.columns:
            logger.warning(f"Could not find any source column for target '{target_col}' from candidates: {source_candidates}. '{target_col}' will be missing.")


    # Handle cell type: Prioritize specific columns
    cell_type_candidates_in_order = [
        'Granular cell type', # From original log
        'Granular_Allele_Specific_Cell_Types_Predicted_Cell_Type', # From previous script version
        'Broad cell type',    # From original log
        'Cell types level 3', # From original log
        'Cell types level 2'  # From original log
    ]
    chosen_cell_type_col_source = None
    if 'cell_type' not in obs_prep.columns: # Only if 'cell_type' isn't already perfectly named
        for candidate in cell_type_candidates_in_order:
            if candidate in obs_prep.columns:
                obs_prep.rename(columns={candidate: 'cell_type'}, inplace=True)
                logger.info(f"Using '{candidate}' as primary 'cell_type'.")
                chosen_cell_type_col_source = candidate
                # Rename other candidates to avoid conflicts if they were not chosen
                for other_candidate in cell_type_candidates_in_order:
                    if other_candidate != chosen_cell_type_col_source and other_candidate in obs_prep.columns:
                        safe_alt_name = f"{other_candidate.lower().replace(' ', '_').replace('-', '_')}_alt"
                        obs_prep.rename(columns={other_candidate: safe_alt_name}, inplace=True)
                        logger.info(f"Renamed alternative cell type column '{other_candidate}' to '{safe_alt_name}'.")
                break
    
    if 'cell_type' not in obs_prep.columns:
        logger.error("CRITICAL: No suitable cell type column found. Cannot proceed.")
        raise ValueError("Missing 'cell_type' column.")

    # Handle 'sample_id' for each cell (this should be unique cell identifier)
    if obs_df.index.name is not None and obs_df.index.name.lower() == 'barcode':
        obs_prep['sample_id'] = obs_df.index.astype(str)
        if 'barcode_original' not in obs_prep.columns:
             obs_prep['barcode_original'] = obs_df.index.astype(str)
    elif 'original_cell_id' in obs_prep.columns: # This comes from 'Sample ID' or 'barcode' (if column) via source_to_target_map
         obs_prep['sample_id'] = obs_prep['original_cell_id'].astype(str)
    else: # Fallback if no clear unique cell ID column was mapped
        logger.warning("Suitable unique cell identifier ('barcode' index, 'Sample ID' column) not found/mapped to 'original_cell_id'. Creating new unique sample_ids for cells.")
        obs_prep['sample_id'] = [f"cell_{i}" for i in range(len(obs_prep))]
    
    # Standardize 'sex'
    if 'sex_original' in obs_prep.columns:
        sex_map_str = {'male': 'male', 'female': 'female'} # GTEx SC uses 'Male'/'Female'
        obs_prep['sex'] = obs_prep['sex_original'].astype(str).str.lower().map(sex_map_str).fillna('unknown')
        logger.info("Mapped string 'sex_original' to 'sex' column.")
    elif 'sex' not in obs_prep.columns: # If 'sex_original' wasn't found and 'sex' doesn't exist
        logger.warning("Column for 'sex' or its source ('Sex') not found. 'sex' will be 'unknown'.")
        obs_prep['sex'] = 'unknown'

    # Standardize 'age'
    if 'age_original' in obs_prep.columns: # From 'Age_bin'
        obs_prep['age'] = obs_prep['age_original'].astype(str).fillna('')
        logger.info("Copied 'age_original' (from 'Age_bin') to 'age' column as string.")
    elif 'age' not in obs_prep.columns:
        logger.warning("Column for 'age' or its source ('Age_bin') not found. 'age' will be empty.")
        obs_prep['age'] = ''

    # Ensure default columns needed by standardize_metadata
    if 'dataset' not in obs_prep.columns: obs_prep['dataset'] = 'GTEx-snRNAseq'
    if 'data_type' not in obs_prep.columns: obs_prep['data_type'] = 'snRNA-seq'
    if 'expression_unit' not in obs_prep.columns: obs_prep['expression_unit'] = 'UMI counts'

    # Final check for critical columns for pseudobulking
    if 'subject_id' not in obs_prep.columns:
        logger.error("CRITICAL: 'subject_id' (donor) column still not found. Cannot proceed.")
        raise ValueError("Missing 'subject_id' column for pseudobulking.")
    
    logger.info("Metadata preparation complete.")
    # Before returning, ensure 'tissue' is a Series
    if 'tissue' in obs_prep.columns and isinstance(obs_prep['tissue'], pd.DataFrame):
        logger.error(f"CRITICAL FAILURE in prepare_metadata: obs_prep['tissue'] is a DataFrame before returning. Columns: {obs_prep['tissue'].columns.tolist()}")
        # Attempt one last fix or raise error
        if len(obs_prep['tissue'].columns) > 0:
            obs_prep['tissue'] = obs_prep['tissue'].iloc[:,0]
        else:
            raise ValueError("obs_prep['tissue'] became an empty DataFrame.")
            
    return obs_prep

def pseudobulk_anndata(adata_sc: ad.AnnData, group_by_cols: list) -> ad.AnnData:
    logger.info(f"Starting pseudobulking by: {group_by_cols}")

    if not all(col in adata_sc.obs.columns for col in group_by_cols):
        missing = [col for col in group_by_cols if col not in adata_sc.obs.columns]
        logger.error(f"Missing columns for pseudobulking: {missing}")
        raise ValueError(f"Missing columns in .obs for grouping: {missing}")

    if adata_sc.X.dtype == 'object':
        logger.warning("Expression data (.X) is of object type. Attempting conversion to numeric.")
        try:
            if sp.issparse(adata_sc.X):
                adata_sc.X = adata_sc.X.astype(np.float32)
            else:
                adata_sc.X = adata_sc.X.astype(np.float32)
        except Exception as e:
            logger.error(f"Failed to convert .X to numeric: {e}. Pseudobulking might fail.")
            raise

    try:
        sc_version_parts = sc.__version__.split('.')
        sc_major = int(sc_version_parts[0])
        sc_minor = int(sc_version_parts[1])
        if sc_major > 1 or (sc_major == 1 and sc_minor >= 7): # Scanpy 1.7+
            logger.info("Using scanpy.experimental.pp.aggregate_profiles for pseudobulking.")
            # Ensure grouping columns are categorical for aggregate_profiles
            for col in group_by_cols:
                if not pd.api.types.is_categorical_dtype(adata_sc.obs[col]):
                    logger.info(f"Converting column '{col}' to categorical for aggregate_profiles.")
                    adata_sc.obs[col] = adata_sc.obs[col].astype('category')
            
            # aggregate_profiles works best if X is CSR
            if not sp.isspmatrix_csr(adata_sc.X):
                logger.info("Converting adata_sc.X to CSR format for aggregate_profiles.")
                adata_sc.X = adata_sc.X.tocsr()

            adata_pb = sc.experimental.pp.aggregate_profiles(adata_sc, group_by_cols, mode='mean')
            
            if 'n_obs' in adata_pb.obs.columns: # Scanpy 1.7 to <1.9
                adata_pb.obs.rename(columns={'n_obs': 'n_cells_in_pseudobulk'}, inplace=True)
            elif 'N' in adata_pb.obs.columns: # Scanpy 1.9+
                 adata_pb.obs.rename(columns={'N': 'n_cells_in_pseudobulk'}, inplace=True)


            temp_obs_list = []
            # The index of adata_pb.obs is a MultiIndex if multiple group_by_cols, or simple Index if one
            for group_identifier in adata_pb.obs.index:
                group_tuple = group_identifier if isinstance(group_identifier, tuple) else (group_identifier,)
                
                mask = pd.Series(True, index=adata_sc.obs.index)
                group_dict_current = {}
                for i, col_name in enumerate(group_by_cols):
                    mask &= (adata_sc.obs[col_name].astype(str) == str(group_tuple[i]))
                    group_dict_current[col_name] = group_tuple[i]

                if not mask.any():
                    logger.warning(f"No cells found for group {group_tuple}. This should not happen if aggregate_profiles worked. Skipping obs meta for this group.")
                    continue

                first_cell_meta = adata_sc.obs[mask].iloc[0].copy()
                obs_entry = {col: group_dict_current[col] for col in group_by_cols}
                
                cols_to_keep_consistent = ['subject_id', 'tissue', 'sex', 'age', 'dataset', 'data_type',
                                           'species', 'species_ontology', 'tissue_ontology',
                                           'developmental_stage_ontology', 'assay_ontology', 'expression_unit']
                
                for col_meta in cols_to_keep_consistent:
                    if col_meta in first_cell_meta and col_meta not in obs_entry: # Don't overwrite groupby keys
                        obs_entry[col_meta] = first_cell_meta[col_meta]
                
                # Add n_cells from the aggregated data
                if 'n_cells_in_pseudobulk' in adata_pb.obs.columns:
                    obs_entry['n_cells_in_pseudobulk'] = adata_pb.obs.loc[group_identifier, 'n_cells_in_pseudobulk']
                else:
                    obs_entry['n_cells_in_pseudobulk'] = mask.sum() # Fallback if column name changed or missing
                temp_obs_list.append(obs_entry)
            
            new_obs_df = pd.DataFrame(temp_obs_list, index=adata_pb.obs.index)
            adata_pb.obs = new_obs_df
            adata_pb.obs.index.name = "pseudobulk_sample_id"

            logger.info("Pseudobulking with sc.experimental.pp.aggregate_profiles complete.")
            return adata_pb
            
    except Exception as e_agg:
        logger.warning(f"scanpy.experimental.pp.aggregate_profiles not used (version {sc.__version__} or error: {e_agg}). Using pandas groupby.")

    # Manual pandas groupby (fallback)
    logger.info("Using pandas groupby for pseudobulking.")
    expr_df = pd.DataFrame(adata_sc.X.toarray() if sp.issparse(adata_sc.X) else adata_sc.X,
                           index=adata_sc.obs.index,
                           columns=adata_sc.var.index)

    for col in group_by_cols: # Add group_by columns to expression dataframe
        expr_df[col] = adata_sc.obs[col].values
    
    logger.info("Calculating mean expression per group...")
    # Use observed=True if group_by_cols are categorical to avoid issues with unused categories
    are_categorical = all(pd.api.types.is_categorical_dtype(adata_sc.obs[col]) for col in group_by_cols)
    pseudobulk_expr = expr_df.groupby(group_by_cols, observed=are_categorical).mean()
    logger.info(f"Pseudobulk expression matrix shape: {pseudobulk_expr.shape}")

    pseudobulk_obs_list = []
    # pseudobulk_expr.index will be a MultiIndex if len(group_by_cols) > 1
    for group_identifier_tuple in pseudobulk_expr.index:
        group_identifier_tuple = group_identifier_tuple if isinstance(group_identifier_tuple, tuple) else (group_identifier_tuple,)
        
        mask = pd.Series(True, index=adata_sc.obs.index)
        group_dict = {}
        for i, col_name in enumerate(group_by_cols):
            mask &= (adata_sc.obs[col_name].astype(str) == str(group_identifier_tuple[i]))
            group_dict[col_name] = group_identifier_tuple[i]
        
        if not mask.any():
            logger.warning(f"No cells found for group {group_identifier_tuple} during pandas groupby. Skipping obs meta.")
            continue

        first_cell_meta = adata_sc.obs[mask].iloc[0].copy()
        obs_entry = {col: group_dict[col] for col in group_by_cols}
        
        cols_to_keep = ['subject_id', 'tissue', 'sex', 'age', 'dataset', 'data_type', 
                        'species', 'species_ontology', 'tissue_ontology', 
                        'developmental_stage_ontology', 'assay_ontology', 'expression_unit']
        
        for col in cols_to_keep:
            if col in first_cell_meta and col not in obs_entry:
                obs_entry[col] = first_cell_meta[col]
        
        obs_entry['n_cells_in_pseudobulk'] = mask.sum()
        pseudobulk_obs_list.append(obs_entry)

    pseudobulk_obs_df = pd.DataFrame(pseudobulk_obs_list)
    
    # Create the index for pseudobulk_obs_df based on the groupby keys
    if len(group_by_cols) > 1:
        # Match the MultiIndex from pseudobulk_expr directly if possible for index creation
        new_index_from_groups = pd.MultiIndex.from_tuples(
            [tuple(row[col] for col in group_by_cols) for _, row in pseudobulk_obs_df.iterrows()],
            names=group_by_cols
        )
        pseudobulk_obs_df.index = new_index_from_groups.map(lambda x: '_'.join(map(str,x)))
    else:
        pseudobulk_obs_df.index = pseudobulk_obs_df[group_by_cols[0]].astype(str)
    
    pseudobulk_obs_df.index.name = "pseudobulk_sample_id"

    # Align pseudobulk_obs_df with pseudobulk_expr index AFTER obs_df index is created
    # Convert pseudobulk_expr.index (which is MultiIndex) to string for matching obs_df string index
    pseudobulk_expr_str_index = pseudobulk_expr.index.map(lambda x: '_'.join(map(str,x)) if isinstance(x, tuple) else str(x))
    pseudobulk_expr.index = pseudobulk_expr_str_index # Make expr_df index also string

    # Reindex obs_df to match the order and content of pseudobulk_expr's new string index
    pseudobulk_obs_df = pseudobulk_obs_df.reindex(pseudobulk_expr.index)


    adata_pb = ad.AnnData(X=pseudobulk_expr.values,
                          obs=pseudobulk_obs_df,
                          var=adata_sc.var.copy())

    logger.info("Pandas groupby pseudobulking complete.")
    return adata_pb


def main():
    args = parse_arguments()

    logger.info(f"Loading single-cell data from: {args.input_h5ad}")
    try:
        adata_sc = sc.read_h5ad(args.input_h5ad)
    except Exception as e:
        logger.error(f"Failed to load input H5AD file {args.input_h5ad}: {e}")
        sys.exit(1)
        
    logger.info(f"Loaded single-cell AnnData: {adata_sc.n_obs} cells x {adata_sc.n_vars} genes")
    logger.info(f"Original .obs columns: {adata_sc.obs.columns.tolist()}")
    logger.info(f"Original .var columns: {adata_sc.var.columns.tolist()}")
    logger.info(f"Original .var index (sample): {adata_sc.var.index[:5].tolist()}")


    # --- 0. Standardize Gene IDs in adata_sc ---
    try:
        gene_ref_df = load_gene_reference_mapping(args.gene_reference_mapping_csv)
        adata_sc = standardize_adata_gene_ids(adata_sc, gene_ref_df)
        logger.info(f"After gene standardization: {adata_sc.n_obs} cells x {adata_sc.n_vars} genes")
        logger.info(f"Standardized .var index (sample): {adata_sc.var.index[:5].tolist()}")
    except Exception as e:
        logger.error(f"Error during gene ID standardization: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)

    # --- 1. Prepare and Standardize Cell-Level Metadata ---
    # ... (this part remains largely the same) ...
    try:
        adata_sc.obs = prepare_metadata_for_standardization(adata_sc.obs)
        
        logger.info("Standardizing cell-level metadata using rnaseq_utils...")
        mappings = load_mappings() # Uses internal METADATA_JSON_DIR
        
        logger.info(f"Before standardize_metadata - adata_sc.obs['tissue'] type: {type(adata_sc.obs['tissue'])}")
        if isinstance(adata_sc.obs['tissue'], pd.DataFrame):
            logger.info(f"Before standardize_metadata - adata_sc.obs['tissue'] columns: {adata_sc.obs['tissue'].columns.tolist()}")
        logger.info(f"Before standardize_metadata - adata_sc.obs['tissue'] head:\n{adata_sc.obs['tissue'].head()}")

        adata_sc.obs = standardize_metadata(
            adata_sc.obs, 
            "GTEx-snRNAseq", 
            mappings=mappings 
        )
        logger.info("Cell-level metadata standardization complete.")
        logger.info(f"Standardized .obs columns: {adata_sc.obs.columns.tolist()}")
    except Exception as e:
        logger.error(f"Error during metadata preparation/standardization: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    # ... (logging unique values) ...

    # --- 2. Perform Pseudobulking ---
    # ... (this part remains largely the same) ...
    group_by_cols = ['subject_id', 'cell_type']
    try:
        adata_pb = pseudobulk_anndata(adata_sc, group_by_cols) # Now uses adata_sc with standardized genes
    except Exception as e:
        logger.error(f"Error during pseudobulking: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    logger.info(f"Pseudobulked AnnData: {adata_pb.n_obs} pseudobulk samples x {adata_pb.n_vars} genes")


    # --- 3. Calculate Gene Variation ---
    # ... (this part remains largely the same) ...
    if adata_pb.n_obs > 1:
        logger.info("Calculating gene expression variance across pseudobulk samples...")
        if sp.issparse(adata_pb.X):
            adata_pb.var['pseudobulk_variance'] = np.var(adata_pb.X.toarray(), axis=0)
        else:
            adata_pb.var['pseudobulk_variance'] = np.var(adata_pb.X, axis=0)
        logger.info("Gene variance calculation complete.")
    else:
        logger.warning("Only one or zero pseudobulk samples generated. Skipping variance calculation.")
        adata_pb.var['pseudobulk_variance'] = np.nan if adata_pb.n_vars > 0 else np.array([])


    # --- 4. Add UNS Metadata ---
    # ... (this part remains largely the same) ...
    logger.info("Loading and applying dataset-specific metadata for pseudobulk AnnData...")
    gtex_sc_cfg = load_dataset_specific_metadata(args.metadata_json_dir, "gtex_scrnaseq")
    if gtex_sc_cfg:
        adata_pb = apply_dataset_specific_metadata(adata_pb, gtex_sc_cfg)
        logger.info("Applied 'gtex_scrnaseq_metadata.json'.")
    else:
        logger.warning("Could not load 'gtex_scrnaseq_metadata.json'. UNS metadata will be minimal.")
        adata_pb.uns['dataset_info'] = adata_pb.uns.get('dataset_info', {})
        adata_pb.uns['dataset_info']['source_dataset_label'] = 'GTEx-snRNAseq-pseudobulk'
        adata_pb.uns['dataset_info']['data_type'] = 'snRNA-seq-pseudobulk'
        adata_pb.uns['dataset_info']['expression_unit'] = 'Mean UMI counts per cell type'
        adata_pb.uns['harmonized_gencode_version'] = '24' 
        adata_pb.uns['harmonized_reference_genome'] = 'hg38'

    adata_pb.uns['pseudobulking_parameters'] = {
        'source_h5ad': os.path.basename(args.input_h5ad),
        'grouped_by': group_by_cols,
        'aggregation_method': 'sum_then_mean_via_pseudobulk', # Clarify it was sum for identical ENSG, then mean for pseudobulk
        'gene_id_standardization_source': os.path.basename(args.gene_reference_mapping_csv)
    }
    adata_pb.uns['creation_date'] = pd.Timestamp.now().strftime('%Y-%m-%d')

    # --- 5. Save Pseudobulk AnnData ---
    # ... (this part remains largely the same) ...
    logger.info(f"Saving pseudobulked data to: {args.output_h5ad} using {SAVE_METHOD} method.")
    if not robust_save_anndata(adata_pb, args.output_h5ad):
        logger.error(f"Failed to save the pseudobulk AnnData to {args.output_h5ad}.")
        sys.exit(1)
    else:
        logger.info(f"Successfully saved pseudobulk AnnData to {args.output_h5ad}")

    logger.info("GTEx single-cell processing complete.")

if __name__ == '__main__':
    main()