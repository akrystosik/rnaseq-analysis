#!/usr/bin/env python3
"""
Wrapper script to save AnnData objects with proper string conversion.
Usage: python anndata_save_wrapper.py input_file output_file
"""

import sys
import os
import scanpy as sc
import pandas as pd
import numpy as np
import logging
import json # For complex object serialization to string

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('anndata_save_wrapper')

def convert_to_serializable(obj):
    """Convert dict values to HDF5-compatible native Python types or JSON strings."""
    if isinstance(obj, dict):
        return {convert_to_serializable(k): convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)): # Handle tuples as lists
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, (np.integer, int)):
        return int(obj)
    elif isinstance(obj, (np.floating, float)):
        if np.isnan(obj):
            return None  # Represent NaN as None, HDF5 handles None
        return float(obj)
    elif isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        if obj.dtype == 'object':
            # If array of objects, recursively serialize each element
            return [convert_to_serializable(x) for x in obj.tolist()]
        elif obj.size == 0: # Handle empty array
            return []
        else:
            # For numeric/boolean arrays, tolist() is generally safe
            # NaNs in float arrays become Python float NaNs, handled by HDF5
            return obj.tolist()
    elif isinstance(obj, (pd.Series, pd.DataFrame)):
        logger.warning(f"Attempting to serialize pandas {type(obj).__name__} in uns. Converting to JSON string.")
        try:
            # For DataFrames, convert to list of records (dicts), then serialize
            if isinstance(obj, pd.DataFrame):
                return json.dumps([convert_to_serializable(record) for record in obj.to_dict(orient='records')], default=str)
            # For Series, convert to list, then serialize
            return json.dumps(convert_to_serializable(obj.tolist()), default=str)
        except Exception as e:
            logger.error(f"Could not serialize pandas {type(obj).__name__} to JSON: {e}. Falling back to str().")
            return str(obj)
    elif pd.isna(obj): # Check for pandas NA types (like pd.NA, NaT) AFTER specific type checks
        return None # Represent as None
    elif isinstance(obj, str):
        return obj
    else:
        # Fallback for other types: attempt to convert to string.
        # If it's a complex custom object, this might be a lossy conversion.
        logger.debug(f"Converting type {type(obj)} to string for serialization.")
        try:
            # For some specific complex types, a targeted conversion or JSON dump might be better
            # Example: if obj is some specific BioPython object, etc.
            # For now, a general attempt to JSON dump if it's not basic, then string
            if not isinstance(obj, (int, float, bool, str, bytes)) and obj is not None:
                 try:
                     return json.dumps(obj, default=str) # Try JSON first for unknown complex types
                 except TypeError:
                     return str(obj) # Final fallback to string
            return str(obj)
        except Exception as e_str:
            logger.error(f"Failed to convert object of type {type(obj)} to string: {e_str}. Returning empty string.")
            return ""


def save_adata_safely(input_file, output_file):
    """Load AnnData and save it with proper string conversion."""
    try:
        logger.info(f"Loading AnnData from {input_file}")
        adata = sc.read_h5ad(input_file)

        logger.info(f"Loaded AnnData with {adata.n_obs} samples and {adata.n_vars} genes")

        # Make a copy of the original uns
        original_uns = adata.uns.copy()

        # Convert all uns values to serializable types
        logger.info("Converting uns dictionary to serializable types")
        adata.uns = convert_to_serializable(original_uns)

        # Save the AnnData object
        logger.info(f"Saving AnnData to {output_file}")
        adata.write(output_file)

        # Verify the save worked
        logger.info("Verifying saved file")
        test_adata = sc.read_h5ad(output_file)
        logger.info(f"Verification successful: {test_adata.n_obs} samples, {test_adata.n_vars} genes")

        return True

    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        logger.error(traceback.format_exc()) # Print full traceback for debugging
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python anndata_save_wrapper.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    success = save_adata_safely(input_file, output_file)

    if success:
        print(f"Successfully saved AnnData to {output_file}")
        sys.exit(0)
    else:
        print("Failed to save AnnData")
        sys.exit(1)