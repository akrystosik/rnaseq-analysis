#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/metadata_consistency_validation.py
#!/usr/bin/env python3
"""
Metadata Consistency Validation Script

This script validates the consistency and completeness of metadata in standardized 
AnnData objects, ensuring proper standardization of categorical values and handling of
missing values.

Key validations:
- Completeness of core metadata fields
- Consistency of categorical values (sex, tissue types, etc.)
- Proper standardization of tissue/cell line identifiers
- Cross-dataset metadata consistency for matching samples/donors
"""

import os
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import logging
import argparse
import json
import matplotlib.pyplot as plt
import seaborn as sns

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('metadata_validator')

# Define paths
DEFAULT_STD_DATA_DIR = Path("/path/to/standardized_data")
DEFAULT_REFERENCE_METADATA = Path("/path/to/reference_metadata")
DEFAULT_OUTPUT_DIR = Path("./validation_results")

# Core metadata fields that should be present in all datasets
CORE_METADATA_FIELDS = [
    'sample_id',       # Unique identifier for each sample
    'subject_id',      # Identifier for the subject/donor
    'sex',             # Biological sex (standardized as 'male', 'female', 'unknown')
    'age',             # Age at sample collection
    'tissue',          # Tissue or brain region
    'dataset',         # Source dataset identifier
    'data_type',       # Type of data (RNA-seq, microarray, etc.)
    'expression_unit'  # Unit of measurement (TPM, FPKM, counts, etc.)
]

# Expected values for categorical fields
EXPECTED_CATEGORIES = {
    'sex': ['male', 'female', 'unknown'],
    'data_type': ['RNA-seq', 'microarray', 'single-cell', 'bulk'],
    'expression_unit': ['TPM', 'FPKM', 'counts', 'normalized_intensity']
}

def load_standardized_dataset(file_path):
    """Load a standardized dataset from an h5ad file."""
    logger.info(f"Loading standardized dataset from {file_path}")
    try:
        adata = ad.read_h5ad(file_path)
        logger.info(f"Loaded dataset with {adata.n_obs} samples and {adata.n_vars} genes")
        return adata
    except Exception as e:
        logger.error(f"Error loading dataset: {e}")
        return None

def load_reference_metadata(dataset_name, ref_dir):
    """Load reference metadata for validation."""
    ref_file = Path(ref_dir) / f"{dataset_name}_metadata.csv"
    
    if not ref_file.exists():
        logger.warning(f"Reference metadata file not found for {dataset_name}: {ref_file}")
        return None
    
    try:
        ref_metadata = pd.read_csv(ref_file)
        logger.info(f"Loaded reference metadata for {dataset_name} with {len(ref_metadata)} entries")
        return ref_metadata
    except Exception as e:
        logger.error(f"Error loading reference metadata: {e}")
        return None

def validate_metadata_completeness(adata, dataset_name):
    """Validate completeness of metadata fields."""
    results = {
        "dataset": dataset_name,
        "sample_count": adata.n_obs,
        "field_presence": {},
        "missing_core_fields": [],
        "field_completeness": {}
    }
    
    # Check presence of core metadata fields
    for field in CORE_METADATA_FIELDS:
        if field in adata.obs.columns:
            results["field_presence"][field] = True
            
            # Check for missing values
            missing_count = adata.obs[field].isna().sum()
            missing_percent = missing_count / adata.n_obs * 100
            
            results["field_completeness"][field] = {
                "present_count": adata.n_obs - missing_count,
                "missing_count": int(missing_count),
                "missing_percent": float(missing_percent)
            }
        else:
            results["field_presence"][field] = False
            results["missing_core_fields"].append(field)
    
    # Calculate overall metadata completeness score
    core_fields_present = sum(results["field_presence"].values())
    results["core_fields_present"] = core_fields_present
    results["core_fields_percent"] = core_fields_present / len(CORE_METADATA_FIELDS) * 100
    
    # Calculate value completeness score (considering missing values)
    if results["field_completeness"]:
        value_completeness = sum(
            field_stats["present_count"] / adata.n_obs 
            for field_stats in results["field_completeness"].values()
        ) / len(results["field_completeness"])
        
        results["value_completeness_percent"] = value_completeness * 100
    else:
        results["value_completeness_percent"] = 0
    
    return results

def validate_categorical_consistency(adata, dataset_name):
    """Validate consistency of categorical metadata fields."""
    results = {
        "dataset": dataset_name,
        "categorical_fields": {},
        "standardization_issues": []
    }
    
    # Check categorical fields against expected values
    for field, expected_values in EXPECTED_CATEGORIES.items():
        if field in adata.obs.columns:
            # Get actual values
            actual_values = adata.obs[field].dropna().unique().tolist()
            
            # Check if all values are within expected set
            unexpected_values = [val for val in actual_values if val not in expected_values]
            
            results["categorical_fields"][field] = {
                "expected_values": expected_values,
                "actual_values": actual_values,
                "unexpected_values": unexpected_values,
                "is_consistent": len(unexpected_values) == 0
            }
            
            # Flag inconsistencies
            if unexpected_values:
                results["standardization_issues"].append({
                    "field": field,
                    "issue": "unexpected_values",
                    "details": f"Found unexpected values: {', '.join(str(v) for v in unexpected_values)}"
                })
    
    # Check tissue field specifically
    if 'tissue' in adata.obs.columns:
        # Check for consistency in case (capitalization)
        tissue_values = adata.obs['tissue'].dropna().unique()
        
        # Check for possible duplicates due to capitalization or spacing
        normalized_tissues = {}
        for tissue in tissue_values:
            normalized = str(tissue).lower().strip()
            if normalized in normalized_tissues:
                results["standardization_issues"].append({
                    "field": "tissue",
                    "issue": "inconsistent_casing",
                    "details": f"Possible duplicate tissues: '{tissue}' and '{normalized_tissues[normalized]}'"
                })
            else:
                normalized_tissues[normalized] = tissue
    
    # Calculate overall consistency score
    consistent_fields = sum(
        1 for field_info in results["categorical_fields"].values() 
        if field_info["is_consistent"]
    )
    
    if results["categorical_fields"]:
        results["categorical_consistency_percent"] = consistent_fields / len(results["categorical_fields"]) * 100
    else:
        results["categorical_consistency_percent"] = 0
    
    return results

def validate_against_reference(adata, ref_metadata, dataset_name):
    """Validate metadata against reference sources."""
    if ref_metadata is None:
        return {
            "dataset": dataset_name,
            "reference_validation": "skipped",
            "reason": "No reference metadata available"
        }
    
    results = {
        "dataset": dataset_name,
        "matched_samples": 0,
        "field_consistency": {},
        "major_discrepancies": []
    }
    
    # Identify key fields for matching samples
    if 'sample_id' in ref_metadata.columns and 'sample_id' in adata.obs.columns:
        match_field = 'sample_id'
    elif 'subject_id' in ref_metadata.columns and 'subject_id' in adata.obs.columns:
        match_field = 'subject_id'
    else:
        return {
            "dataset": dataset_name,
            "reference_validation": "skipped",
            "reason": "No common identifier field for matching samples"
        }
    
    # Match samples between standardized data and reference
    matched_samples = 0
    field_matches = {}
    
    # Get common fields for comparison
    common_fields = [col for col in ref_metadata.columns if col in adata.obs.columns]
    
    # Remove the match field from comparison fields
    if match_field in common_fields:
        common_fields.remove(match_field)
    
    # Initialize field match counters
    for field in common_fields:
        field_matches[field] = {'matched': 0, 'mismatched': 0, 'mismatches': []}
    
    # Compare each sample
    for std_idx, std_sample_id in enumerate(adata.obs[match_field]):
        # Find this sample in the reference data
        ref_rows = ref_metadata[ref_metadata[match_field] == std_sample_id]
        
        if len(ref_rows) == 0:
            continue
            
        matched_samples += 1
        ref_row = ref_rows.iloc[0]
        
        # Compare values for each common field
        for field in common_fields:
            std_value = adata.obs[field].iloc[std_idx]
            ref_value = ref_row[field]
            
            # Handle categorical fields like sex differently - allow case-insensitive comparison
            if field == 'sex':
                if str(std_value).lower() == str(ref_value).lower():
                    field_matches[field]['matched'] += 1
                else:
                    field_matches[field]['mismatched'] += 1
                    # Record this mismatch
                    field_matches[field]['mismatches'].append({
                        match_field: std_sample_id,
                        'standardized': std_value,
                        'reference': ref_value
                    })
                    
                    # Flag major discrepancies
                    if field == 'sex' and pd.notna(std_value) and pd.notna(ref_value):
                        results["major_discrepancies"].append({
                            "field": field,
                            "sample": std_sample_id,
                            "standardized": std_value,
                            "reference": ref_value
                        })
            else:
                # For other fields, do direct comparison
                if pd.isna(std_value) and pd.isna(ref_value):
                    # Both are NaN, consider a match
                    field_matches[field]['matched'] += 1
                elif pd.notna(std_value) and pd.notna(ref_value) and std_value == ref_value:
                    field_matches[field]['matched'] += 1
                else:
                    field_matches[field]['mismatched'] += 1
                    # Record this mismatch
                    field_matches[field]['mismatches'].append({
                        match_field: std_sample_id,
                        'standardized': std_value,
                        'reference': ref_value
                    })
    
    # Calculate match percentages
    results["matched_samples"] = matched_samples
    
    for field, counts in field_matches.items():
        if matched_samples > 0:
            match_percent = counts['matched'] / matched_samples * 100
        else:
            match_percent = 0
            
        results["field_consistency"][field] = {
            "matched": counts['matched'],
            "mismatched": counts['mismatched'],
            "match_percent": match_percent,
            # Limit number of mismatches recorded in the output
            "mismatches": counts['mismatches'][:10] if len(counts['mismatches']) > 10 else counts['mismatches']
        }
    
    # Calculate overall consistency score
    if field_matches:
        avg_match_percent = np.mean([
            stats["match_percent"] for stats in results["field_consistency"].values()
        ])
        results["overall_consistency_percent"] = avg_match_percent
    else:
        results["overall_consistency_percent"] = 0
    
    return results

def validate_cross_dataset_consistency(datasets_dict):
    """Validate metadata consistency across different datasets."""
    # Skip if fewer than 2 datasets
    if len(datasets_dict) < 2:
        return {
            "cross_dataset_validation": "skipped",
            "reason": "Need at least 2 datasets for cross-validation"
        }
    
    results = {
        "dataset_count": len(datasets_dict),
        "common_donors": [],
        "tissue_encoding_consistency": {},
        "sex_encoding_consistency": {}
    }
    
    # Find donors that appear in multiple datasets
    donor_datasets = {}
    
    for dataset_name, adata in datasets_dict.items():
        # Find donor/subject column
        subject_col = next((col for col in ['subject_id', 'donor_id', 'donor'] if col in adata.obs.columns), None)
        
        if subject_col:
            for donor_id in adata.obs[subject_col].unique():
                if donor_id not in donor_datasets:
                    donor_datasets[donor_id] = []
                donor_datasets[donor_id].append(dataset_name)
    
    # Filter for donors in multiple datasets
    multi_dataset_donors = {
        donor_id: dataset_list 
        for donor_id, dataset_list in donor_datasets.items() 
        if len(dataset_list) > 1 and pd.notna(donor_id) and donor_id != ""
    }
    
    if not multi_dataset_donors:
        return {
            "cross_dataset_validation": "skipped",
            "reason": "No donors found in multiple datasets"
        }
    
    results["common_donors"] = [
        {"donor_id": donor_id, "datasets": dataset_list}
        for donor_id, dataset_list in multi_dataset_donors.items()
    ]
    
    # Check metadata consistency for common donors
    for donor_id, dataset_list in multi_dataset_donors.items():
        donor_info = {"donor_id": donor_id, "datasets": {}}
        
        for dataset_name in dataset_list:
            adata = datasets_dict[dataset_name]
            subject_col = next((col for col in ['subject_id', 'donor_id', 'donor'] if col in adata.obs.columns), None)
            
            if not subject_col:
                continue
                
            # Get samples for this donor
            donor_samples = adata.obs[adata.obs[subject_col] == donor_id]
            
            if donor_samples.empty:
                continue
                
            donor_info["datasets"][dataset_name] = {
                "sample_count": len(donor_samples)
            }
            
            # Check sex consistency if available
            if 'sex' in donor_samples.columns:
                sex_values = donor_samples['sex'].unique()
                donor_info["datasets"][dataset_name]["sex"] = sex_values[0] if len(sex_values) == 1 else sex_values.tolist()
            
            # Check tissue listing if available
            if 'tissue' in donor_samples.columns:
                tissues = donor_samples['tissue'].unique().tolist()
                donor_info["datasets"][dataset_name]["tissues"] = tissues
        
        # Check for inconsistencies across datasets
        # Sex consistency
        sex_values = {}
        for dataset_name, info in donor_info["datasets"].items():
            if "sex" in info:
                sex = info["sex"]
                if isinstance(sex, list):
                    for s in sex:
                        sex_values[s] = sex_values.get(s, 0) + 1
                else:
                    sex_values[sex] = sex_values.get(sex, 0) + 1
        
        if len(sex_values) > 1:
            # Found inconsistent sex values
            results["sex_encoding_consistency"][donor_id] = {
                "issue": "inconsistent_values",
                "values": sex_values
            }
    
    # Check tissue naming consistency across datasets
    all_tissues = {}
    
    for dataset_name, adata in datasets_dict.items():
        if 'tissue' in adata.obs.columns:
            dataset_tissues = adata.obs['tissue'].unique()
            all_tissues[dataset_name] = sorted(t for t in dataset_tissues if pd.notna(t) and t != "")
    
    # Find potential tissue naming inconsistencies
    tissue_mapping = {}
    
    for dataset1, tissues1 in all_tissues.items():
        for tissue1 in tissues1:
            normalized1 = tissue1.lower().strip()
            
            for dataset2, tissues2 in all_tissues.items():
                if dataset1 == dataset2:
                    continue
                    
                for tissue2 in tissues2:
                    normalized2 = tissue2.lower().strip()
                    
                    # Check for similar tissue names
                    if normalized1 == normalized2 and tissue1 != tissue2:
                        # Found potential naming inconsistency
                        key = f"{normalized1}"
                        if key not in tissue_mapping:
                            tissue_mapping[key] = {}
                        
                        if dataset1 not in tissue_mapping[key]:
                            tissue_mapping[key][dataset1] = []
                        
                        if dataset2 not in tissue_mapping[key]:
                            tissue_mapping[key][dataset2] = []
                        
                        if tissue1 not in tissue_mapping[key][dataset1]:
                            tissue_mapping[key][dataset1].append(tissue1)
                        
                        if tissue2 not in tissue_mapping[key][dataset2]:
                            tissue_mapping[key][dataset2].append(tissue2)
    
    results["tissue_encoding_consistency"] = tissue_mapping
    
    return results

def generate_metadata_plots(datasets_dict, output_dir):
    """Generate visualization plots for metadata."""
    plot_dir = output_dir / "metadata_plots"
    os.makedirs(plot_dir, exist_ok=True)
    
    # 1. Sex distribution by dataset
    plt.figure(figsize=(12, 6))
    sex_data = []
    
    for dataset_name, adata in datasets_dict.items():
        if 'sex' in adata.obs.columns:
            # Count sex distribution
            sex_counts = adata.obs['sex'].value_counts()
            
            # Add to data for plotting
            for sex, count in sex_counts.items():
                sex_data.append({
                    'Dataset': dataset_name,
                    'Sex': sex,
                    'Count': count
                })
    
    if sex_data:
        sex_df = pd.DataFrame(sex_data)
        
        # Plot
        ax = sns.barplot(x='Dataset', y='Count', hue='Sex', data=sex_df)
        plt.title('Sex Distribution by Dataset')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        # Save plot
        sex_plot_file = plot_dir / "sex_distribution.png"
        plt.savefig(sex_plot_file)
        plt.close()
    
    # 2. Tissue count by dataset
    plt.figure(figsize=(10, 6))
    tissue_counts = []
    
    for dataset_name, adata in datasets_dict.items():
        if 'tissue' in adata.obs.columns:
            unique_tissues = len(adata.obs['tissue'].unique())
            tissue_counts.append({
                'Dataset': dataset_name,
                'Unique Tissues': unique_tissues
            })
    
    if tissue_counts:
        tissue_df = pd.DataFrame(tissue_counts)
        
        # Plot
        ax = sns.barplot(x='Dataset', y='Unique Tissues', data=tissue_df)
        plt.title('Number of Unique Tissues by Dataset')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        # Save plot
        tissue_plot_file = plot_dir / "tissue_counts.png"
        plt.savefig(tissue_plot_file)
        plt.close()
    
    # 3. Metadata completeness heatmap
    plt.figure(figsize=(12, 8))
    completeness_data = []
    
    for dataset_name, adata in datasets_dict.items():
        for field in CORE_METADATA_FIELDS:
            if field in adata.obs.columns:
                non_missing = adata.obs[field].notna().sum()
                percent_complete = non_missing / len(adata.obs) * 100
            else:
                percent_complete = 0
                
            completeness_data.append({
                'Dataset': dataset_name,
                'Field': field,
                'Percent Complete': percent_complete
            })
    
    if completeness_data:
        completeness_df = pd.DataFrame(completeness_data)
        
        # Convert to wide format for heatmap
        heatmap_data = completeness_df.pivot(index='Dataset', columns='Field', values='Percent Complete')
        
        # Plot
        ax = sns.heatmap(heatmap_data, annot=True, cmap='YlGnBu', fmt='.1f')
        plt.title('Metadata Completeness (% of samples with value)')
        plt.tight_layout()
        
        # Save plot
        completeness_plot_file = plot_dir / "metadata_completeness.png"
        plt.savefig(completeness_plot_file)
        plt.close()
    
    return plot_dir

def main():
    parser = argparse.ArgumentParser(description='Validate metadata consistency in standardized datasets')
    parser.add_argument('--std-dir', type=str, default=str(DEFAULT_STD_DATA_DIR), 
                        help='Directory containing standardized h5ad files')
    parser.add_argument('--ref-dir', type=str, default=str(DEFAULT_REFERENCE_METADATA), 
                        help='Directory containing reference metadata files')
    parser.add_argument('--datasets', type=str, default='all', 
                        help='Comma-separated list of datasets to validate (or "all")')
    parser.add_argument('--output-dir', type=str, default=str(DEFAULT_OUTPUT_DIR), 
                        help='Output directory for validation results')
    
    args = parser.parse_args()
    
    # Find standardized datasets
    std_data_dir = Path(args.std_dir)
    std_files = list(std_data_dir.glob("*_standardized.h5ad"))
    
    if not std_files:
        logger.error(f"No standardized datasets found in {std_data_dir}")
        return
    
    logger.info(f"Found {len(std_files)} standardized datasets")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Determine which datasets to validate
    if args.datasets.lower() == 'all':
        datasets_to_validate = [f.stem.replace('_standardized', '') for f in std_files]
    else:
        datasets_to_validate = [d.strip() for d in args.datasets.split(',')]
    
    # Load all datasets for validation
    datasets_dict = {}
    
    for dataset_name in datasets_to_validate:
        std_file = std_data_dir / f"{dataset_name}_standardized.h5ad"
        
        if not std_file.exists():
            logger.warning(f"Standardized file not found for {dataset_name}")
            continue
        
        # Load standardized dataset
        adata = load_standardized_dataset(std_file)
        if adata is not None:
            datasets_dict[dataset_name] = adata
    
    if not datasets_dict:
        logger.error("No datasets could be loaded")
        return
    
    # Individual dataset validations
    completeness_results = []
    categorical_results = []
    reference_results = []
    
    for dataset_name, adata in datasets_dict.items():
        logger.info(f"Validating metadata for {dataset_name}")
        
        # Validate metadata completeness
        completeness = validate_metadata_completeness(adata, dataset_name)
        completeness_results.append(completeness)
        
        # Validate categorical consistency
        categorical = validate_categorical_consistency(adata, dataset_name)
        categorical_results.append(categorical)
        
        # Load reference metadata if available
        ref_metadata = load_reference_metadata(dataset_name, args.ref_dir)
        
        # Validate against reference
        reference = validate_against_reference(adata, ref_metadata, dataset_name)
        reference_results.append(reference)
    
    # Cross-dataset validation
    cross_dataset_results = validate_cross_dataset_consistency(datasets_dict)
    
    # Generate visualizations
    plot_dir = generate_metadata_plots(datasets_dict, output_dir)
    
    # Combine all results
    all_results = {
        "metadata_completeness": completeness_results,
        "categorical_consistency": categorical_results,
        "reference_validation": reference_results,
        "cross_dataset_validation": cross_dataset_results,
        "visualization_plots": str(plot_dir)
    }
    
    # Save full results as JSON
    results_file = output_dir / "metadata_validation_results.json"
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    logger.info(f"Saved full validation results to {results_file}")
    
    # Create summary CSV
    summary_data = []
    
    for i, dataset_name in enumerate(datasets_dict.keys()):
        completeness = completeness_results[i]
        categorical = categorical_results[i]
        reference = reference_results[i]
        
        summary_data.append({
            "Dataset": dataset_name,
            "Sample Count": completeness["sample_count"],
            "Core Fields Present (%)": completeness["core_fields_percent"],
            "Value Completeness (%)": completeness["value_completeness_percent"],
            "Categorical Consistency (%)": categorical.get("categorical_consistency_percent", "N/A"),
            "Reference Consistency (%)": reference.get("overall_consistency_percent", "N/A"),
            "Missing Core Fields": ", ".join(completeness["missing_core_fields"]),
            "Standardization Issues": len(categorical.get("standardization_issues", [])),
            "Major Discrepancies": len(reference.get("major_discrepancies", []))
        })
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_dir / "metadata_validation_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    
    logger.info(f"Saved validation summary to {summary_file}")
    
    # Print brief summary to console
    print("\nMetadata Validation Summary:")
    for row in summary_data:
        print(f"\n{row['Dataset']} ({row['Sample Count']} samples):")
        print(f"  Core Fields Present: {row['Core Fields Present (%)']: .1f}%")
        print(f"  Value Completeness: {row['Value Completeness (%)']: .1f}%")
        print(f"  Categorical Consistency: {row['Categorical Consistency (%)']: .1f}%")
        if row['Missing Core Fields']:
            print(f"  Missing Core Fields: {row['Missing Core Fields']}")
        if row['Standardization Issues'] > 0:
            print(f"  Standardization Issues: {row['Standardization Issues']}")
        if row['Major Discrepancies'] > 0:
            print(f"  Major Discrepancies: {row['Major Discrepancies']}")

if __name__ == "__main__":
    main()