#!/usr/bin/env python3
"""
Create complete sample_id to ethnicity CSV mapping using all available data sources.
Now includes GTEx controlled-access race/ethnicity data!
"""

import scanpy as sc
import pandas as pd
import gzip
import json
import os

def load_gtex_phenotype_data():
    """Load GTEx controlled-access phenotype data with race/ethnicity."""
    gtex_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/gtex/metadata/phs000424.v10.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz"
    
    print("Loading GTEx controlled-access phenotype data...")
    
    # Read the compressed file directly with pandas, skipping comment lines
    import io
    
    with gzip.open(gtex_file, 'rt') as f:
        content = f.read()
    
    # Remove comment lines
    lines = content.split('\n')
    data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
    
    # Rejoin and create DataFrame
    clean_content = '\n'.join(data_lines)
    df = pd.read_csv(io.StringIO(clean_content), sep='\t')
    
    # Load race/ethnicity mappings from GTEx metadata file (correct for GTEx data!)
    try:
        gtex_metadata_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/gtex_metadata.json'
        if os.path.exists(gtex_metadata_file):
            with open(gtex_metadata_file, 'r') as f:
                gtex_metadata = json.load(f)
                
                # Get race mapping from GTEx metadata
                race_mapping_str = gtex_metadata.get('source_race_to_standard_label_map', {})
                race_mapping = {int(k): v for k, v in race_mapping_str.items()}
                
                # Get ethnicity mapping from GTEx metadata  
                ethnicity_mapping_str = gtex_metadata.get('source_ethnicity_to_standard_label_map', {})
                ethnicity_mapping = {int(k): v for k, v in ethnicity_mapping_str.items()}
                
                print(f"Loaded GTEx mappings from metadata file: {len(race_mapping)} race codes, {len(ethnicity_mapping)} ethnicity codes")
        else:
            raise FileNotFoundError("GTEx metadata file not found")
    except Exception as e:
        print(f"Warning: Could not load GTEx mappings from metadata file: {e}")
        print("Using fallback hardcoded mappings")
        
        # Fallback hardcoded mappings (based on GTEx data dictionary)
        race_mapping = {
            1: 'asian',
            2: 'black or african american',
            3: 'white',  # CRITICAL FIX: GTEx code 3 = white (not native hawaiian!)
            4: 'american indian or alaska native',
            5: 'native hawaiian or other pacific islander',
            98: 'unknown or not reported',
            99: 'unknown or not reported'
        }
        
        ethnicity_mapping = {
            1: 'hispanic or latino',
            2: 'not hispanic or latino', 
            98: 'unknown or not reported',
            99: 'unknown or not reported'
        }
    
    # Process the data
    gtex_ethnicity = {}
    for _, row in df.iterrows():
        subject_id = row['SUBJID']
        race_code = row['RACE']
        ethnicity_code = row['ETHNCTY']
        
        # Determine final ethnicity (prioritize Hispanic/Latino, then race)
        if ethnicity_code == 1:  # Hispanic or Latino
            final_ethnicity = 'hispanic or latino'
        else:
            final_ethnicity = race_mapping.get(race_code, 'unknown or not reported')
        
        gtex_ethnicity[subject_id] = final_ethnicity
    
    print(f"  Loaded ethnicity data for {len(gtex_ethnicity):,} GTEx subjects")
    
    # Show distribution
    ethnicity_counts = pd.Series(list(gtex_ethnicity.values())).value_counts()
    print(f"  GTEx ethnicity distribution:")
    for ethnicity, count in ethnicity_counts.items():
        print(f"    {ethnicity}: {count}")
    
    return gtex_ethnicity

def create_complete_ethnicity_mapping():
    """Create complete sample_id to ethnicity mapping with GTEx data."""
    
    # Load GTEx ethnicity data
    gtex_ethnicity_map = load_gtex_phenotype_data()
    
    datasets = {
        'ADNI': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/adni_standardized_preprocessed.h5ad',
        'ENCODE': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/encode_standardized_preprocessed.h5ad',
        'GTEx': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/gtex_standardized_preprocessed.h5ad',
        'MAGE': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/preprocessed_data/latest/mage_standardized_preprocessed.h5ad'
    }
    
    all_mappings = []
    
    for dataset_name, file_path in datasets.items():
        print(f"\nProcessing {dataset_name}...")
        
        adata = sc.read_h5ad(file_path)
        
        # Get sample IDs
        if 'subject_id' in adata.obs.columns:
            sample_ids = adata.obs['subject_id']
        elif 'donor_id' in adata.obs.columns:
            sample_ids = adata.obs['donor_id']
        elif 'sample_id' in adata.obs.columns:
            sample_ids = adata.obs['sample_id']
        else:
            sample_ids = adata.obs.index
        
        # Get ethnicity information and preserve original data - ENHANCED LOGIC
        original_data = []
        
        if dataset_name == 'ENCODE':
            # For ENCODE, use the raw 'ethnicity' column which has actual data
            if 'ethnicity' in adata.obs.columns:
                original_ethnicity = adata.obs['ethnicity'].copy()
                ethnicity = original_ethnicity.replace('European', 'white')
                original_data = original_ethnicity.tolist()
                print(f"  ENCODE: Using raw 'ethnicity' column")
            else:
                ethnicity = ['unknown'] * len(adata.obs)
                original_data = ['unknown'] * len(adata.obs)
                
        elif dataset_name == 'GTEx':
            # For GTEx, use the controlled-access phenotype data
            ethnicity_list = []
            original_data_list = []
            matched_count = 0
            for sample_id in sample_ids:
                # Extract subject ID from sample ID (format: GTEX-XXXXX-...)
                if isinstance(sample_id, str) and sample_id.startswith('GTEX-'):
                    subject_id = sample_id.split('-')[0] + '-' + sample_id.split('-')[1]  # GTEX-XXXXX
                    if subject_id in gtex_ethnicity_map:
                        mapped_ethnicity = gtex_ethnicity_map[subject_id]
                        ethnicity_list.append(mapped_ethnicity)
                        # For GTEx, the original data is the race/ethnicity codes that were mapped
                        original_data_list.append(mapped_ethnicity)  # This is already the processed form
                        matched_count += 1
                    else:
                        ethnicity_list.append('unknown or not reported')
                        original_data_list.append('unknown or not reported')
                else:
                    ethnicity_list.append('unknown or not reported')
                    original_data_list.append('unknown or not reported')
            
            ethnicity = ethnicity_list
            original_data = original_data_list
            print(f"  GTEx: Matched {matched_count:,}/{len(adata.obs):,} samples to controlled-access data")
            
        elif dataset_name == 'MAGE':
            # For MAGE, use proper population code mapping from metadata
            # Load MAGE metadata for population mappings
            try:
                mage_metadata_file = '/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/metadata/json/mage_metadata.json'
                with open(mage_metadata_file, 'r') as f:
                    mage_metadata = json.load(f)
                    pop_to_race_map = mage_metadata.get('pop_to_race_map', {})
                    print(f"  MAGE: Loaded {len(pop_to_race_map)} population mappings from metadata")
            except Exception as e:
                print(f"  MAGE: Warning - could not load metadata: {e}, using fallback")
                pop_to_race_map = {}
            
            # Preserve original population codes
            if 'population_code_1000g' in adata.obs.columns:
                original_data = adata.obs['population_code_1000g'].tolist()
                print(f"  MAGE: Found population_code_1000g for original data")
                
                # Map population codes to specific ethnicities using metadata
                ethnicity_list = []
                for pop_code in original_data:
                    if pop_code in pop_to_race_map:
                        ethnicity_list.append(pop_to_race_map[pop_code])
                    else:
                        ethnicity_list.append('unknown or not reported')
                ethnicity = ethnicity_list
                
            else:
                print(f"  MAGE: No population_code_1000g found, checking other columns")
                # Fall back to other available data
                if 'self_reported_ethnicity' in adata.obs.columns:
                    original_data = adata.obs['self_reported_ethnicity'].tolist()
                    ethnicity = adata.obs['self_reported_ethnicity']
                else:
                    original_data = ['unknown'] * len(adata.obs)
                    ethnicity = ['unknown'] * len(adata.obs)
        else:
            # For other datasets (ADNI), preserve original if available
            if 'self_reported_ethnicity' in adata.obs.columns:
                ethnicity = adata.obs['self_reported_ethnicity']
                original_data = ethnicity.tolist()  # For ADNI, this is likely already the original
            else:
                ethnicity = ['unknown'] * len(adata.obs)
                original_data = ['unknown'] * len(adata.obs)
        
        # Create mappings for this dataset with original data preserved
        for sample_id, eth, orig in zip(sample_ids, ethnicity, original_data):
            all_mappings.append({
                'sample_id': str(sample_id),
                'dataset': dataset_name,
                'ethnicity': str(eth),
                'original_ethnicity_race_population_data': str(orig)
            })
        
        print(f"  Added {len(adata.obs):,} samples from {dataset_name}")
    
    # Create DataFrame and normalize ethnicity categories
    df = pd.DataFrame(all_mappings)
    
    # Normalize ethnicity categories for consistency
    print(f"\n🔧 Normalizing ethnicity categories...")
    
    # Count changes made
    unknown_consolidated = (df['ethnicity'] == 'unknown').sum()
    if unknown_consolidated > 0:
        df.loc[df['ethnicity'] == 'unknown', 'ethnicity'] = 'unknown or not reported'
        print(f"  Consolidated {unknown_consolidated:,} 'unknown' → 'unknown or not reported'")
    
    df = df.sort_values(['dataset', 'sample_id'])
    
    output_file = 'sample_ethnicity_mapping_complete.csv'
    df.to_csv(output_file, index=False)
    
    print(f"\n✅ Complete ethnicity mapping saved to: {output_file}")
    print(f"📊 Total samples: {len(df):,}")
    
    # Show summary by dataset
    print("\n📋 Summary by dataset:")
    summary = df.groupby('dataset').agg({
        'sample_id': 'count',
        'ethnicity': lambda x: x.nunique()
    }).rename(columns={'sample_id': 'total_samples', 'ethnicity': 'unique_ethnicities'})
    print(summary.to_string())
    
    # Show ethnicity distribution
    print("\n🌍 Overall ethnicity distribution:")
    ethnicity_counts = df['ethnicity'].value_counts()
    for ethnicity, count in ethnicity_counts.head(15).items():
        print(f"  {ethnicity}: {count:,}")
    
    # Calculate improvement
    known_ethnicity = df[~df['ethnicity'].isin(['unknown or not reported', 'unknown'])]['ethnicity'].count()
    print(f"\n🎯 MAJOR IMPROVEMENT:")
    print(f"  Known ethnicity: {known_ethnicity:,}/{len(df):,} samples ({known_ethnicity/len(df)*100:.1f}%)")
    print(f"  Previous coverage: 6.6%")
    print(f"  New coverage: {known_ethnicity/len(df)*100:.1f}% (improvement: +{known_ethnicity/len(df)*100-6.6:.1f}%)")
    
    return output_file

if __name__ == '__main__':
    create_complete_ethnicity_mapping()