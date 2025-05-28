#!/usr/bin/env python3
"""
Create complete sample_id to ethnicity CSV mapping using all available data sources.
Now includes GTEx controlled-access race/ethnicity data!
"""

import scanpy as sc
import pandas as pd
import gzip

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
    
    # Map race codes to standard labels (from data dictionary)
    race_mapping = {
        1: 'asian',
        2: 'black or african american', 
        3: 'white',
        4: 'american indian or alaska native',
        98: 'unknown or not reported',
        99: 'unknown or not reported'
    }
    
    # Map ethnicity codes
    ethnicity_mapping = {
        0: 'not hispanic or latino',
        1: 'hispanic or latino',
        97: 'ceph',  # Special CEPH population
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
        'ADNI': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137/adni_standardized_preprocessed.h5ad',
        'ENCODE': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137/encode_standardized_preprocessed.h5ad',
        'GTEx': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137/gtex_standardized_preprocessed.h5ad',
        'MAGE': '/mnt/czi-sci-ai/intrinsic-variation-gene-ex-2/rnaseq/preprocessed_data/run_20250524_183137/mage_standardized_preprocessed.h5ad'
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
        
        # Get ethnicity information - ENHANCED LOGIC
        if dataset_name == 'ENCODE':
            # For ENCODE, use the raw 'ethnicity' column which has actual data
            if 'ethnicity' in adata.obs.columns:
                ethnicity = adata.obs['ethnicity'].copy()
                ethnicity = ethnicity.replace('European', 'white')
                print(f"  ENCODE: Using raw 'ethnicity' column")
            else:
                ethnicity = ['unknown'] * len(adata.obs)
                
        elif dataset_name == 'GTEx':
            # For GTEx, use the controlled-access phenotype data
            ethnicity_list = []
            matched_count = 0
            for sample_id in sample_ids:
                # Extract subject ID from sample ID (format: GTEX-XXXXX-...)
                if isinstance(sample_id, str) and sample_id.startswith('GTEX-'):
                    subject_id = sample_id.split('-')[0] + '-' + sample_id.split('-')[1]  # GTEX-XXXXX
                    if subject_id in gtex_ethnicity_map:
                        ethnicity_list.append(gtex_ethnicity_map[subject_id])
                        matched_count += 1
                    else:
                        ethnicity_list.append('unknown or not reported')
                else:
                    ethnicity_list.append('unknown or not reported')
            
            ethnicity = ethnicity_list
            print(f"  GTEx: Matched {matched_count:,}/{len(adata.obs):,} samples to controlled-access data")
            
        else:
            # For other datasets (ADNI, MAGE), use processed ethnicity
            if 'self_reported_ethnicity' in adata.obs.columns:
                ethnicity = adata.obs['self_reported_ethnicity']
            else:
                ethnicity = ['unknown'] * len(adata.obs)
        
        # Create mappings for this dataset
        for sample_id, eth in zip(sample_ids, ethnicity):
            all_mappings.append({
                'sample_id': str(sample_id),
                'dataset': dataset_name,
                'ethnicity': str(eth)
            })
        
        print(f"  Added {len(adata.obs):,} samples from {dataset_name}")
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(all_mappings)
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