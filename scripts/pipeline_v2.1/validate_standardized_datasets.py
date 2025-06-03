#!/usr/bin/env python3
"""
Standalone Validation Script for Standardized RNA-seq Datasets

This script validates standardized RNA-seq datasets against the required standards
and generates a comprehensive report.

Usage:
  python validate_standardized_datasets.py --input-dir /path/to/standardized/data \
                                          --output-file /path/to/validation_report.json
"""

import os
import sys
import json
import argparse
import logging
import scanpy as sc
import pandas as pd
from pathlib import Path
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('dataset_validator')

def extract_dataset_type(file_path):
    """Extract dataset type from filename."""
    file_name = os.path.basename(file_path)
    parts = file_name.split('_')
    if len(parts) > 0:
        return parts[0].lower()
    return "unknown"

def validate_dataset(file_path):
    """
    Validate a standardized dataset against requirements.
    
    Args:
        file_path: Path to the standardized dataset file
        
    Returns:
        Dictionary with validation results
    """
    try:
        file_path = Path(file_path)
        dataset_type = extract_dataset_type(file_path)
        logger.info(f"Validating {dataset_type} dataset: {file_path}")
        
        if not file_path.exists():
            logger.error(f"File not found: {file_path}")
            return {
                'dataset': dataset_type,
                'file': str(file_path),
                'status': 'error',
                'message': 'File not found',
                'validations': {}
            }
        
        # Load the dataset
        adata = sc.read_h5ad(file_path)
        
        # Initialize validation results
        validation_results = {
            'dataset': dataset_type,
            'file': str(file_path),
            'status': 'passed',
            'n_samples': adata.n_obs,
            'n_genes': adata.n_vars,
            'validations': {}
        }
        
        # Check harmonized GENCODE version
        gencode_version = adata.uns.get('harmonized_gencode_version', None)
        if gencode_version is not None:
            gencode_version = str(gencode_version).replace('v', '')
            if gencode_version == '24':
                validation_results['validations']['gencode_version'] = {
                    'status': 'passed',
                    'value': gencode_version
                }
            else:
                validation_results['validations']['gencode_version'] = {
                    'status': 'failed',
                    'value': gencode_version,
                    'expected': '24'
                }
                validation_results['status'] = 'failed'
        else:
            validation_results['validations']['gencode_version'] = {
                'status': 'missing',
                'expected': '24'
            }
            validation_results['status'] = 'failed'
        
        # Check harmonized reference genome
        genome_version = adata.uns.get('harmonized_reference_genome', None)
        if genome_version is not None:
            if genome_version in ['hg38', 'GRCh38']:
                validation_results['validations']['reference_genome'] = {
                    'status': 'passed',
                    'value': genome_version
                }
            else:
                validation_results['validations']['reference_genome'] = {
                    'status': 'failed',
                    'value': genome_version,
                    'expected': 'hg38/GRCh38'
                }
                validation_results['status'] = 'failed'
        else:
            validation_results['validations']['reference_genome'] = {
                'status': 'missing',
                'expected': 'hg38/GRCh38'
            }
            validation_results['status'] = 'failed'
        
        # Check observation metadata fields
        metadata_fields = {
            'tissue': {
                'ontology_field': 'tissue_ontology',
                'ontology_prefix': 'UBERON:',
                'importance': 'critical'
            },
            'sex': {
                'values': ['male', 'female', 'unknown'],
                'importance': 'important'
            },
            'species': {
                'ontology_field': 'species_ontology',
                'ontology_prefix': 'NCBITaxon:',
                'importance': 'important'
            },
            'data_type': {
                'values': ['RNA-seq', 'microarray'],
                'importance': 'important'
            },
            'assay_ontology': {
                'ontology_prefix': 'EFO:',
                'importance': 'important'
            },
            'cell_type': {
                'ontology_field': 'cell_type_ontology_term_id',
                'ontology_prefix': 'CL:',
                'importance': 'important'
            }
        }


        
        for field, config in metadata_fields.items():
            if field in adata.obs.columns:
                # Check for missing values
                missing_count = adata.obs[field].isna().sum()
                missing_percentage = (missing_count / adata.n_obs) * 100
                
                # Check for ontology field
                if 'ontology_field' in config and config['ontology_field'] in adata.obs.columns:
                    # Check if ontology values are valid
                    ontology_field = config['ontology_field']
                    ontology_prefix = config.get('ontology_prefix', '')
                    
                    # Count values with correct prefix
                    valid_values = adata.obs[ontology_field].astype(str).str.startswith(ontology_prefix)
                    valid_count = valid_values.sum()
                    valid_percentage = (valid_count / adata.n_obs) * 100
                    
                    validation_results['validations'][field] = {
                        'status': 'passed' if valid_percentage >= 90 else 'warning' if valid_percentage >= 70 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'valid_percentage': float(valid_percentage),
                        'has_ontology': True
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
                
                # Check for enumerated values
                elif 'values' in config:
                    valid_values = adata.obs[field].isin(config['values'])
                    valid_count = valid_values.sum()
                    valid_percentage = (valid_count / adata.n_obs) * 100
                    
                    validation_results['validations'][field] = {
                        'status': 'passed' if valid_percentage >= 90 else 'warning' if valid_percentage >= 70 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'valid_percentage': float(valid_percentage),
                        'has_ontology': False
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
                
                # Simple presence check
                else:
                    validation_results['validations'][field] = {
                        'status': 'passed' if missing_percentage <= 10 else 'warning' if missing_percentage <= 30 else 'failed',
                        'missing_percentage': float(missing_percentage),
                        'has_ontology': False
                    }
                    
                    if validation_results['validations'][field]['status'] == 'failed' and config['importance'] == 'critical':
                        validation_results['status'] = 'failed'
                    elif validation_results['validations'][field]['status'] == 'warning' and validation_results['status'] == 'passed':
                        validation_results['status'] = 'warning'
            else:
                validation_results['validations'][field] = {
                    'status': 'missing',
                    'importance': config['importance']
                }
                
                if config['importance'] == 'critical':
                    validation_results['status'] = 'failed'
        
        # Additional validation: Check gene IDs format
        if adata.n_vars > 0:
            # Check if gene IDs follow Ensembl format
            sample_genes = adata.var_names[:10].tolist()
            ensembl_format = all(str(g).startswith('ENSG') for g in sample_genes)
            
            validation_results['validations']['gene_id_format'] = {
                'status': 'passed' if ensembl_format else 'failed',
                'value': 'Ensembl' if ensembl_format else 'Unknown'
            }
            
            if not ensembl_format and validation_results['status'] == 'passed':
                validation_results['status'] = 'warning'
        
        logger.info(f"Validation completed for {dataset_type}: {validation_results['status']}")
        return validation_results
        
    except Exception as e:
        logger.error(f"Error validating dataset {file_path}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {
            'dataset': extract_dataset_type(file_path),
            'file': str(file_path),
            'status': 'error',
            'message': str(e),
            'validations': {}
        }

def generate_report(validation_results, output_file):
    """
    Generate a validation report from the results.
    
    Args:
        validation_results: List of validation result dictionaries
        output_file: Path to save the report
        
    Returns:
        Path to the report file
    """
    try:
        output_file = Path(output_file)
        
        # Create summary statistics
        summary = {
            'timestamp': datetime.now().isoformat(),
            'datasets_validated': len(validation_results),
            'datasets_passed': sum(1 for r in validation_results if r['status'] == 'passed'),
            'datasets_warning': sum(1 for r in validation_results if r['status'] == 'warning'),
            'datasets_failed': sum(1 for r in validation_results if r['status'] == 'failed'),
            'datasets_error': sum(1 for r in validation_results if r['status'] == 'error'),
            'total_samples': sum(r.get('n_samples', 0) for r in validation_results),
            'total_genes': sum(r.get('n_genes', 0) for r in validation_results),
            'dataset_results': validation_results
        }
        
        # Ensure output directory exists
        output_file.parent.mkdir(exist_ok=True, parents=True)
        
        # Write report to file
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Also generate a more human-readable HTML report
        html_file = output_file.with_suffix('.html')
        
        # Simple HTML report
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>RNA-seq Standardization Validation Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                h1, h2, h3 {{ color: #333; }}
                .summary {{ background-color: #f5f5f5; padding: 10px; border-radius: 5px; }}
                .dataset {{ margin-bottom: 20px; padding: 10px; border: 1px solid #ddd; }}
                .passed {{ background-color: #dff0d8; }}
                .warning {{ background-color: #fcf8e3; }}
                .failed {{ background-color: #f2dede; }}
                .error {{ background-color: #f5f5f5; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <h1>RNA-seq Standardization Validation Report</h1>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            
            <div class="summary">
                <h2>Summary</h2>
                <p>Datasets validated: {summary['datasets_validated']}</p>
                <p>Datasets passed: {summary['datasets_passed']}</p>
                <p>Datasets with warnings: {summary['datasets_warning']}</p>
                <p>Datasets failed: {summary['datasets_failed']}</p>
                <p>Datasets with errors: {summary['datasets_error']}</p>
                <p>Total samples: {summary['total_samples']}</p>
                <p>Total genes: {summary['total_genes']}</p>
            </div>
            
            <h2>Dataset Results</h2>
        """
        
        # Add sections for each dataset
        for result in validation_results:
            status_class = result['status']
            dataset_name = result['dataset'].upper()
            
            html_content += f"""
            <div class="dataset {status_class}">
                <h3>{dataset_name} - {result['status'].upper()}</h3>
                <p>File: {result.get('file', 'N/A')}</p>
                <p>Samples: {result.get('n_samples', 'N/A')}</p>
                <p>Genes: {result.get('n_genes', 'N/A')}</p>
            """
            
            if result['status'] == 'error':
                html_content += f"""
                <p>Error: {result.get('message', 'Unknown error')}</p>
                """
            else:
                html_content += """
                <h4>Validation Results</h4>
                <table>
                    <tr>
                        <th>Field</th>
                        <th>Status</th>
                        <th>Details</th>
                    </tr>
                """
                
                for field, details in result.get('validations', {}).items():
                    status = details.get('status', 'unknown')
                    
                    if status == 'passed':
                        details_str = f"Value: {details.get('value', 'N/A')}"
                    elif status == 'failed':
                        details_str = f"Value: {details.get('value', 'N/A')}, Expected: {details.get('expected', 'N/A')}"
                    elif status == 'missing':
                        details_str = f"Expected: {details.get('expected', 'N/A')}"
                    else:
                        # For fields with percentages
                        if 'missing_percentage' in details:
                            details_str = f"Missing: {details['missing_percentage']:.1f}%"
                        
                        if 'valid_percentage' in details:
                            details_str += f", Valid: {details['valid_percentage']:.1f}%"
                    
                    html_content += f"""
                    <tr>
                        <td>{field}</td>
                        <td>{status.upper()}</td>
                        <td>{details_str}</td>
                    </tr>
                    """
                
                html_content += """
                </table>
                """
            
            html_content += """
            </div>
            """
        
        html_content += """
        </body>
        </html>
        """
        
        # Write HTML report
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Report saved to {output_file}")
        logger.info(f"HTML report saved to {html_file}")
        
        return output_file
        
    except Exception as e:
        logger.error(f"Error generating report: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Validate Standardized RNA-seq Datasets')
    
    parser.add_argument('--input-dir', required=True, help='Directory containing standardized datasets')
    parser.add_argument('--output-file', help='Path to save the validation report JSON file')
    parser.add_argument('--file-pattern', default='*_standardized_v2.h5ad', 
                       help='File pattern to match standardized datasets')
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Find standardized datasets
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        sys.exit(1)
    
    # Determine output file
    if args.output_file:
        output_file = Path(args.output_file)
    else:
        output_file = input_dir / f"validation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    # Find datasets to validate
    h5ad_files = list(input_dir.glob(args.file_pattern))
    if not h5ad_files:
        logger.error(f"No standardized datasets found matching pattern '{args.file_pattern}' in {input_dir}")
        sys.exit(1)
    
    logger.info(f"Found {len(h5ad_files)} datasets to validate")
    
    # Validate each dataset
    validation_results = []
    for h5ad_file in h5ad_files:
        result = validate_dataset(h5ad_file)
        validation_results.append(result)
    
    # Generate report
    report_file = generate_report(validation_results, output_file)
    if report_file:
        logger.info(f"Validation completed. Report saved to {report_file}")
    else:
        logger.error("Failed to generate validation report")
        sys.exit(1)

if __name__ == '__main__':
    main()