
#!/usr/bin/env python3
"""
RNA-seq Analysis Pipeline Setup Script
-------------------------------------

This script sets up the RNA-seq analysis pipeline:
- Creates the necessary directory structure
- Copies the default configuration
- Ensures required packages are installed
- Verifies access to data directories
"""

import os
import sys
import shutil
import json
import argparse
from pathlib import Path
import subprocess
import pkg_resources

export PYTHONPATH=$PYTHONPATH:/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq_analysis_new/src


def check_packages():
    """Check for required packages"""
    required_packages = [
        'numpy',
        'pandas',
        'matplotlib',
        'seaborn',
        'requests'
    ]
    
    optional_packages = [
        'mygene'  # For enhanced gene mapping
    ]
    
    missing = []
    
    for package in required_packages:
        try:
            pkg_resources.get_distribution(package)
        except pkg_resources.DistributionNotFound:
            missing.append(package)
    
    if missing:
        print("Missing required packages:")
        for package in missing:
            print(f"  - {package}")
        print("\nInstall with:")
        print(f"  pip install {' '.join(missing)}")
        return False
    
    print("All required packages are installed.")
    
    # Check optional packages
    missing_optional = []
    for package in optional_packages:
        try:
            pkg_resources.get_distribution(package)
        except pkg_resources.DistributionNotFound:
            missing_optional.append(package)
    
    if missing_optional:
        print("\nOptional packages not installed:")
        for package in missing_optional:
            print(f"  - {package}")
        print("\nInstall with:")
        print(f"  pip install {' '.join(missing_optional)}")
        print("\nThese packages are optional but provide additional functionality.")
    
    return True

def create_directory_structure(base_dir):
    """Create the directory structure"""
    base_path = Path(base_dir)
    
    # Create main directories
    directories = [
        base_path,
        base_path / "raw_data" / "gene_quantification",
        base_path / "analysis",
        base_path / "metadata",
        base_path / "config"
    ]
    
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)
        print(f"Created directory: {directory}")
    
    return True

def create_config(base_dir):
    """Create the default configuration file"""
    base_path = Path(base_dir)
    config_file = base_path / "config" / "settings.json"
    
    # Default configuration
    config = {
        "base_dir": str(base_path),
        "raw_data_dir": str(base_path / "raw_data" / "gene_quantification"),
        "analysis_dir": str(base_path / "analysis"),
        "metadata_dir": str(base_path / "metadata"),
        "housekeeping_genes": [
            "ACTB", "GAPDH", "TUBB", "B2M", "PPIA", "TBP", "HPRT1", "RPL13A"
        ],
        "rna_seq_types": ["total", "polyA+", "polyA-"],
        "cell_lines": {
            "K562": {"experiment_id": "ENCSR792OIJ", "rna_type": "total"},
            "Panc1": {"experiment_id": "ENCSR128CYL", "rna_type": "total"},
            "NCI-H460": {"experiment_id": "ENCSR164OCT", "rna_type": "total"},
            "A549": {"experiment_id": "ENCSR414IGI", "rna_type": "total"},
            "K562_AEL": {"experiment_id": "ENCSR000AEL", "rna_type": "total"},
            "HepG2_total": {"experiment_id": "ENCSR245ATJ", "rna_type": "total"},
            "HepG2_polyA_plus": {"experiment_id": "ENCSR000CPD", "rna_type": "polyA+"},
            "HepG2_polyA_minus": {"experiment_id": "ENCSR000CPE", "rna_type": "polyA-"},
            "HepG2_total_ZGR": {"experiment_id": "ENCSR181ZGR", "rna_type": "total"},
            "Caki2_584JXD": {"experiment_id": "ENCSR584JXD", "rna_type": "total"},
            "A549_polyA_minus": {"experiment_id": "ENCSR000CQC", "rna_type": "polyA-"},
            "A549_polyA_plus": {"experiment_id": "ENCSR000CON", "rna_type": "polyA+"},
            "GM23248_BPP": {"experiment_id": "ENCSR797BPP", "rna_type": "total"}
        },
        "known_gene_mappings": {
            "ACTB": "ENSG00000075624",
            "GAPDH": "ENSG00000111640",
            "TUBB": "ENSG00000196230",
            "B2M": "ENSG00000166710",
            "PPIA": "ENSG00000196262",
            "TBP": "ENSG00000112592",
            "HPRT1": "ENSG00000165704",
            "RPL13A": "ENSG00000142541"
        }
    }
    
    # Save configuration
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"Created configuration file: {config_file}")
    return True

def copy_scripts(base_dir, script_dir):
    """Copy script files to the base directory"""
    base_path = Path(base_dir)
    scripts_path = base_path / "scripts"
    scripts_path.mkdir(exist_ok=True)
    
    # Source script directory
    src_dir = Path(script_dir)
    
    # Copy scripts
    script_files = [
        "analyze_cell_lines.py",
        "build_gene_mapping.py"
    ]
    
    for script in script_files:
        src_file = src_dir / script
        dst_file = scripts_path / script
        
        if src_file.exists():
            shutil.copy2(src_file, dst_file)
            # Make executable
            os.chmod(dst_file, 0o755)
            print(f"Copied script: {dst_file}")
        else:
            print(f"Warning: Could not find script {src_file}")
    
    return True

def verify_directories(base_dir):
    """Verify that all directories are accessible"""
    base_path = Path(base_dir)
    
    directories = [
        base_path,
        base_path / "raw_data" / "gene_quantification",
        base_path / "analysis",
        base_path / "metadata",
        base_path / "config"
    ]
    
    for directory in directories:
        if not directory.exists():
            print(f"Error: Directory {directory} does not exist")
            return False
        
        # Check if writable
        try:
            test_file = directory / ".test_write"
            with open(test_file, 'w') as f:
                f.write("test")
            test_file.unlink()  # Remove test file
        except (IOError, PermissionError) as e:
            print(f"Error: Directory {directory} is not writable: {e}")
            return False
    
    print("All directories are accessible and writable.")
    return True

def create_module_structure(base_dir, script_dir):
    """Create the module structure for the pipeline"""
    base_path = Path(base_dir)
    
    # Create module directories
    module_dirs = [
        base_path / "core",
        base_path / "utils",
        base_path / "visualization",
        base_path / "qc"
    ]
    
    for directory in module_dirs:
        directory.mkdir(exist_ok=True)
        # Create __init__.py
        init_file = directory / "__init__.py"
        if not init_file.exists():
            with open(init_file, 'w') as f:
                f.write("# Module initialization\n")
        print(f"Created module directory: {directory}")
    
    # Copy module files from script directory
    src_dir = Path(script_dir)
    
    # Map of source files to destination directories
    module_files = {
        "core/dataset_manager.py": base_path / "core",
        "core/gene_identifier.py": base_path / "core",
        "core/expression_analyzer.py": base_path / "core",
        "qc/basic_qc.py": base_path / "qc",
        "visualization/heatmap.py": base_path / "visualization"
    }
    
    for src_file, dst_dir in module_files.items():
        full_src = src_dir.parent / src_file
        full_dst = dst_dir / Path(src_file).name
        
        if full_src.exists():
            shutil.copy2(full_src, full_dst)
            print(f"Copied module file: {full_dst}")
        else:
            print(f"Warning: Could not find module file {full_src}")
    
    return True

def main():
    """Set up the RNA-seq analysis pipeline"""
    parser = argparse.ArgumentParser(description='Set up RNA-seq analysis pipeline')
    
    parser.add_argument('--base-dir', required=True,
                      help='Base directory for the pipeline')
    parser.add_argument('--script-dir', default=os.path.dirname(os.path.abspath(__file__)),
                      help='Directory containing the pipeline scripts')
    parser.add_argument('--skip-package-check', action='store_true',
                      help='Skip checking for required packages')
    
    args = parser.parse_args()
    
    print("\n===== Setting up RNA-seq Analysis Pipeline =====")
    
    # Check packages
    if not args.skip_package_check:
        if not check_packages():
            print("\nPlease install required packages and run setup again.")
            return 1
    
    # Create directory structure
    if not create_directory_structure(args.base_dir):
        print("\nFailed to create directory structure.")
        return 1
    
    # Create configuration
    if not create_config(args.base_dir):
        print("\nFailed to create configuration.")
        return 1
    
    # Create module structure
    if not create_module_structure(args.base_dir, args.script_dir):
        print("\nFailed to create module structure.")
        return 1
    
    # Copy scripts
    if not copy_scripts(args.base_dir, args.script_dir):
        print("\nFailed to copy scripts.")
        return 1
    
    # Verify directories
    if not verify_directories(args.base_dir):
        print("\nDirectory verification failed.")
        return 1
    
    print("\n===== RNA-seq Analysis Pipeline Setup Complete =====")
    print(f"Pipeline directory: {args.base_dir}")
    print("\nNext steps:")
    print("1. Build gene mapping database:")
    print(f"   python {args.base_dir}/scripts/build_gene_mapping.py --use-biomart")
    print("\n2. Run analysis:")
    print(f"   python {args.base_dir}/scripts/analyze_cell_lines.py")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())