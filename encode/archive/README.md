# RNA-seq Gene Quantification Analysis Pipeline

A modular, maintainable pipeline for analyzing gene expression from RNA-seq experiments, with a focus on cross-dataset comparisons.

## Overview

This pipeline analyzes gene quantification data from ENCODE RNA-seq experiments to identify comparable datasets across cell lines and tissues. It uses a data-driven approach to expression analysis, focusing on housekeeping genes to establish robust comparisons.

## Features

- **Robust Gene Identification**: Improved string handling and Ensembl ID processing for accurate gene identification
- **Cross-Dataset Analysis**: Compare expression patterns between cell lines and tissue samples
- **Data-Driven Clustering**: Identify groups of related datasets based on expression similarity
- **Quality Control**: Automatically filter low-quality datasets
- **Visualization**: Generate informative heatmaps for gene expression and dataset correlations

## Directory Structure

```
/rnaseq_analysis/
├── README.md
├── config/               # Configuration settings
├── src/                  # Source code modules
│   ├── core/             # Core functionality
│   │   ├── dataset_manager.py    # Dataset handling
│   │   ├── gene_identifier.py    # Gene identification
│   │   └── expression_analyzer.py # Expression analysis
│   ├── qc/               # Quality control
│   │   └── basic_qc.py
│   ├── visualization/    # Visualization tools
│   │   └── heatmap.py
│   └── utils/            # Utility functions
├── scripts/              # Analysis scripts
│   ├── analyze_cell_lines.py
│   ├── build_gene_mapping.py
│   ├── generate_visualizations.py
│   └── setup.py
├── logs/                 # Log files
├── metadata/             # Cached metadata and gene mappings
├── raw_data/             # Raw data files
│   └── gene_quantification/
└── analysis/             # Analysis results
    └── YYYYMMDD/         # Date-stamped directories
        ├── results/
        └── visualizations/
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/rnaseq-analysis.git
   cd rnaseq-analysis
   ```

2. Set up the environment:
   ```bash
   pip install numpy pandas matplotlib seaborn requests
   ```

3. Optional but recommended packages:
   ```bash
   pip install mygene
   ```

4. Initialize the pipeline structure:
   ```bash
   mkdir -p src/config src/core src/qc src/visualization src/utils scripts logs metadata raw_data/gene_quantification analysis
   ```

## Usage

### 1. Setup the Pipeline

```bash
python scripts/setup.py --base-dir /path/to/rnaseq_analysis
```

### 2. Build Gene Mapping Database

```bash
python scripts/build_gene_mapping.py --use-biomart
```

### 3. Run Analysis

Basic analysis:
```bash
python scripts/analyze_cell_lines.py
```

Include ENTEx tissue datasets:
```bash
python scripts/analyze_cell_lines.py --include-entex
```

With log transformation:
```bash
python scripts/analyze_cell_lines.py --use-log
```

Analyze specific datasets:
```bash
python scripts/analyze_cell_lines.py --datasets K562 A549 HepG2_total
```

### 4. Generate Visualizations

If the pipeline visualizations don't appear, you can use the dedicated visualization script:
```bash
python scripts/simple_viz.py
```

### 5. Run the Complete Pipeline

Use the run_pipeline.sh script to execute the entire workflow:
```bash
./run_pipeline_fixed.sh --include-entex --use-log
```

## Key Improvements

This pipeline improves on the original version by:

1. **Modular Design**: Separating concerns for better maintainability
2. **Robust Gene Identification**: Handling different file formats and gene ID versions
3. **Data-Driven Analysis**: Letting the data guide the analysis instead of relying on complex marker definitions
4. **Improved Visualization**: Clear, informative visualizations of expression patterns

## Adding New Datasets

To add new cell lines or tissue samples, edit the `config/settings.json` file or use the DatasetManager:

```python
from src.core.dataset_manager import DatasetManager
dm = DatasetManager()
dm.add_cell_line("NewCell", "ENCSR123456", "total")
```

## Troubleshooting

If you encounter issues with file paths, create a link to the configuration file:
```bash
mkdir -p src/config
ln -s /path/to/rnaseq_analysis/config/settings.json src/config/settings.json
```

If visualizations don't appear, try the standalone visualization script:
```bash
python scripts/simple_viz.py
```