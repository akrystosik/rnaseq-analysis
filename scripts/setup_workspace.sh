#!/bin/bash

# Setup script for RNA-seq analysis workspace in Kubernetes environment
# Creates directory structure and installs required packages

# Set base directory
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
SCRIPT_DIR="${BASE_DIR}/scripts"

echo "Setting up RNA-seq analysis workspace..."

# Determine Python path
PYTHON_PATH=$(which python)
echo "Using Python at: ${PYTHON_PATH}"
$PYTHON_PATH --version

# Install pip if needed
echo "Ensuring pip is installed..."
if ! $PYTHON_PATH -m pip --version > /dev/null 2>&1; then
    echo "Installing pip..."
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
    $PYTHON_PATH get-pip.py
    rm get-pip.py
fi

# Install required system dependencies
echo "Installing system dependencies..."
apt-get update && apt-get install -y \
    build-essential \
    gcc \
    python3-dev \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev

# Install required Python packages using the correct Python interpreter
echo "Installing Python packages..."
$PYTHON_PATH -m pip install --no-cache-dir \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    requests \
    scipy \
    pyBigWig \
    scanpy \
    h5py \
    tqdm \
    tables \
    statsmodels \
    openpyxl

# Install anndata separately to ensure it completes successfully
echo "Installing anndata package..."
$PYTHON_PATH -m pip install --no-cache-dir anndata

# Verify installations
echo "Verifying installations..."
$PYTHON_PATH -c "import numpy; print(f'NumPy {numpy.__version__} installed successfully')"
$PYTHON_PATH -c "import pandas; print(f'Pandas {pandas.__version__} installed successfully')"
$PYTHON_PATH -c "import matplotlib; print(f'Matplotlib {matplotlib.__version__} installed successfully')"
$PYTHON_PATH -c "import seaborn; print(f'Seaborn {seaborn.__version__} installed successfully')"
$PYTHON_PATH -c "import requests; print(f'Requests {requests.__version__} installed successfully')"
$PYTHON_PATH -c "import scipy; print(f'SciPy {scipy.__version__} installed successfully')"
$PYTHON_PATH -c "import pyBigWig; print('pyBigWig installed successfully')"
$PYTHON_PATH -c "import h5py; print(f'h5py {h5py.__version__} installed successfully')"
$PYTHON_PATH -c "import tqdm; print(f'tqdm {tqdm.__version__} installed successfully')"
$PYTHON_PATH -c "import openpyxl; print(f'openpyxl {openpyxl.__version__} installed successfully')"

# Verify anndata installation specifically
echo "Verifying anndata installation..."
$PYTHON_PATH -c "
try:
    import anndata
    print(f'anndata {anndata.__version__} installed successfully')
except ImportError as e:
    print(f'Warning: anndata package not installed properly. Error: {e}')
    exit(1)
"

echo "Setup complete!"