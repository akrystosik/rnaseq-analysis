#!/bin/bash
# Script to prepare the repository for GitHub

# Set up Git if not already initialized
if [ ! -d "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/.git" ]; then
    cd /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq
    git init
    echo "Git repository initialized"
else
    echo "Git repository already exists"
fi

# Create .gitignore file
cat > /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/.gitignore << 'GITIGNORE'
# Large data files
*.h5ad
*.gct.gz
*.tsv
*.csv.gz

# Log files
logs/*

# Temporary directories
/temp/
__pycache__/
*.py[cod]
*$py.class

# Jupyter notebook checkpoints
.ipynb_checkpoints/

# OS specific files
.DS_Store
Thumbs.db
GITIGNORE

# Stage files for commit
cd /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq
git add scripts/
git add metadata/json/*.json
git add README.md CHANGELOG.md USAGE.md .gitignore

# Initial commit
git commit -m "Initial commit: RNA-seq standardization pipeline"

# Set up remote if not already configured
git remote -v | grep -q "origin" || git remote add origin https://github.com/akrystosik/rnaseq-standardization-pipeline.git

echo "Repository is ready to push to GitHub"
echo "Run the following commands to push:"
echo "  cd /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq"
echo "  git branch -M main"
echo "  git push -u origin main"
