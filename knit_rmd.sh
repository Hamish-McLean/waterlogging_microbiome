#!/bin/bash
# filepath: /home/hmclean/waterlogging_microbiome/analysis/knit_rmd.sh

# Conda environment
CONDA_ENV="waterlogging_microbiome"

# Default pipeline name
DEFAULT_PIPELINE="DADA2"

# Check if a pipeline name is provided, otherwise use the default
PIPELINE=${1:-$DEFAULT_PIPELINE}

# Activate the Conda environment
echo "Activating Conda environment: $CONDA_ENV"
source "$(conda info --base)/etc/profile.d/conda.sh" || { echo "Conda not found"; exit 1; }
conda activate "$CONDA_ENV" || { echo "Failed to activate Conda environment: $CONDA_ENV"; exit 1; }

# Knit the R Markdown file with the specified pipeline parameter
echo "Knitting R Markdown file with pipeline: $PIPELINE"
Rscript -e "rmarkdown::render(\
  'microbiome_analysis.rmd', \
  params = list(pipeline = '$PIPELINE'), \
  output_file = paste0('microbiome_analysis_', '$PIPELINE', '.html')\
)" || { echo "Failed to knit R Markdown file"; exit 1; }

# Deactivate the Conda environment
echo "Deactivating Conda environment: $CONDA_ENV"
conda deactivate
