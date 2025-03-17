#!/bin/bash
#SBATCH --account=def-yourpi
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out

# Load required modules
module load apptainer

# Pull the EDTA container in advance
apptainer pull docker://oushujun/edta:2.0.0

# Ensure directories exist
mkdir -p edta_output

# Run the analysis
python yeast_te_analyzer.py \
    --genome r64/GCF_000146045.2_R64_genomic.fasta \
    --gff r64/GCF_000146045.2_R64_genomic.gff \
    --output yeast_te_results

# Clean up temporary files
rm -f edta_2.0.0.sif 