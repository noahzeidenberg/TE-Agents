#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out

# Load required modules
module purge
module load StdEnv/2020
module load gcc/9.3.0
module load python/3.11
module load blast+/2.14.0
module load hmmer/3.3.2
module load trf/4.09.1
module load perl/5.30.2
module load muscle
module load cd-hit
module load gt

# Source the tools configuration
source /scratch/$USER/tools/tools_config.sh

# Ensure directories exist
mkdir -p edta_output

# Run the analysis
python yeast_te_analyzer.py \
    --genome r64/GCF_000146045.2_R64_genomic.fasta \
    --gff r64/GCF_000146045.2_R64_genomic.gff \
    --output yeast_te_results

# Clean up
rm -f edta.sif 