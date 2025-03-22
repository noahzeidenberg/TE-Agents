#!/bin/bash
#SBATCH --account=def-yourpi
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=%x-%j.out

# Load required modules
module purge
module load StdEnv/2020
module load gcc/9.3.0
module load python/3.11
module load blast+/2.14.0
module load hmmer/3.3.2
module load perl/5.30.2

# Set up working directory
TEST_DIR="$SLURM_TMPDIR/te_test"
mkdir -p $TEST_DIR
cd $TEST_DIR

# First verify the genome file exists
if [ ! -f "/scratch/nzeidenb/TE-Agents/r64/GCF_000146045.2_R64_genomic.fasta" ]; then
    echo "Error: Genome file not found at expected location"
    exit 1
fi

# Run extraction script
python extract_chr3.py

# Verify the extraction worked
if [ ! -f "chr3.fasta" ]; then
    echo "Error: Chromosome extraction failed"
    exit 1
fi

# Create results directory if it doesn't exist
mkdir -p /scratch/nzeidenb/TE-Agents/te_test_results

# Just count sequences and basic stats for now
echo "Basic sequence statistics:" > chr3_stats.txt
grep -c ">" chr3.fasta >> chr3_stats.txt
wc -l chr3.fasta >> chr3_stats.txt

# Copy results back
cp chr3.fasta chr3_stats.txt /scratch/nzeidenb/TE-Agents/te_test_results/

echo "Test analysis complete. Results are in /scratch/nzeidenb/TE-Agents/te_test_results/" 