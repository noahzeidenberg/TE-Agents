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

# Activate virtual environment if you have one
# source ~/env/bin/activate  # Uncomment if needed

# Set up working directory
TEST_DIR="$SLURM_TMPDIR/te_test"
mkdir -p $TEST_DIR
cd $TEST_DIR

# Extract chromosome III using Python
cat > extract_chr3.py << 'EOF'
from Bio import SeqIO
from pathlib import Path

def extract_chromosome(input_file, output_file, chr_id="NC_001135.5"):
    """Extract chromosome III from the genome file"""
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for record in SeqIO.parse(fin, 'fasta'):
            if chr_id in record.id:
                SeqIO.write(record, fout, 'fasta')
                print(f"Extracted {record.id}, length: {len(record.seq)} bp")
                return True
    return False

if __name__ == "__main__":
    genome_file = "GCF_000146045.2_R64_genomic.fasta"
    output_file = "chr3.fasta"
    
    if extract_chromosome(genome_file, output_file):
        print(f"Successfully extracted chromosome III to {output_file}")
    else:
        print("Failed to find chromosome III")
EOF

# Copy genome file to temp directory
cp /scratch/nzeidenb/TE-Agents/r64/GCF_000146045.2_R64_genomic.fasta .

# Run extraction script
python extract_chr3.py

# Use RepeatMasker for initial TE identification
# Note: Using basic parameters for test run
RepeatMasker -species Saccharomyces -pa $SLURM_CPUS_PER_TASK chr3.fasta

# Copy results back
mkdir -p /scratch/nzeidenb/TE-Agents/te_test_results
cp chr3.fasta.out /scratch/nzeidenb/TE-Agents/te_test_results/
cp chr3.fasta.tbl /scratch/nzeidenb/TE-Agents/te_test_results/

echo "Test analysis complete. Results are in /scratch/nzeidenb/TE-Agents/te_test_results/" 