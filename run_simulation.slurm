#!/bin/bash
#SBATCH --account=def-username  # Replace with your account
#SBATCH --time=24:00:00        # Adjust based on your needs
#SBATCH --cpus-per-task=32     # Request 32 CPU cores
#SBATCH --gpus-per-node=1      # Request 1 GPU
#SBATCH --mem=64G              # Request 64GB RAM
#SBATCH --output=%x-%j.out     # Standard output
#SBATCH --error=%x-%j.err      # Standard error

# Load required modules
module load python/3.9
module load cuda/11.4
module load scipy-stack

# Create and activate virtual environment
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

# Install required packages
pip install --no-index --upgrade pip
pip install numpy pandas matplotlib seaborn tqdm torch networkx

# Run the simulation
python main.py 