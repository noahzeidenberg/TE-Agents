#!/bin/bash
#SBATCH --account=def-username  # Replace with your account
#SBATCH --time=24:00:00        # Adjust based on your needs
#SBATCH --cpus-per-task=32     # Request 32 CPU cores
#SBATCH --gpus-per-node=1      # Request 1 GPU
#SBATCH --mem=64G              # Request 64GB RAM
#SBATCH --output=logs/%x-%j.out     # Standard output
#SBATCH --error=%x-%j.err      # Standard error

# Load required modules
module load python
module load cuda
module load scipy-stack

# Set up virtual environment with system packages
virtualenv --system-site-packages $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

# Configure pip to use the Compute Canada wheels
WHEELHOUSE=/cvmfs/soft.computecanada.ca/custom/python/wheelhouse/gentoo/avx2
WHEELHOUSE_NOAVX=/cvmfs/soft.computecanada.ca/custom/python/wheelhouse/gentoo/generic
export PIP_CONFIG_FILE=/cvmfs/soft.computecanada.ca/config/python/pip.conf
export PIP_FIND_LINKS="$WHEELHOUSE $WHEELHOUSE_NOAVX"
export PIP_NO_INDEX=1

# Install required packages
pip install --no-index --upgrade pip
pip install numpy pandas matplotlib seaborn tqdm torch networkx

# Run the simulation
echo "Starting simulation..."
python main.py 