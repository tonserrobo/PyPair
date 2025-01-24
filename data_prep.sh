#!/bin/bash
#SBATCH --job-name=data-prep
#SBATCH --output=data-prep-%j.out
#SBATCH --error=data-prep-%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --partition=k2-medpri
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

echo "Starting data preparation job..."

# 1) Load the necessary modules and initialize Conda
module load apps/anaconda3/2024.10/bin

# Set Conda environment directories (ensures installation on scratch space)
mkdir -p /mnt/scratch2/users/$USER/conda
export CONDA_PKGS_DIRS=/mnt/scratch2/users/$USER/conda/pkgs
export CONDA_ENVS_PATH=/mnt/scratch2/users/$USER/conda/envs

# Initialize Conda properly to avoid the 'conda init' issue
source $(conda info --base)/etc/profile.d/conda.sh

# 2) Create/activate environment (if not already existing)
if conda env list | grep -q "^pypair "; then
    echo "Environment 'pypair' found."
    conda activate pypair
else
    echo "Environment 'pypair' not found, creating..."
    conda create -n pypair python=3.8 -y
    conda activate pypair
    echo "Environment 'pypair' created."
fi

# 3) Install required packages in the 'pypair' environment from requirements.txt
if [ -f requirements.txt ]; then
    echo "Installing dependencies from requirements.txt..."
    pip install --no-cache-dir -r requirements.txt
else
    echo "ERROR: requirements.txt not found in current directory."
    exit 1
fi

# 4) Run the data preparation script
echo "Running multimodel_training.py..."
python multimodel_training.py

if [ $? -eq 0 ]; then
    echo "Data preparation completed successfully."
else
    echo "ERROR: Data preparation script failed."
    exit 1
fi
