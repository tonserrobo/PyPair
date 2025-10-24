#!/bin/bash
#SBATCH --job-name=random-forest-classifier
#SBATCH --output=random-forest-classifier-%j.out
#SBATCH --error=random-forest-classifier-%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --partition=k2-medpri
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

echo "Starting Random Forest training job..."
module load apps/anaconda3/2024.10/bin

# activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pypair || { echo "Failed to activate Conda environment"; exit 1; }

# run job
python random_foreset.py
if [ $? -ne 0 ]; then
    echo "Random Forest encountered an error." >&2
    exit 1
fi

echo "Random Forest training complete"
conda deactivate