#!/bin/bash
#SBATCH --job-name=XGB-classifier
#SBATCH --output=XGB-classifier-%j.out
#SBATCH --error=XGB-classifier-%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --partition=k2-medpri
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

echo "Starting XGB Classifier training job..."
module load apps/anaconda3/2024.10/bin

# activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pypair || { echo "Failed to activate Conda environment"; exit 1; }

# run job
python XGB_classifier.py
if [ $? -ne 0 ]; then
    echo "XGB Classifier encountered an error." >&2
    exit 1
fi

echo "XGB Classifier training complete"
conda deactivate