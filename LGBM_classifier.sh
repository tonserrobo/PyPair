#!/bin/bash
#SBATCH --job-name=LGBM-classifier
#SBATCH --output=LGBM-classifier-%j.out
#SBATCH --error=LGBM-classifier-%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --partition=k2-medpri
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

echo "Starting LGBM Classifier training job..."
module load apps/anaconda3/2024.10/bin

# activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pypair || { echo "Failed to activate Conda environment"; exit 1; }

# run job
python LGBM_classifier.py
if [ $? -ne 0 ]; then
    echo "LGBM Classifier encountered an error." >&2
    exit 1
fi

echo "LGBM Classifier training complete"