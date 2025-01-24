#!/bin/bash
#SBATCH --job-name=kneighbours-classifier
#SBATCH --output=kneighbours-classifier-%j.out
#SBATCH --error=kneighbours-classifier-%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --partition=k2-medpri
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

echo "Starting KNeighbours Classifier training job..."
module load apps/anaconda3/2024.10/bin

# activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pypair || { echo "Failed to activate Conda environment"; exit 1; }

# run job
python kNeighbours_classifier.py
if [ $? -ne 0 ]; then
    echo "KNeighbours Classifier encountered an error." >&2
    exit 1
fi

echo "kNeighbours Classifier training complete"
conda deactivate