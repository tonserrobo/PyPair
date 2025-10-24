#!/bin/bash
#SBATCH --job-name=nearest-centroid
#SBATCH --output=nearest-centroid-%j.out
#SBATCH --error=nearest-centroid-%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --partition=k2-medpri
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

echo "Starting nearest centroid training job..."
module load apps/anaconda3/2024.10/bin

# activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pypair || { echo "Failed to activate Conda environment"; exit 1; }

# Run the script
python nearest_centroid.py
if [ $? -ne 0 ]; then
    echo "Nearest centroid encountered an error." >&2
    exit 1
fi

echo "Nearest Centroid training complete."
conda deactivate
