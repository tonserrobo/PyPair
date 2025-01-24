#!/bin/bash
#SBATCH --job-name=GaussianNaiveBayes-classifier
#SBATCH --output=GaussianNaiveBayes-classifier-%j.out
#SBATCH --error=GaussianNaiveBayes-classifier-%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --partition=k2-medpri
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

echo "Starting Gaussian Naive Bayes training job..."
module load apps/anaconda3/2024.10/bin

# activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pypair || { echo "Failed to activate Conda environment"; exit 1; }

# run job
python gaussianNaiveBayes_classifier.py
if [ $? -ne 0 ]; then
    echo "Gaussian Naive Bayes Classifier encountered an error." >&2
    exit 1
fi

echo "Gaussian Naive Bayes Classifier training complete"
conda deactivate