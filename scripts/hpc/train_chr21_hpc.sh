#!/bin/bash
#SBATCH --job-name=pypair_chr21
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/chr21_training_%j.out
#SBATCH --error=logs/chr21_training_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

# PyPair: chr21 Genome-Scale Training on NI-HPC
# Trains winning models (Nearest Centroid, KNN, GNB) on 80M k-mers
# with batch profiling (50k-500k) and non-inferiority testing

echo "========================================================================"
echo "PyPair: chr21 Genome-Scale Training"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64GB"
echo "Start: $(date)"
echo "========================================================================"

# Environment setup
module purge
module load python/3.10
module load gcc/11.2.0

# Activate virtual environment (or create if needed)
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

source venv/bin/activate

# Install dependencies
echo "Installing dependencies..."
pip install --upgrade pip
pip install pandas numpy scipy scikit-learn joblib
pip install vaex vaex-hdf5
pip install memory-profiler psutil
pip install researchpy statsmodels
pip install shap matplotlib seaborn

# Navigate to project root
cd $SLURM_SUBMIT_DIR

# Create output directories
mkdir -p results/chr21
mkdir -p results/shap/chr21
mkdir -p models/model_weights
mkdir -p logs

# Run training
echo ""
echo "========================================================================"
echo "Starting chr21 Training Pipeline"
echo "========================================================================"
echo "Dataset: data/reference_genomes/chr21_kmers.hdf5 (80M k-mers)"
echo "Models: Nearest Centroid, KNN, Gaussian NB"
echo "Batch sizes: 50k, 100k, 200k, 500k"
echo "========================================================================"

python scripts/train_chr21_winners.py

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Job Complete"
echo "========================================================================"
echo "Exit code: $EXIT_CODE"
echo "End: $(date)"
echo "========================================================================"

# Deactivate environment
deactivate

exit $EXIT_CODE
