#!/bin/bash
#SBATCH --job-name=pypair_genome
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=logs/genome_training_%j.out
#SBATCH --error=logs/genome_training_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

# PyPair: Whole Human Genome Training on NI-HPC
# Trains winning models on full human genome (3 billion bases → ~3B k-mers)
# Production-scale validation for clinical deployment

echo "========================================================================"
echo "PyPair: Whole Human Genome Training"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 256GB"
echo "Start: $(date)"
echo "========================================================================"
echo "WARNING: This is a LARGE-SCALE job (~3 billion k-mers)"
echo "         Estimated runtime: 24-48 hours"
echo "         Recommend testing on chr21 first!"
echo "========================================================================"

# Environment setup
module purge
module load python/3.10
module load gcc/11.2.0

# Activate virtual environment
source venv/bin/activate

# Navigate to project root
cd $SLURM_SUBMIT_DIR

# Create output directories
mkdir -p results/whole_genome
mkdir -p results/shap/whole_genome
mkdir -p models/model_weights/production
mkdir -p logs

# Configuration
GENOME_FILE="data/reference_genomes/human_genome_grch38.hdf5"  # Update path
TRAIN_FRACTION=0.80
MAX_SAMPLES=3000000000  # ~3 billion k-mers

echo ""
echo "========================================================================"
echo "Dataset Configuration"
echo "========================================================================"
echo "Genome file: $GENOME_FILE"
echo "Train fraction: $TRAIN_FRACTION"
echo "Max samples: $MAX_SAMPLES"
echo "========================================================================"

# Check if genome file exists
if [ ! -f "$GENOME_FILE" ]; then
    echo "ERROR: Genome file not found: $GENOME_FILE"
    echo "Please run k-mer extraction first:"
    echo "  sbatch scripts/hpc/extract_genome_kmers.sh"
    exit 1
fi

# Run training (only on winning models for production)
echo ""
echo "========================================================================"
echo "Training Production Models"
echo "========================================================================"
echo "Models: Nearest Centroid, KNN (based on S1-S4 + chr21 validation)"
echo "Batch sizes: 100k, 500k, 1M, 5M"
echo "Non-inferiority testing: vs BWA (0.003 sec/seed)"
echo "========================================================================"

# For whole genome, we'll only train the top 2 models (NC + KNN)
# GNB may struggle with 3B samples

python scripts/train_whole_genome_winners.py \
    --genome-file "$GENOME_FILE" \
    --train-fraction $TRAIN_FRACTION \
    --batch-sizes 100000,500000,1000000,5000000 \
    --output-dir results/whole_genome \
    --save-models

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Job Complete"
echo "========================================================================"
echo "Exit code: $EXIT_CODE"
echo "End: $(date)"
echo "Runtime: $(($SECONDS / 3600))h $(($SECONDS % 3600 / 60))m"
echo "========================================================================"

# Deactivate environment
deactivate

exit $EXIT_CODE
