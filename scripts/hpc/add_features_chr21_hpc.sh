#!/bin/bash
#SBATCH --job-name=pypair_chr21_features
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=logs/chr21_features_%j.out
#SBATCH --error=logs/chr21_features_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

# PyPair: chr21 Extended Feature Engineering on NI-HPC
# Adds S1-S4 compatible features to chr21 k-mers (80M rows)
# Required before training models on chr21

echo "========================================================================"
echo "PyPair: chr21 Extended Feature Engineering"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 32GB"
echo "Start: $(date)"
echo "========================================================================"
echo "Processing: 80,175,810 k-mers"
echo "Features to add:"
echo "  - Nucleotide counts (first, A, C, G, T)"
echo "  - Reverse complement encoding (rc_seeds)"
echo "  - Run-length encoding (runlength_*)"
echo "========================================================================"

# Environment setup
module purge
module load python/3.10
module load gcc/11.2.0

# Activate virtual environment
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

source venv/bin/activate

# Install dependencies
echo "Installing dependencies..."
pip install --upgrade pip
pip install pandas numpy
pip install vaex vaex-hdf5

# Navigate to project root
cd $SLURM_SUBMIT_DIR

# Create output directories
mkdir -p data/reference_genomes
mkdir -p logs

# Check input file exists
INPUT_FILE="data/reference_genomes/chr21_kmers.hdf5"
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: chr21 k-mer file not found: $INPUT_FILE"
    echo "Please run k-mer extraction first:"
    echo "  sbatch scripts/hpc/extract_chromosome_data.sh"
    exit 1
fi

# Run feature engineering
echo ""
echo "========================================================================"
echo "Starting Feature Engineering"
echo "========================================================================"

python scripts/add_extended_features_chr21.py

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Job Complete"
echo "========================================================================"
echo "Exit code: $EXIT_CODE"
echo "End: $(date)"
echo "Runtime: $(($SECONDS / 60))m $(($SECONDS % 60))s"

if [ $EXIT_CODE -eq 0 ]; then
    OUTPUT_FILE="data/reference_genomes/chr21_kmers_extended.hdf5"
    if [ -f "$OUTPUT_FILE" ]; then
        SIZE_GB=$(du -h "$OUTPUT_FILE" | cut -f1)
        echo "Output file: $OUTPUT_FILE ($SIZE_GB)"
        echo ""
        echo "Next step: Train models on chr21"
        echo "  sbatch scripts/hpc/train_chr21_hpc.sh"
    fi
fi

echo "========================================================================"

# Deactivate environment
deactivate

exit $EXIT_CODE
