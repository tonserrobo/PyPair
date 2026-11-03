#!/bin/bash
#SBATCH --job-name=pypair_scaling
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/pypair_scaling_%j.out
#SBATCH --error=logs/pypair_scaling_%j.err

################################################################################
# PyPair: Scaling Analysis Only (KNN vs GNB)
################################################################################
#
# Tests hypothesis: GNB will outperform KNN at genome scale due to lower
# memory requirements.
#
# Reference Genome: GRCh38 (Human Genome Assembly)
#
# Data Scales Tested:
#   - S1 (44.8k samples)
#   - S1+S2 (~90k samples)
#   - S1+S2+S3 (~135k samples)
#   - S1+S2+S3+S4 (~180k samples)
#
# Usage:
#   sbatch scripts/hpc/pypair_scaling_only.sh
#
################################################################################

set -e
mkdir -p logs

echo "=============================================================================="
echo "PyPair: Scaling Analysis (KNN vs GNB Memory Comparison)"
echo "=============================================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Reference: GRCh38"
echo "Hypothesis: GNB scales better than KNN for genome-scale data"
echo "=============================================================================="

cd $SLURM_SUBMIT_DIR

# Load modules
module purge
module load python/3.10

# Run scaling analysis
echo "Running scaling analysis..."
python scripts/scaling_analysis.py

# Generate scaling plots
echo "Generating scaling visualizations..."
python scripts/create_scaling_plots.py

echo "=============================================================================="
echo "Scaling Analysis Complete"
echo "=============================================================================="
echo "Results:"
echo "  - results/scaling_analysis.csv"
echo "  - docs/figures/scaling_*.png"
echo "=============================================================================="
