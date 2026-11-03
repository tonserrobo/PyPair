#!/bin/bash
#SBATCH --job-name=pypair_validation
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/pypair_validation_%j.out
#SBATCH --error=logs/pypair_validation_%j.err

################################################################################
# PyPair: S1-S4 Validation Only
################################################################################
#
# Runs comprehensive model validation across S1-S4 baseline datasets
# WITHOUT scaling analysis (faster execution).
#
# Reference Genome: GRCh38 (Human Genome Assembly)
#
# Usage:
#   sbatch scripts/hpc/pypair_validation_only.sh
#
################################################################################

set -e
mkdir -p logs

echo "=============================================================================="
echo "PyPair: S1-S4 Validation (GRCh38 Reference)"
echo "=============================================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=============================================================================="

cd $SLURM_SUBMIT_DIR

# Load modules
module purge
module load python/3.10

# Run validation
echo "Running comprehensive model analysis (S1-S4)..."
python scripts/comprehensive_model_analysis.py

echo "Running performance profiling..."
python scripts/profile_seed_mapping_fast.py

echo "Generating visualizations..."
python scripts/create_viva_plots.py

echo "=============================================================================="
echo "Validation Complete"
echo "=============================================================================="
echo "End time: $(date)"
echo "Results:"
echo "  - results/comprehensive_model_analysis.csv"
echo "  - results/seed_mapping_profile_fast.csv"
echo "  - docs/figures/*.png"
echo "=============================================================================="
