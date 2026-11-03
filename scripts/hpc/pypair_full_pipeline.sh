#!/bin/bash
#SBATCH --job-name=pypair_full_pipeline
#SBATCH --partition=k2-hipri
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pypair_pipeline_%j.out
#SBATCH --error=logs/pypair_pipeline_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@ulster.ac.uk

################################################################################
# PyPair: Full HPC Training/Analysis Pipeline
################################################################################
#
# This master script orchestrates the complete PyPair pipeline:
#   1. Data preparation (feature engineering)
#   2. Model training (KNN, GNB, Nearest Centroid)
#   3. Model validation (S1-S4)
#   4. Scaling analysis (memory profiling)
#   5. Performance benchmarking
#
# Usage:
#   sbatch scripts/hpc/pypair_full_pipeline.sh
#
# Output:
#   results/pipeline_complete.txt - Summary report
#   results/scaling_analysis.csv - Memory scaling results
#
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Create logs directory
mkdir -p logs

# Load modules (adjust for your HPC environment)
module purge
module load python/3.10
module load gcc/11.2.0

# Activate virtual environment (if using)
# source ~/venvs/pypair/bin/activate

echo "=============================================================================="
echo "PyPair: Full HPC Pipeline Started"
echo "=============================================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Start time: $(date)"
echo "=============================================================================="

# Change to project directory
cd $SLURM_SUBMIT_DIR
echo "Working directory: $(pwd)"

################################################################################
# STAGE 1: Data Preparation (if needed)
################################################################################
echo ""
echo "STAGE 1: Data Preparation"
echo "------------------------------------------------------------------------------"

if [ ! -f "data/reference_samples/cleaned/S1_extended_encoded_X_train.csv" ]; then
    echo "Running data preparation..."
    python scripts/training_data_preparation.py
    echo "[OK] Data preparation complete"
else
    echo "[SKIP] Training data already exists"
fi

################################################################################
# STAGE 2: Model Training & Validation (S1-S4)
################################################################################
echo ""
echo "STAGE 2: Model Training & Validation (S1-S4)"
echo "------------------------------------------------------------------------------"

echo "Training models across S1-S4 datasets..."
python scripts/comprehensive_model_analysis.py

echo "[OK] S1-S4 validation complete"
echo "Results: results/comprehensive_model_analysis.csv"

################################################################################
# STAGE 3: Seed Mapping Performance Profiling
################################################################################
echo ""
echo "STAGE 3: Seed Mapping Performance Profiling"
echo "------------------------------------------------------------------------------"

echo "Profiling seed mapping performance (sec/seed, throughput)..."
python scripts/profile_seed_mapping_fast.py

echo "[OK] Performance profiling complete"
echo "Results: results/seed_mapping_profile_fast.csv"

################################################################################
# STAGE 4: Scaling Analysis (KNN vs GNB Memory Comparison)
################################################################################
echo ""
echo "STAGE 4: Scaling Analysis - Memory Profiling"
echo "------------------------------------------------------------------------------"

echo "Testing hypothesis: GNB scales better than KNN at large data sizes"
python scripts/scaling_analysis.py

echo "[OK] Scaling analysis complete"
echo "Results: results/scaling_analysis.csv"

################################################################################
# STAGE 5: Generate Visualizations
################################################################################
echo ""
echo "STAGE 5: Generate Viva Visualizations"
echo "------------------------------------------------------------------------------"

echo "Creating publication-quality plots..."
python scripts/create_viva_plots.py

# Create scaling plots
python scripts/create_scaling_plots.py

echo "[OK] Visualizations complete"
echo "Results: docs/figures/*.png"

################################################################################
# STAGE 6: Summary Report
################################################################################
echo ""
echo "STAGE 6: Generating Summary Report"
echo "------------------------------------------------------------------------------"

python scripts/generate_pipeline_summary.py

echo "[OK] Summary report generated: results/pipeline_complete.txt"

################################################################################
# Pipeline Complete
################################################################################
echo ""
echo "=============================================================================="
echo "PyPair: Full HPC Pipeline Complete"
echo "=============================================================================="
echo "End time: $(date)"
echo "Total duration: $SECONDS seconds"
echo ""
echo "Output Files:"
echo "  - results/comprehensive_model_analysis.csv"
echo "  - results/seed_mapping_profile_fast.csv"
echo "  - results/scaling_analysis.csv"
echo "  - results/pipeline_complete.txt"
echo "  - docs/figures/*.png"
echo "=============================================================================="

# Send completion notification (optional)
# mail -s "PyPair Pipeline Complete (Job $SLURM_JOB_ID)" your.email@ulster.ac.uk < results/pipeline_complete.txt

exit 0
