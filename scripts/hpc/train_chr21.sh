#!/bin/bash
#SBATCH --job-name=pypair_gnb_chr21
#SBATCH --partition=k2-medpri
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/chr21_train_%j.out
#SBATCH --error=logs/chr21_train_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

#================================================
# PyPair: Train GNB on chr21
# Expected: ~2 minutes, F1 = 95.2%
#================================================

echo "=========================================="
echo "PyPair: Chr21 GNB Training"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=========================================="

# Load environment
module load apps/anaconda3/2024.10/bin
conda activate pypair

# Navigate to project
cd /scratch/k1234567/PyPair

# Verify data file exists
DATA_FILE="data/reference_genomes/chr21_kmers.hdf5"
if [ ! -f "$DATA_FILE" ]; then
    echo "ERROR: Data file not found: $DATA_FILE"
    echo "Please run extract_chromosome_data.py first"
    exit 1
fi

echo ""
echo "Data file found: $DATA_FILE"
echo "Size: $(du -h $DATA_FILE | cut -f1)"
echo ""

# Train model
echo "=========================================="
echo "Training GNB on chr21..."
echo "=========================================="

python scripts/train_on_chromosomes.py \
    --data $DATA_FILE \
    --output models/model_weights/gnb_chr21.joblib \
    --test-size 0.2 \
    --random-state 42

TRAIN_STATUS=$?

if [ $TRAIN_STATUS -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Training completed successfully!"
    echo "=========================================="

    # Run SHAP analysis
    echo ""
    echo "=========================================="
    echo "Running SHAP explainability analysis..."
    echo "=========================================="

    python scripts/utils/shap_analysis.py \
        --model models/model_weights/gnb_chr21.joblib \
        --data $DATA_FILE \
        --output results/shap/gnb_chr21 \
        --samples 1000 \
        --model-type kernel

    SHAP_STATUS=$?

    if [ $SHAP_STATUS -eq 0 ]; then
        echo ""
        echo "SHAP analysis complete!"
        echo "Results: results/shap/gnb_chr21/"
    else
        echo "WARNING: SHAP analysis failed (non-fatal)"
    fi

    # Print summary
    echo ""
    echo "=========================================="
    echo "SUMMARY"
    echo "=========================================="
    echo "Model: models/model_weights/gnb_chr21.joblib"
    echo "Size: $(du -h models/model_weights/gnb_chr21.joblib | cut -f1)"
    echo "Metrics: models/model_weights/gnb_chr21_metrics.json"
    echo ""
    cat models/model_weights/gnb_chr21_metrics.json
    echo ""
    echo "=========================================="
else
    echo ""
    echo "ERROR: Training failed with exit code $TRAIN_STATUS"
    exit $TRAIN_STATUS
fi

echo "End time: $(date)"
echo "=========================================="
echo "Job complete!"
