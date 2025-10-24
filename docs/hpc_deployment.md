# HPC Deployment Guide: NIHPC Cluster

## Overview

PyPair is designed for deployment on the **Northern Ireland High Performance Computing (NIHPC)** cluster. This guide covers environment setup, job submission, and best practices for training models at scale.

## NIHPC Cluster Specifications

### Available Resources

- **Partition**: `k2-medpri` (medium priority)
- **Compute Nodes**: 100+ nodes with varying specs
- **Memory per Node**: Up to 512GB RAM
- **CPUs per Node**: 40-80 cores (Intel Xeon)
- **GPUs**: NVIDIA A100, V100 available on `gpu` partition
- **Storage**:
  - Home: `/users/k1234567/` (100GB quota)
  - Scratch: `/scratch/k1234567/` (10TB, purged after 30 days)
- **Module System**: Environment Modules (Lmod)
- **Scheduler**: SLURM Workload Manager

### Network & File System

- **Shared File System**: Lustre (parallel I/O optimized)
- **Home Directory**: NFS-mounted, slow for large I/O
- **Scratch**: High-performance parallel file system (use for training data)

## Environment Setup

### 1. Initial Login

```bash
# SSH to NIHPC login node
ssh k1234567@login.hpc.ac.uk

# Check quota
quota -s

# Navigate to scratch (for large datasets)
cd /scratch/k1234567/
mkdir -p PyPair
cd PyPair
```

### 2. Clone Repository

```bash
# Clone from GitHub (or transfer via scp)
git clone https://github.com/your-org/PyPair.git
cd PyPair

# Alternatively, transfer from local machine
# scp -r PyPair/ k1234567@login.hpc.ac.uk:/scratch/k1234567/
```

### 3. Load Modules

```bash
# Load Anaconda module
module load apps/anaconda3/2024.10/bin

# Verify Python version
python --version  # Should be 3.10.x or 3.11.x

# List available modules
module avail
```

### 4. Create Conda Environment

```bash
# Create environment with Python 3.10 (Vaex compatibility)
conda create -n pypair python=3.10 -y
conda activate pypair

# Verify activation
which python
# Output: /opt/gridware/apps/anaconda3/2024.10/envs/pypair/bin/python
```

### 5. Install PyPair

#### Option A: Using UV (Recommended)

```bash
# Install UV
pip install uv

# Install PyPair with all dependencies
uv pip install -e ".[hpc]"

# Verify installation
python -c "import vaex; print(vaex.__version__)"
python -c "import polars; print(polars.__version__)"
```

#### Option B: Using pip

```bash
# Install from requirements.txt
pip install -r requirements.txt

# Install PyPair in development mode
pip install -e .
```

### 6. Verify Installation

```bash
# Test imports
python -c "
from accelerated_GNB.GNB_vaex.GNB_HDF5_optimised import VaexGaussianNB
from accelerated_GNB.GNB_polars.GNB_polars import GaussianNBPolars
print('Installation successful!')
"
```

## Data Preparation

### 1. Transfer Data to Scratch

```bash
# From local machine, transfer large datasets
scp -r data/reference_genomes/ k1234567@login.hpc.ac.uk:/scratch/k1234567/PyPair/data/

# Or use rsync for resumable transfers
rsync -avz --progress data/reference_genomes/ \
    k1234567@login.hpc.ac.uk:/scratch/k1234567/PyPair/data/reference_genomes/
```

### 2. Convert Data to HDF5 (for Vaex)

Create a data preparation job:

```bash
# data_prep.sh
#!/bin/bash
#SBATCH --job-name=pypair_data_prep
#SBATCH --partition=k2-medpri
#SBATCH --mem=40GB
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/data_prep_%j.out
#SBATCH --error=logs/data_prep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

module load apps/anaconda3/2024.10/bin
conda activate pypair

cd /scratch/k1234567/PyPair

# Convert CSV to HDF5 for Vaex
python -c "
import vaex
import glob

for csv_file in glob.glob('data/reference_samples/*.csv'):
    print(f'Converting {csv_file}...')
    df = vaex.from_csv(csv_file, convert=True, chunk_size=1_000_000)
    hdf5_file = csv_file.replace('.csv', '.hdf5')
    df.export_hdf5(hdf5_file)
    print(f'Saved to {hdf5_file}')
"

echo "Data preparation complete!"
```

Submit job:
```bash
mkdir -p logs
sbatch data_prep.sh
```

## Training Jobs

### 1. Gaussian Naive Bayes (Vaex) - Primary Model

**gaussianNaiveBayes_classifier.sh**:

```bash
#!/bin/bash
#SBATCH --job-name=pypair_gnb_vaex
#SBATCH --partition=k2-medpri
#SBATCH --mem=8GB  # Vaex uses minimal RAM
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/gnb_vaex_%j.out
#SBATCH --error=logs/gnb_vaex_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

module load apps/anaconda3/2024.10/bin
conda activate pypair

cd /scratch/k1234567/PyPair

# Run training
python -u << 'PYTHON_SCRIPT'
import sys
import vaex
import numpy as np
from sklearn.metrics import f1_score, precision_score, recall_score
import joblib
from accelerated_GNB.GNB_vaex.GNB_HDF5_optimised import VaexGaussianNB

print("Loading training data...")
train_df = vaex.open("data/reference_samples/S1_extended_encoded_x_train.hdf5")
y_train = vaex.open("data/reference_samples/S1_extended_encoded_y_train.hdf5")

print(f"Training samples: {len(train_df)}")
print(f"Features: {len(train_df.columns)}")

# Train model
print("Training Gaussian Naive Bayes...")
model = VaexGaussianNB()
model.fit(train_df, y_train["label"])

print("Model training complete!")

# Evaluate on test set
print("Loading test data...")
test_df = vaex.open("data/reference_samples/S1_extended_encoded_X_test.hdf5")
y_test = vaex.open("data/reference_samples/S1_extended_encoded_y_test.hdf5")["label"]

print("Predicting (chunked)...")
y_pred = model.predict_chunked(test_df, chunk_size=10000)

# Compute metrics
f1 = f1_score(y_test.values, y_pred, average='weighted')
precision = precision_score(y_test.values, y_pred, average='weighted')
recall = recall_score(y_test.values, y_pred, average='weighted')

print(f"\n=== Results ===")
print(f"F1 Score: {f1:.4f}")
print(f"Precision: {precision:.4f}")
print(f"Recall: {recall:.4f}")

# Save model
output_path = "models/model_weights/gaussian_nb_vaex_model.joblib"
joblib.dump(model, output_path)
print(f"Model saved to {output_path}")

PYTHON_SCRIPT

echo "Job complete!"
```

### 2. Full Genome Training (Vaex)

**train_full_genome.sh**:

```bash
#!/bin/bash
#SBATCH --job-name=pypair_gnb_full_genome
#SBATCH --partition=k2-medpri
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/gnb_full_genome_%j.out
#SBATCH --error=logs/gnb_full_genome_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=t.robinson@ulster.ac.uk

module load apps/anaconda3/2024.10/bin
conda activate pypair

cd /scratch/k1234567/PyPair

python -u << 'PYTHON_SCRIPT'
from accelerated_GNB.GNB_vaex.GNB_HDF5_optimised import VaexGaussianNB
import joblib
import time

print("Training on full human reference genome (GRCh38)...")
start_time = time.time()

model = VaexGaussianNB()
model.fit("/scratch/k1234567/reference_genomes/GRCh38_full_kmers.hdf5")

elapsed = time.time() - start_time
print(f"Training completed in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")

# Save model
joblib.dump(model, "models/model_weights/gnb_full_genome_grch38.joblib")
print("Model saved!")

# Test prediction speed
print("Testing prediction speed...")
import vaex
test_df = vaex.open("data/test_reads_sample.hdf5")
start = time.time()
predictions = model.predict_chunked(test_df, chunk_size=10000)
pred_time = time.time() - start
print(f"Predicted {len(predictions)} reads in {pred_time:.2f}s")
print(f"Average: {pred_time/len(predictions)*1000:.4f} ms per read")

PYTHON_SCRIPT
```

### 3. GPU Training (CUDA)

**gnb_cuda_train.sh**:

```bash
#!/bin/bash
#SBATCH --job-name=pypair_cuda
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/gnb_cuda_%j.out
#SBATCH --error=logs/gnb_cuda_%j.err

module load cuda/11.8
module load apps/anaconda3/2024.10/bin
conda activate pypair

cd /scratch/k1234567/PyPair/accelerated_GNB/Cpp

# Compile CUDA code
nvcc -std=c++11 -O3 GNB_cuda.cu -o gnb_cuda

# Run GPU prediction on large test set
./gnb_cuda \
    --model ../../models/model_weights/gnb_weights.bin \
    --input ../../data/test_reads_1M.csv \
    --output ../../results/predictions_cuda.txt \
    --batch-size 100000

echo "GPU job complete!"
```

## Job Orchestration

### Master Script: run_all.sh

Orchestrate all training jobs with dependencies:

```bash
#!/bin/bash
# run_all.sh - Submit all PyPair training jobs with dependencies

echo "Submitting PyPair training pipeline..."

# 1. Data preparation (runs first)
jobid_data=$(sbatch --parsable data_prep.sh)
echo "Data prep job: $jobid_data"

# 2. GNB Vaex training (depends on data prep)
jobid_gnb=$(sbatch --parsable --dependency=afterok:$jobid_data \
                   gaussianNaiveBayes_classifier.sh)
echo "GNB Vaex job: $jobid_gnb"

# 3. Full genome training (depends on GNB completion)
jobid_full=$(sbatch --parsable --dependency=afterok:$jobid_gnb \
                    train_full_genome.sh)
echo "Full genome job: $jobid_full"

# 4. Comparative models (run in parallel after data prep)
jobid_rf=$(sbatch --parsable --dependency=afterok:$jobid_data random_forest.sh)
jobid_lgbm=$(sbatch --parsable --dependency=afterok:$jobid_data LGBM_classifier.sh)
jobid_xgb=$(sbatch --parsable --dependency=afterok:$jobid_data XGB_classifier.sh)

echo "Random Forest job: $jobid_rf"
echo "LGBM job: $jobid_lgbm"
echo "XGBoost job: $jobid_xgb"

# 5. Benchmarking (depends on all models)
jobid_bench=$(sbatch --parsable \
    --dependency=afterok:$jobid_gnb:$jobid_rf:$jobid_lgbm:$jobid_xgb \
    benchmark_models.sh)
echo "Benchmark job: $jobid_bench"

echo "All jobs submitted!"
echo "Monitor with: squeue -u $USER"
```

Submit pipeline:
```bash
chmod +x run_all.sh
./run_all.sh
```

## Monitoring Jobs

### Check Job Status

```bash
# View your jobs
squeue -u $USER

# View specific job
squeue -j <job_id>

# View detailed job info
scontrol show job <job_id>

# View completed job statistics
sacct -j <job_id> --format=JobID,JobName,Partition,State,Elapsed,MaxRSS
```

### View Job Output

```bash
# Tail live output
tail -f logs/gnb_vaex_<job_id>.out

# View errors
cat logs/gnb_vaex_<job_id>.err

# Monitor resource usage
sstat -j <job_id> --format=AveCPU,MaxRSS,MaxVMSize
```

### Cancel Jobs

```bash
# Cancel single job
scancel <job_id>

# Cancel all your jobs
scancel -u $USER

# Cancel jobs by name
scancel -n pypair_gnb_vaex
```

## Performance Optimization

### 1. Memory Optimization

**Vaex** - Use minimal memory allocation:
```bash
#SBATCH --mem=8GB  # Sufficient for out-of-core processing
```

**Polars** - Allocate based on dataset size:
```bash
#SBATCH --mem=32GB  # For S1-S4 datasets
#SBATCH --mem=64GB  # For larger combined datasets
```

### 2. CPU Optimization

Utilize multiple cores with parallel processing:

```python
# In Python script
import os
os.environ['OMP_NUM_THREADS'] = '8'  # Match --cpus-per-task

# Polars automatically uses all available cores
import polars as pl
pl.Config.set_global_string_cache(True)

# Vaex also parallelizes automatically
import vaex
vaex.set_parallel(True)
```

### 3. I/O Optimization

Store data on scratch for fast access:

```bash
# Good: Fast parallel file system
/scratch/k1234567/PyPair/data/

# Bad: Slow NFS-mounted home directory
/users/k1234567/PyPair/data/
```

Use HDF5 for efficient random access:
```python
# Convert CSV to HDF5 once
vaex.from_csv("large_file.csv", convert=True)
# Subsequent reads are 10-100x faster
df = vaex.open("large_file.hdf5")
```

## Troubleshooting

### Common Issues

#### 1. Out of Memory (OOM)

**Symptoms**:
```
slurmstepd: error: Detected 1 oom-kill event(s) in StepId=12345.batch
```

**Solution**:
- Increase `--mem` allocation
- Use Vaex instead of Polars/pandas
- Enable chunked processing

#### 2. Timeout

**Symptoms**:
```
CANCELLED AT 2025-01-24T12:00:00 DUE TO TIME LIMIT
```

**Solution**:
- Increase `--time` limit
- Optimize training code
- Use checkpointing for long jobs

#### 3. Module Load Failures

**Symptoms**:
```
ModuleCmd_Load.c(213):ERROR:105: Unable to locate a modulefile for 'anaconda3'
```

**Solution**:
```bash
# List available modules
module spider anaconda

# Load exact version
module load apps/anaconda3/2024.10/bin
```

#### 4. Conda Environment Not Found

**Symptoms**:
```
CondaEnvironmentError: environment does not exist: pypair
```

**Solution**:
```bash
# Recreate environment
conda create -n pypair python=3.10
conda activate pypair
pip install -e .
```

## Best Practices

### 1. Use Scratch for Large Files

```bash
# Store training data on scratch
/scratch/k1234567/PyPair/data/

# Keep code in home directory
/users/k1234567/PyPair/src/

# Save models to scratch (backed up separately)
/scratch/k1234567/PyPair/models/
```

### 2. Request Appropriate Resources

Don't over-allocate (wastes cluster resources):

```bash
# Good: Right-sized for Vaex
#SBATCH --mem=8GB
#SBATCH --time=02:00:00

# Bad: Wasteful over-allocation
#SBATCH --mem=256GB
#SBATCH --time=72:00:00
```

### 3. Use Job Arrays for Multiple Datasets

Train on S1-S4 simultaneously:

```bash
#!/bin/bash
#SBATCH --array=1-4
#SBATCH --job-name=gnb_array

DATASET="S${SLURM_ARRAY_TASK_ID}"
python train_model.py --dataset $DATASET
```

### 4. Enable Email Notifications

```bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=t.robinson@ulster.ac.uk
```

### 5. Log Everything

```bash
#SBATCH --output=logs/job_%j_%x.out
#SBATCH --error=logs/job_%j_%x.err

# In Python
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'logs/training_{os.getenv("SLURM_JOB_ID")}.log'),
        logging.StreamHandler()
    ]
)
```

## Production Deployment

### Automated Pipeline

Create a cron job for regular model updates:

```bash
# crontab entry (run weekly)
0 2 * * 0 /scratch/k1234567/PyPair/scripts/weekly_retrain.sh
```

**weekly_retrain.sh**:
```bash
#!/bin/bash
cd /scratch/k1234567/PyPair
git pull origin main
./run_all.sh
```

### Model Versioning

```bash
# Save models with timestamps
model_version=$(date +%Y%m%d_%H%M%S)
joblib.dump(model, f"models/gnb_vaex_{model_version}.joblib")

# Keep last 10 versions
ls -t models/gnb_vaex_*.joblib | tail -n +11 | xargs rm
```

## Summary Checklist

- [ ] SSH access to NIHPC configured
- [ ] Conda environment created with Python 3.10
- [ ] PyPair installed with UV/pip
- [ ] Training data transferred to `/scratch/`
- [ ] Data converted to HDF5 format
- [ ] Job scripts reviewed and customized
- [ ] Email notifications configured
- [ ] Test job submitted and completed successfully
- [ ] Full training pipeline executed via `run_all.sh`
- [ ] Benchmark results reviewed
- [ ] Production model deployed

---

**Document Version**: 1.0
**Last Updated**: 2025-01-24
**Cluster**: NIHPC (Northern Ireland HPC)
**Contact**: t.robinson@ulster.ac.uk
