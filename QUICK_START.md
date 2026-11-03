# PyPair Quick Start Guide

This is a condensed guide to get you up and running with PyPair quickly.

## 1. Environment Setup (5 minutes)

### Local Development

```bash
# Navigate to project
cd "c:\Users\Tony\OneDrive - Ulster University\Projects\PyPair"

# Install UV package manager
pip install uv

# Create virtual environment with Python 3.10
uv venv --python 3.10

# Activate environment
# Windows:
.venv\Scripts\activate
# Linux/Mac:
source .venv/bin/activate

# Install PyPair
uv pip install -e .

# Verify installation
python -c "import vaex; import polars; print('Success!')"
```

### HPC (NIHPC)

```bash
# SSH to cluster
ssh k1234567@login.hpc.ac.uk

# Load modules
module load apps/anaconda3/2024.10/bin

# Create environment
conda create -n pypair python=3.10 -y
conda activate pypair

# Install PyPair
cd /scratch/k1234567/PyPair
pip install uv
uv pip install -e ".[hpc]"
```

## 2. Test Installation (2 minutes)

```bash
# Test imports
python << 'EOF'
import vaex
import polars as pl
import numpy as np
from accelerated_GNB.GNB_polars.GNB_polars import GaussianNBPolars
from accelerated_GNB.GNB_vaex.GNB_HDF5_optimised import VaexGaussianNB

print(f"✓ Vaex {vaex.__version__}")
print(f"✓ Polars {pl.__version__}")
print(f"✓ NumPy {np.__version__}")
print("✓ All imports successful!")
EOF
```

## 3. Convert Data to HDF5 (One-Time, 10 minutes)

```python
# convert_data.py
import vaex
import glob

for csv_file in glob.glob("data/reference_samples/*.csv"):
    print(f"Converting {csv_file}...")
    df = vaex.from_csv(csv_file, convert=True, chunk_size=1_000_000)
    hdf5_file = csv_file.replace(".csv", ".hdf5")
    df.export_hdf5(hdf5_file)
    print(f"✓ Saved to {hdf5_file}")
```

Run it:
```bash
python convert_data.py
```

## 4. Train Your First Model (2 minutes)

### Polars (Fast, In-Memory)

```python
import polars as pl
from accelerated_GNB.GNB_polars.GNB_polars import GaussianNBPolars
import joblib

# Load data
X_train = pl.read_csv("data/reference_samples/S1_extended_encoded_x_train.csv")
y_train = pl.read_csv("data/reference_samples/S1_extended_encoded_y_train.csv")["label"]

# Train model
model = GaussianNBPolars()
model.fit(X_train, y_train)

# Save model
joblib.dump(model, "models/gnb_polars_test.joblib")
print("✓ Model trained and saved!")
```

### Vaex (Memory-Efficient)

```python
import vaex
from accelerated_GNB.GNB_vaex.GNB_HDF5_optimised import VaexGaussianNB
import joblib

# Load data (memory-mapped, no RAM usage)
X_train = vaex.open("data/reference_samples/S1_extended_encoded_x_train.hdf5")
y_train = vaex.open("data/reference_samples/S1_extended_encoded_y_train.hdf5")["label"]

# Train model
model = VaexGaussianNB()
model.fit(X_train, y_train)

# Save model
joblib.dump(model, "models/gnb_vaex_test.joblib")
print("✓ Model trained and saved!")
```

## 5. Run Predictions (1 minute)

```python
import vaex
import joblib
from sklearn.metrics import f1_score, accuracy_score

# Load model
model = joblib.load("models/gnb_vaex_test.joblib")

# Load test data
X_test = vaex.open("data/reference_samples/S1_extended_encoded_X_test.hdf5")
y_test = vaex.open("data/reference_samples/S1_extended_encoded_y_test.hdf5")["label"]

# Predict
y_pred = model.predict_chunked(X_test, chunk_size=10000)

# Evaluate
f1 = f1_score(y_test.values, y_pred, average='weighted')
acc = accuracy_score(y_test.values, y_pred)

print(f"F1 Score: {f1:.4f}")
print(f"Accuracy: {acc:.4f}")
```

## 6. Run Benchmark (5 minutes)

```bash
python models/benchmark.py \
    --implementation vaex \
    --dataset S1 \
    --runs 1 \
    --output results/quick_benchmark.csv

# View results
cat results/quick_benchmark.csv
```

## 7. HPC Job Submission (2 minutes)

```bash
# On NIHPC cluster
cd /scratch/k1234567/PyPair

# Submit data preparation
sbatch data_prep.sh

# Check status
squeue -u $USER

# View output (after job completes)
cat logs/data_prep_*.out
```

## 8. Full Genome Training (HPC, 4 hours)

```bash
# Submit full genome training job
sbatch train_full_genome.sh

# Monitor progress
tail -f logs/gnb_full_genome_*.out

# Expected output:
# Training on full human reference genome (GRCh38)...
# Training completed in 2847.3 seconds (47.5 minutes)
# Model saved!
```

## Common Commands

### Development

```bash
# Activate environment
# Windows:
.venv\Scripts\activate
# Linux/Mac:
source .venv/bin/activate

# Run tests (when implemented)
pytest tests/

# Format code
black src/ tests/
ruff check src/ tests/

# Type checking
mypy src/
```

### HPC

```bash
# Check job queue
squeue -u $USER

# View job details
scontrol show job <job_id>

# Cancel job
scancel <job_id>

# View logs
tail -f logs/*.out

# Check resource usage
sacct -j <job_id> --format=JobID,Elapsed,MaxRSS,State
```

## Troubleshooting

### Issue: Import errors

```bash
# Verify Python version
python --version  # Should be 3.10.x

# Verify packages
pip list | grep -E "vaex|polars|numpy"

# Reinstall
uv pip install -e . --force-reinstall
```

### Issue: Out of memory

```python
# Use Vaex instead of Polars
import vaex  # Handles unlimited data size
X = vaex.open("large_file.hdf5")  # Memory-mapped, no RAM
```

### Issue: HDF5 file not found

```bash
# Convert CSV to HDF5 first
python convert_data.py

# Verify files exist
ls -lh data/reference_samples/*.hdf5
```

### Issue: SLURM job fails

```bash
# Check error log
cat logs/job_name_<job_id>.err

# Common fixes:
# 1. Increase memory: #SBATCH --mem=16GB
# 2. Increase time: #SBATCH --time=04:00:00
# 3. Check conda env: conda activate pypair
```

## Key Files

| File | Purpose |
|------|---------|
| [README.md](README.md) | Full project documentation |
| [pyproject.toml](pyproject.toml) | Package configuration |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) | Complete project review |
| [docs/overfitting_analysis.md](docs/overfitting_analysis.md) | Why GNB is optimal |
| [docs/gnb_implementations.md](docs/gnb_implementations.md) | Implementation details |
| [docs/hpc_deployment.md](docs/hpc_deployment.md) | HPC deployment guide |
| [docs/benchmarking.md](docs/benchmarking.md) | Benchmarking guide |

## Next Steps

1. ✅ Environment setup
2. ✅ Test installation
3. ✅ Convert data to HDF5
4. ✅ Train first model
5. ⏳ Run full benchmarks
6. ⏳ Train on full genome (HPC)
7. ⏳ Publish results

## Getting Help

- **Documentation**: See [README.md](README.md) and [docs/](docs/)
- **Issues**: Check [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) section 9 (Risks & Mitigations)
- **Contact**: t.robinson@ulster.ac.uk

---

**Last Updated**: 2025-01-24
**Time to Complete**: ~30 minutes (excluding HPC jobs)
**Status**: Ready to use
