# PyPair: HPC Deployment Guide

## Quick Start

**Single command to run full pipeline**:
```bash
sbatch scripts/hpc/pypair_full_pipeline.sh
```

This executes:
1. Data preparation
2. Model training (S1-S4)
3. Performance profiling
4. Scaling analysis (KNN vs GNB)
5. Visualization generation

**Expected runtime**: 24-48 hours
**Memory required**: 64GB
**CPUs required**: 16

---

## Pipeline Stages

### Stage 1: Data Preparation
```bash
# If needed (usually already done)
python scripts/training_data_preparation.py
```

**Output**:
- `data/reference_samples/cleaned/S1_extended_encoded_*.csv`
- `data/reference_samples/cleaned/S2_extended_encoded_*.csv`
- `data/reference_samples/cleaned/S3_extended_encoded_*.csv`
- `data/reference_samples/cleaned/S4_extended_encoded_*.csv`

### Stage 2: Model Training & Validation
```bash
python scripts/comprehensive_model_analysis.py
```

**Models trained**:
- Gaussian NB ✅
- KNN ✅
- Nearest Centroid ✅
- Random Forest (memory limited)
- Decision Tree (fails to learn)
- XGBoost (fails to learn)
- LGBM (fails to learn)

**Output**:
- `results/comprehensive_model_analysis.csv`
- `models/model_weights/*.joblib`

### Stage 3: Performance Profiling
```bash
python scripts/profile_seed_mapping_fast.py
```

**Metrics captured**:
- Time per seed (ms)
- Throughput (seeds/sec)
- BWA speedup comparison

**Output**:
- `results/seed_mapping_profile_fast.csv`
- `results/seed_mapping_summary_fast.txt`

### Stage 4: Scaling Analysis
```bash
python scripts/scaling_analysis.py
```

**Hypothesis tested**:
> GNB will outperform KNN at genome scale due to lower memory requirements

**Data scales**:
1. S1 (44.8k samples)
2. S1+S2 (~90k samples)
3. S1+S2+S3 (~135k samples)
4. S1+S2+S3+S4 (~180k samples)

**Output**:
- `results/scaling_analysis.csv`

### Stage 5: Visualizations
```bash
python scripts/create_viva_plots.py
python scripts/create_scaling_plots.py
```

**Output**:
- `docs/figures/00_summary_infographic.png`
- `docs/figures/01_model_performance.png`
- `docs/figures/02_overfitting_gap.png`
- `docs/figures/03_bwa_speedup.png`
- `docs/figures/04_batch_performance.png`
- `docs/figures/05_cross_dataset_consistency.png`
- `docs/figures/scaling_01_memory.png`
- `docs/figures/scaling_02_memory_per_sample.png`
- `docs/figures/scaling_03_training_time.png`
- `docs/figures/scaling_04_comparison.png`
- `docs/figures/scaling_05_genome_projection.png`

---

## Individual HPC Job Submissions

### Run S1-S4 Validation Only
```bash
sbatch scripts/hpc/pypair_validation_only.sh
```

### Run Scaling Analysis Only
```bash
sbatch scripts/hpc/pypair_scaling_only.sh
```

### Run Chr21 Training
```bash
# Step 1: Add extended features (if needed)
sbatch scripts/hpc/add_features_chr21_hpc.sh

# Step 2: Train models
sbatch scripts/hpc/train_chr21_hpc.sh
```

---

## Memory Requirements by Scale

| Data Scale | Samples | KNN Memory | GNB Memory | Winner |
|------------|---------|------------|------------|--------|
| S1 | 44.8k | ~4 MB | ~0.88 MB | GNB (4.5x smaller) |
| S1+S2 | ~90k | ~8 MB | ~0.88 MB | GNB (9x smaller) |
| S1-S4 | ~180k | ~16 MB | ~0.88 MB | GNB (18x smaller) |
| chr21 | ~80M | ~7 GB | ~0.88 MB | GNB (8000x smaller) |
| Genome | ~3B | ~260 GB | ~0.88 MB | **GNB only viable** |

**Conclusion**: KNN memory grows linearly O(n), GNB stays constant O(c×d)

---

## HPC Configuration (NI-HPC Ulster)

### Partition Information
```
Partition: k2-hipri
Max nodes: 20
Cores per node: 32
Memory per node: 256GB
Max walltime: 168 hours
```

### Module Loading
```bash
module purge
module load python/3.10
module load gcc/11.2.0
```

### Python Environment Setup
```bash
# Create virtual environment (first time only)
python3 -m venv ~/venvs/pypair
source ~/venvs/pypair/bin/activate

# Install dependencies
pip install --upgrade pip
pip install pandas numpy scikit-learn joblib matplotlib seaborn psutil
pip install biopython pysam  # For FASTA/BAM handling
```

---

## Job Monitoring

### Check job status
```bash
squeue -u $USER
```

### View job output (live)
```bash
tail -f logs/pypair_pipeline_<JOB_ID>.out
```

### Cancel job
```bash
scancel <JOB_ID>
```

### Check resource usage
```bash
sacct -j <JOB_ID> --format=JobID,JobName,Partition,State,Elapsed,MaxRSS,MaxVMSize
```

---

## Troubleshooting

### Problem: Out of Memory
**Solution**: Increase `#SBATCH --mem=128G` or reduce dataset size

### Problem: Job Times Out
**Solution**: Increase `#SBATCH --time=72:00:00`

### Problem: Module Not Found
**Solution**: Check module availability with `module avail python`

### Problem: Permission Denied
**Solution**: Make script executable: `chmod +x scripts/hpc/*.sh`

---

## Expected Results

### Performance Metrics (from S1-S4 validation)

**KNN**:
- Accuracy: 99.88% F1
- Speed: 0.084 ms/seed (35.7x faster than BWA)
- Memory: ~16 MB (S1-S4 combined)

**Gaussian NB**:
- Accuracy: 95.29% F1
- Speed: 0.84 ms/seed (3.6x faster than BWA)
- Memory: ~0.88 MB (constant across scales)

**Nearest Centroid**:
- Accuracy: 99.88% F1
- Memory: ~1-2 MB per 10k samples

### Scaling Hypothesis Test

**Expected outcome**:
```
S1 (44.8k):     KNN competitive with GNB
S1-S4 (180k):   KNN still manageable, GNB more efficient
chr21 (80M):    KNN struggles (7GB), GNB maintains (0.88MB)
Genome (3B):    KNN infeasible (260GB), GNB only option
```

**Hypothesis**: **GNB will be the only viable model at genome scale**

---

## File Structure After Pipeline

```
PyPair/
├── results/
│   ├── comprehensive_model_analysis.csv     # S1-S4 validation
│   ├── seed_mapping_profile_fast.csv        # Performance metrics
│   ├── scaling_analysis.csv                 # Scaling experiment
│   └── pipeline_complete.txt                # Summary report
├── models/model_weights/
│   ├── knn_s1.joblib
│   ├── knn_s2.joblib
│   ├── knn_s3.joblib
│   ├── knn_s4.joblib
│   ├── gaussian_nb_s1.joblib
│   ├── gaussian_nb_s2.joblib
│   ├── gaussian_nb_s3.joblib
│   └── gaussian_nb_s4.joblib
├── docs/figures/
│   ├── 00_summary_infographic.png
│   ├── 01_model_performance.png
│   ├── 02_overfitting_gap.png
│   ├── 03_bwa_speedup.png
│   ├── 04_batch_performance.png
│   ├── 05_cross_dataset_consistency.png
│   ├── scaling_01_memory.png
│   ├── scaling_02_memory_per_sample.png
│   ├── scaling_03_training_time.png
│   ├── scaling_04_comparison.png
│   └── scaling_05_genome_projection.png
└── logs/
    └── pypair_pipeline_<JOB_ID>.out
```

---

## Next Steps

### For Viva (Next Week)
1. ✅ Run full pipeline: `sbatch scripts/hpc/pypair_full_pipeline.sh`
2. ✅ Review results in `results/`
3. ✅ View all plots in `docs/figures/`
4. ✅ Prepare talking points from [REORGANIZATION_COMPLETE.md](REORGANIZATION_COMPLETE.md)

### For Future Work
1. Scale to chr21 (80M k-mers)
2. Test genome-scale (3B k-mers) - **GNB only**
3. Optimize batching for production deployment
4. Integrate with clinical WGS pipelines

---

**Status**: Ready for HPC deployment
**Contact**: Tony Robinson (Ulster University)
**Last Updated**: 2025-10-24
