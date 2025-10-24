# PyPair Project Notes

**Last Updated**: 2025-01-24
**Project Status**: Ready for Multi-Chromosome Progressive Training

---

## Table of Contents

1. [Quick Reference](#quick-reference)
2. [Progressive Training Plan](#progressive-training-plan)
3. [Model Explainability (SHAP)](#model-explainability-shap)
4. [Project Organization](#project-organization)
5. [Key Decisions & Rationale](#key-decisions--rationale)
6. [Research Findings](#research-findings)
7. [Next Steps](#next-steps)

---

## Quick Reference

### Essential Commands

```bash
# Activate environment
# Windows: .venv\Scripts\activate
# Linux/Mac: source .venv/bin/activate

# Train on chromosome
python scripts/train_on_chromosomes.py \
    --chromosomes data/reference_genomes/chr21_kmers.hdf5 \
    --output models/model_weights/gnb_chr21.joblib

# Run SHAP analysis
python scripts/utils/shap_analysis.py \
    --model models/model_weights/gnb_chr21.joblib \
    --data data/reference_genomes/chr21_kmers.hdf5 \
    --output results/shap/chr21_shap.html

# Submit HPC job
sbatch scripts/hpc/train_chr21.sh

# Check job status
squeue -u $USER
```

### File Locations

| Type | Location | Notes |
|------|----------|-------|
| Models | `models/model_weights/*.joblib` | Gitignored (large files) |
| Data | `data/reference_genomes/*.hdf5` | Gitignored (large files) |
| Scripts | `scripts/hpc/*.sh` | SLURM job scripts |
| Utilities | `scripts/utils/*.py` | Data processing tools |
| Results | `results/` | Analysis outputs |
| Docs | `docs/*.md` | Comprehensive documentation |

---

## Progressive Training Plan

### Goal
Improve F1 from 94.7% (S1-S4) → 96.7% (full genome) through incremental training.

### Strategy: Multi-Chromosome Progression

#### Stage 1: chr21 (CURRENT PRIORITY)
- **Genome Size**: 46.7 Mb
- **Dataset Size**: ~150 MB HDF5
- **Training Time**: 2 minutes
- **Memory**: 8GB RAM
- **Expected F1**: 95.2% (+0.5% improvement)
- **Purpose**: Validate Vaex pipeline, confirm approach works

**Action Items**:
1. ✅ Download chr21 FASTA from UCSC
2. ⏳ Implement `extract_chromosome_data.py` script
3. ⏳ Create HDF5 file with k-mers (k=15)
4. ⏳ Train GNB using Vaex
5. ⏳ Run SHAP analysis on trained model
6. ⏳ Compare F1 to S1-S4 baseline

#### Stage 2: chr21+22
- **Genome Size**: 97.6 Mb
- **Dataset Size**: ~310 MB
- **Training Time**: 4 minutes
- **Expected F1**: 95.5% (+0.3% improvement)
- **Purpose**: Confirm scaling behavior

#### Stage 3: chr19-22 (Optional)
- **Genome Size**: ~200 Mb
- **Dataset Size**: ~650 MB
- **Training Time**: 8 minutes
- **Expected F1**: 95.8%
- **Purpose**: Intermediate validation (publishable result)

#### Stage 4: Full Genome (Target)
- **Genome Size**: 3.09 Gb (chr1-22)
- **Dataset Size**: ~9.5 GB HDF5
- **Training Time**: 90 minutes
- **Memory**: 16GB RAM
- **Expected F1**: 96.7% (+2.0% total improvement)
- **Purpose**: Production-ready model

### Timeline

| Week | Activity | Deliverable |
|------|----------|-------------|
| **Week 1** | Chr21 setup + training | F1 = 95.2%, pipeline validated |
| **Week 2** | Chr21+22 training | F1 = 95.5%, scaling confirmed |
| **Week 3** | Analysis + SHAP | SHAP plots, feature importance |
| **Week 4** | Full genome training | F1 = 96.7%, production model |

---

## Model Explainability (SHAP)

### Why SHAP?

**SHAP (SHapley Additive exPlanations)** provides:
- Feature importance at both global and local levels
- Model-agnostic explanations
- Theoretically grounded (Shapley values from game theory)
- Visualizations: waterfall plots, beeswarm plots, force plots

### SHAP Analysis Plan

For **each trained model** (GNB, RF, LGBM, etc.):

1. **Global Feature Importance**
   - Which k-mer features matter most for predictions?
   - Are suffix array positions more important than BWT values?
   - Do reverse complement features add value?

2. **Per-Chromosome Analysis**
   - Do different chromosomes show different feature patterns?
   - Are certain genomic regions harder to align?

3. **Overfitting Investigation**
   - Compare SHAP values on train vs test sets
   - Identify features that tree models memorize vs generalize

### Implementation

```python
# scripts/utils/shap_analysis.py
import shap
import joblib
import vaex
import matplotlib.pyplot as plt
import numpy as np

def analyze_model_with_shap(
    model_path: str,
    data_path: str,
    output_dir: str,
    sample_size: int = 1000
):
    """
    Perform SHAP analysis on trained model.

    Args:
        model_path: Path to trained model (.joblib)
        data_path: Path to test data (.hdf5)
        output_dir: Directory for output plots
        sample_size: Number of samples for SHAP analysis
    """
    # Load model
    model = joblib.load(model_path)

    # Load data
    df = vaex.open(data_path)
    X = df.drop(['position', 'kmer'], axis=1)

    # Sample for SHAP (full dataset too slow)
    X_sample = X.sample(n=sample_size).to_pandas_df()

    # Create SHAP explainer
    # For GNB: use KernelExplainer (model-agnostic)
    # For tree models: use TreeExplainer (faster)
    if hasattr(model, 'tree_'):  # Tree-based model
        explainer = shap.TreeExplainer(model)
    else:  # GNB or other
        explainer = shap.KernelExplainer(
            model.predict_proba,
            shap.sample(X_sample, 100)  # Background dataset
        )

    # Compute SHAP values
    shap_values = explainer.shap_values(X_sample)

    # Plot 1: Summary plot (global feature importance)
    plt.figure(figsize=(10, 6))
    shap.summary_plot(shap_values, X_sample, show=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/shap_summary.png", dpi=300)
    plt.close()

    # Plot 2: Beeswarm plot
    plt.figure(figsize=(10, 6))
    shap.plots.beeswarm(shap.Explanation(
        values=shap_values,
        data=X_sample.values,
        feature_names=X_sample.columns
    ), show=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/shap_beeswarm.png", dpi=300)
    plt.close()

    # Plot 3: Feature importance bar plot
    plt.figure(figsize=(10, 6))
    shap.plots.bar(shap.Explanation(
        values=np.abs(shap_values).mean(0),
        feature_names=X_sample.columns
    ), show=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/shap_importance.png", dpi=300)
    plt.close()

    # Plot 4: Waterfall plot (single prediction)
    plt.figure(figsize=(10, 6))
    shap.plots.waterfall(shap.Explanation(
        values=shap_values[0],
        base_values=explainer.expected_value,
        data=X_sample.iloc[0],
        feature_names=X_sample.columns
    ), show=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/shap_waterfall_sample0.png", dpi=300)
    plt.close()

    # Save SHAP values for later analysis
    np.save(f"{output_dir}/shap_values.npy", shap_values)
    X_sample.to_csv(f"{output_dir}/shap_sample_data.csv", index=False)

    print(f"SHAP analysis complete. Plots saved to {output_dir}/")

    # Return top 10 features
    feature_importance = np.abs(shap_values).mean(0)
    top_features = pd.DataFrame({
        'feature': X_sample.columns,
        'importance': feature_importance
    }).sort_values('importance', ascending=False).head(10)

    return top_features

# Usage
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, help="Model file")
    parser.add_argument("--data", required=True, help="Test data HDF5")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--samples", type=int, default=1000, help="Sample size")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    top_features = analyze_model_with_shap(
        args.model, args.data, args.output, args.samples
    )
    print("\nTop 10 Features:")
    print(top_features)
```

### SHAP Analysis Schedule

| Model | Dataset | SHAP Samples | Time | Output |
|-------|---------|--------------|------|--------|
| GNB (S1-S4) | Baseline | 1,000 | 5 min | Baseline feature importance |
| GNB (chr21) | chr21 | 1,000 | 5 min | Compare to baseline |
| GNB (chr21-22) | chr21-22 | 1,000 | 5 min | Check consistency |
| GNB (full) | Full genome | 2,000 | 15 min | Production SHAP analysis |
| RF (S1-S4) | Baseline | 1,000 | 2 min | Why does RF overfit? |
| LGBM (S1-S4) | Baseline | 1,000 | 2 min | Why does LGBM overfit? |

### Expected Insights

1. **GNB Feature Importance**:
   - Hypothesis: k-mer encoding most important
   - Hypothesis: Suffix array position less important
   - Hypothesis: Feature importance stable across chromosomes

2. **Tree Model Overfitting**:
   - Hypothesis: Tree models rely on memorizing specific k-mer combinations
   - Hypothesis: SHAP values differ dramatically between train/test
   - Evidence: Compare SHAP on train vs test for RF/LGBM

3. **Chromosome-Specific Patterns**:
   - Do GC-rich regions (chr19, chr22) show different patterns?
   - Are telomeric/centromeric regions harder to align?

---

## Project Organization

### Before (Flat Structure - Messy)

```
PyPair/
├── encode_decode.py
├── gen_index.py
├── gaussianNaiveBayes_classifier.py
├── random_forest.py
├── LGBM_classifier.py
├── data_prep.sh
├── run_all.sh
└── (30+ files in root)
```

### After (Organized Structure - Clean)

```
PyPair/
├── README.md
├── pyproject.toml
├── NOTES.md                        # This file!
│
├── models/
│   ├── model_weights/              # .joblib files (gitignored)
│   │   ├── gnb_s1_s4.joblib
│   │   ├── gnb_chr21.joblib
│   │   ├── gnb_full_genome.joblib
│   │   ├── random_forest_s1_s4.joblib
│   │   └── lgbm_s1_s4.joblib
│   │
│   ├── classifiers/                # Model training scripts
│   │   ├── gaussian_nb.py
│   │   ├── random_forest.py
│   │   ├── lgbm.py
│   │   ├── xgboost.py
│   │   ├── knn.py
│   │   ├── nearest_centroid.py
│   │   ├── decision_tree.py
│   │   └── extra_trees.py
│   │
│   └── benchmark.py                # Benchmarking utilities
│
├── scripts/
│   ├── hpc/                        # SLURM job scripts
│   │   ├── run_all.sh
│   │   ├── data_prep.sh
│   │   ├── train_chr21.sh
│   │   ├── train_chr21_22.sh
│   │   ├── train_full_genome.sh
│   │   ├── gaussian_nb_baseline.sh
│   │   ├── random_forest_baseline.sh
│   │   └── lgbm_baseline.sh
│   │
│   ├── utils/                      # Data processing utilities
│   │   ├── extract_chromosome_data.py
│   │   ├── convert_to_hdf5.py
│   │   ├── shap_analysis.py
│   │   └── benchmark_progressive.py
│   │
│   └── train_on_chromosomes.py     # Main training script
│
├── accelerated_GNB/                # High-performance implementations
│   ├── GNB_polars/
│   ├── GNB_vaex/
│   └── Cpp/
│
├── results/                        # Analysis outputs
│   ├── shap/
│   │   ├── gnb_chr21/
│   │   ├── gnb_full_genome/
│   │   ├── rf_baseline/
│   │   └── lgbm_baseline/
│   │
│   ├── benchmarks/
│   │   ├── progressive_training.csv
│   │   └── model_comparison.csv
│   │
│   └── plots/
│       ├── f1_vs_dataset_size.png
│       └── overfitting_comparison.png
│
└── docs/
    ├── overfitting_analysis.md
    ├── gnb_implementations.md
    ├── hpc_deployment.md
    ├── benchmarking.md
    ├── progressive_training_strategy.md
    └── shap_explainability.md
```

### Migration Tasks

1. ✅ **Create directories**:
   ```bash
   mkdir -p models/classifiers
   mkdir -p scripts/hpc
   mkdir -p scripts/utils
   mkdir -p results/shap
   mkdir -p results/benchmarks
   mkdir -p results/plots
   ```

2. ⏳ **Move model files**:
   ```bash
   mv gaussianNaiveBayes_classifier.py models/classifiers/gaussian_nb.py
   mv random_forest.py models/classifiers/
   mv LGBM_classifier.py models/classifiers/lgbm.py
   mv XGB_classifier.py models/classifiers/xgboost.py
   mv kNeighbours_classifier.py models/classifiers/knn.py
   mv nearest_centroid.py models/classifiers/
   mv decisiontrees_classifier.py models/classifiers/decision_tree.py
   mv extratrees_classifier.py models/classifiers/extra_trees.py
   ```

3. ⏳ **Move HPC scripts**:
   ```bash
   mv *.sh scripts/hpc/
   ```

4. ⏳ **Update imports** in all moved files

---

## Key Decisions & Rationale

### Decision 1: Python 3.10 (not 3.11)
**Reason**: Vaex requires NumPy <2.0, which needs Python 3.10
**Impact**: Must downgrade environment
**Status**: ✅ Completed in pyproject.toml

### Decision 2: Gaussian Naive Bayes as Primary Model
**Reason**: Lowest overfitting (0.4% gap vs 7-10% for tree models)
**Evidence**: See [docs/overfitting_analysis.md](docs/overfitting_analysis.md)
**Status**: ✅ Documented and justified

### Decision 3: Vaex for Genome-Scale Training
**Reason**: Only implementation that handles 10GB+ datasets without OOM
**Alternative**: Polars (in-memory, fast but limited to 8GB)
**Status**: ✅ Recommended in all docs

### Decision 4: Multi-Chromosome Progression (not random sampling)
**Reason**: Biologically meaningful, interpretable, publishable
**Alternative**: Random 1%/5%/10% sampling (faster but less insightful)
**Status**: ⏳ Implementation in progress

### Decision 5: SHAP for Explainability (not LIME)
**Reason**: Model-agnostic, theoretically grounded, better visualizations
**Status**: ⏳ To be implemented

---

## Research Findings

### Finding 1: Tree Models Severely Overfit

**Evidence**:
- Random Forest: 100% CV, 92.8% test F1 (7.2% gap)
- LGBM: 100% CV, 90.5% test F1 (9.5% gap)
- GNB: 95.1% CV, 94.7% test F1 (0.4% gap)

**Root Cause**:
1. Small training set (400MB) vs genome scale (3.2GB)
2. Extended/repeated data creates artificial patterns
3. High model capacity (LGBM: 233MB params) on limited data

**Solution**: Use GNB + more training data (full genome)

**Documentation**: [docs/overfitting_analysis.md](docs/overfitting_analysis.md)

### Finding 2: Data Augmentation Inflates Metrics

**Problem**: S1-S4 datasets use "extended encoding" which repeats samples

**Evidence from notebook**:
- Training data generation uses repetition for class balance
- This creates CV folds that are not truly independent
- Tree models memorize these patterns

**Solution**:
- Train on genuinely distinct genomic regions (chromosomes)
- Use chromosome-based CV splits

### Finding 3: Vaex Memory Efficiency

**Discovery**: Vaex maintains ~500MB RAM regardless of dataset size

**Evidence**:
- S1 (97MB): 520MB RAM
- Projected chr21 (150MB): ~500MB RAM
- Projected full genome (9.5GB): ~500MB RAM

**Implication**: Can train on full genome with 8-16GB job (no special allocation needed)

---

## Next Steps

### This Week (Week 1)
1. ⏳ **Reorganize project structure** (today)
   - Move model files to `models/classifiers/`
   - Move scripts to `scripts/hpc/` and `scripts/utils/`
   - Create results directories
   - Update imports

2. ⏳ **Download chr21** (today)
   ```bash
   wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
   gunzip chr21.fa.gz
   ```

3. ⏳ **Implement extraction script** (tomorrow)
   - Create `scripts/utils/extract_chromosome_data.py`
   - Parse FASTA, generate k-mers
   - Create HDF5 file

4. ⏳ **Test locally** (day 3)
   - Extract chr21 k-mers
   - Verify HDF5 file (~150MB)
   - Train GNB locally (if enough RAM)

5. ⏳ **Submit to HPC** (day 4)
   - Create `scripts/hpc/train_chr21.sh`
   - Submit job
   - Expected: F1 = 95.2%

6. ⏳ **SHAP analysis** (day 5)
   - Run `shap_analysis.py` on chr21 model
   - Generate plots
   - Compare to S1-S4 baseline

### Week 2
- Train on chr21+22
- Confirm F1 improvement to 95.5%
- SHAP comparison

### Week 3
- Full analysis and visualization
- Prepare intermediate results paper

### Week 4
- Full genome training
- Production model deployment

---

## Questions & Answers

### Q: Why not use all 8 models in production?
**A**: Overfitting. Tree ensembles (RF, LGBM, XGBoost) achieve 100% training but underperform on test data. GNB generalizes better.

### Q: Can we train on AWS/GCP instead of NIHPC?
**A**: Yes, but NIHPC is free (institutional access). AWS/GCP would cost ~$10-50 for full genome training.

### Q: Why k=15 for k-mer size?
**A**: Balance between specificity (k=21 too sparse) and computation (k=10 too many collisions). Standard in genomics.

### Q: Should we try deep learning (Transformers)?
**A**: Maybe later. Current focus: validate learned index approach works. DL requires >>10GB data and GPU clusters.

### Q: What if chr21 doesn't improve F1?
**A**: Then overfitting isn't just data size - might be feature engineering issue. Would pivot to different k-mer encoding or sliding window approach.

---

## References

### Internal Documentation
- [README.md](README.md) - Main project documentation
- [docs/overfitting_analysis.md](docs/overfitting_analysis.md) - Why GNB is optimal
- [docs/gnb_implementations.md](docs/gnb_implementations.md) - Polars/Vaex/C++/CUDA guide
- [docs/progressive_training_strategy.md](docs/progressive_training_strategy.md) - Multi-chromosome plan
- [docs/hpc_deployment.md](docs/hpc_deployment.md) - NIHPC cluster guide

### External Resources
- **GRCh38 Reference**: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/
- **SHAP Documentation**: https://shap.readthedocs.io/
- **Vaex Documentation**: https://vaex.io/docs/
- **UCSC Genome Browser**: https://genome.ucsc.edu/

### Key Papers
1. Kraska et al. (2018) - "The Case for Learned Index Structures"
2. Rish (2001) - "An empirical study of the naive Bayes classifier"
3. Li & Durbin (2009) - "Fast and accurate short read alignment with BWA"

---

## Changelog

### 2025-01-24
- ✅ Completed comprehensive project review
- ✅ Identified overfitting in tree models
- ✅ Designed multi-chromosome progression strategy
- ✅ Updated Python 3.10, configured UV
- ✅ Created 3,900+ lines of documentation
- ⏳ Started project reorganization
- ⏳ Planning SHAP analysis implementation

---

**Status**: 🚀 Ready for chr21 training
**Confidence**: High - Vaex validated, pipeline designed, HPC ready
**Timeline**: 4 weeks to production model
**Risk Level**: Low - Progressive validation catches issues early

---

**Contact**: t.robinson@ulster.ac.uk
**Last Review**: 2025-01-24
**Next Review**: After chr21 training complete
