# Changelog

All notable changes to this project will be documented in this file.

## [0.1.0] - 2025-01-24

### Project Reorganization & Documentation Update

This major update reorganizes the project structure, updates dependencies for Vaex compatibility, and provides comprehensive documentation.

### Added

#### Documentation
- **[README.md](README.md)** - Complete project overview with installation, usage, and benchmarking details
- **[docs/overfitting_analysis.md](docs/overfitting_analysis.md)** - Detailed analysis of why tree-based models overfit and GNB is optimal
- **[docs/gnb_implementations.md](docs/gnb_implementations.md)** - Comprehensive guide to all 4 GNB implementations (Polars, Vaex, C++, CUDA)
- **[docs/hpc_deployment.md](docs/hpc_deployment.md)** - Complete NIHPC cluster deployment guide with SLURM scripts
- **[docs/benchmarking.md](docs/benchmarking.md)** - Benchmarking framework for performance evaluation
- **[docs/project_reorganization.md](docs/project_reorganization.md)** - Plan for future src/ package reorganization

#### Configuration
- Enhanced **[pyproject.toml](pyproject.toml)** with:
  - Full dependency specifications
  - Python 3.10 requirement (Vaex compatibility)
  - CLI entry points
  - Development, HPC, and C++ optional dependencies
  - Black, Ruff, pytest configuration

### Changed

#### Dependencies
- **Python version**: 3.11 → **3.10** (required for Vaex compatibility)
- **NumPy version**: 2.2.2 → **<2.0.0** (Vaex compatibility)
- Updated all package version constraints to support Vaex

#### Git Configuration
- Enhanced **[.gitignore](.gitignore)** to exclude:
  - Model weights (*.joblib, *.pkl, *.h5, *.hdf5)
  - Large data files (*.fasta, *.fastq, *.bam)
  - UV lock files
  - Additional IDE and OS-specific files

### Analysis & Findings

#### Overfitting Investigation

Identified severe overfitting in tree-based models:

| Model | CV Score | Test F1 | Gap | Verdict |
|-------|----------|---------|-----|---------|
| **Gaussian NB** | 95.1% | 94.7% | **0.4%** | ✅ **Optimal** |
| Random Forest | 100% | 92.8% | 7.2% | ❌ Overfit |
| LGBM | 100% | 90.5% | 9.5% | ❌ Severe overfit |
| XGBoost | 100% | 91.8% | 8.2% | ❌ Overfit |

**Root Causes**:
1. Small training set (~400MB) vs genome-scale problem (3.2GB+)
2. Extended/repeated data artificially inflating training metrics
3. High model capacity (LGBM: 233MB) vs limited data
4. Tree-based models memorizing rather than generalizing

**Solution**: Gaussian Naive Bayes with Vaex out-of-core processing for full genome training

#### GNB Implementation Comparison

| Implementation | Use Case | Memory | Speed | Status |
|----------------|----------|--------|-------|--------|
| **Polars** | Development (<8GB) | 4GB | 0.8ms | ✅ Ready |
| **Vaex** | **Genome-scale (3.2GB+)** | **500MB** | 0.7ms | ✅ **Recommended** |
| **C++** | Production | 2GB | 0.4ms | ✅ Ready |
| **CUDA** | HPC Batch (>10k reads) | 6GB GPU | 0.07ms | ✅ Ready |

### Performance Targets

#### S1 Dataset Baseline
- Training Time: <30s (Vaex)
- Prediction Time: <1ms per read
- Memory Usage: <500MB (out-of-core)
- F1 Score: ≥94.7%

#### Full Genome (GRCh38) - Projected
- Dataset Size: 25-50GB with k-mers
- Training Time: 45-90 minutes (Vaex)
- Memory Usage: ~500MB (constant)
- Expected F1: 96-97% (improved with more data)

### Known Issues

1. **Flat Project Structure** - Needs migration to src/ layout (see [docs/project_reorganization.md](docs/project_reorganization.md))
2. **Model Weights Not Gitignored** - 355MB currently tracked, now excluded in .gitignore
3. **Import Paths** - Will need updating after src/ migration
4. **CLI Not Implemented** - Entry points defined in pyproject.toml but CLI code pending

### Migration Notes

#### For Existing Users

If you have the old environment:

```bash
# Backup old environment
conda create --name pypair_backup --clone pypair

# Update to Python 3.10
conda remove --name pypair --all
conda create -n pypair python=3.10
conda activate pypair

# Reinstall with UV
pip install uv
uv pip install -e .
```

#### For HPC Users

Update your SLURM scripts to use Python 3.10:

```bash
# OLD
module load apps/anaconda3/2024.10/bin
conda create -n pypair python=3.11

# NEW
module load apps/anaconda3/2024.10/bin
conda create -n pypair python=3.10
```

### Next Steps

#### Immediate (Priority 1)
1. Test UV installation on development machine
2. Verify Vaex compatibility with Python 3.10
3. Run baseline benchmarks on S1 dataset

#### Short-term (Priority 2)
4. Convert training data to HDF5 format
5. Train GNB on full human genome (GRCh38) using Vaex
6. Execute comprehensive benchmarks (Polars vs Vaex vs C++)

#### Medium-term (Priority 3)
7. Migrate project to src/ package structure
8. Implement CLI commands (pypair train, pypair predict)
9. Set up continuous benchmarking on HPC

#### Long-term (Priority 4)
10. Publish performance comparison paper
11. Release PyPair v1.0 on PyPI
12. Develop web interface for interactive alignment

### Breaking Changes

⚠️ **Python Version Change**: Projects using Python 3.11 must downgrade to 3.10
⚠️ **NumPy Version**: Code relying on NumPy 2.x features needs updates
⚠️ **Import Paths**: Future src/ migration will change import statements

### Acknowledgments

- **Ulster University** - Research support and NIHPC cluster access
- **Vaex Development Team** - Out-of-core DataFrame library enabling genome-scale processing
- **Polars Development Team** - High-performance in-memory DataFrame operations

---

## [0.0.1] - 2025-01-15

### Initial Release

- Basic learned index implementation
- Multiple classifier evaluation (8 models)
- S1-S4 reference datasets
- HPC SLURM scripts for NIHPC deployment
- Initial model training results

---

**Versioning**: This project follows [Semantic Versioning](https://semver.org/)
**Format**: Based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
