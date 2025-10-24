# Project Reorganization Plan

## Current Status

PyPair currently has a **flat root directory structure** with Python modules scattered in the project root. This document outlines the reorganization to a proper Python package structure.

## Current Structure (Problematic)

```
PyPair/
├── encode_decode.py               # Root-level modules (BAD)
├── gen_index.py
├── gen_training_data.py
├── io_processing.py
├── input_encoder.py
├── predict_seeds.py
├── gaussianNaiveBayes_classifier.py
├── random_forest.py
├── (... 10+ more classifier files)
├── accelerated_GNB/               # Good subdirectory
├── models/                        # Good subdirectory
├── data/                          # Good subdirectory
└── research/                      # Good subdirectory
```

**Problems**:
- ❌ No `src/` directory - violates modern Python packaging standards
- ❌ Modules not importable as `pypair.module`
- ❌ Mixing source code with data/config files
- ❌ No clear separation of concerns
- ❌ CLI entry points not defined
- ❌ Difficult to install as package

## Target Structure (Best Practice)

```
PyPair/
├── README.md                      # Main documentation
├── pyproject.toml                 # Package configuration (UV/pip)
├── .python-version                # Python 3.10 requirement
├── .gitignore                     # Updated to exclude large files
│
├── src/
│   └── pypair/                    # Main package
│       ├── __init__.py
│       ├── __version__.py
│       │
│       ├── core/                  # Core utilities
│       │   ├── __init__.py
│       │   ├── encode_decode.py   # DNA encoding/decoding
│       │   ├── io_processing.py   # FASTA/FASTQ I/O
│       │   └── indexing.py        # Suffix array/FM-index (gen_index.py)
│       │
│       ├── data/                  # Data processing
│       │   ├── __init__.py
│       │   ├── generator.py       # Training data generation
│       │   └── encoder.py         # Input encoding
│       │
│       ├── models/                # ML models
│       │   ├── __init__.py
│       │   ├── base.py            # Base classifier interface
│       │   ├── gaussian_nb.py     # Standard GNB (sklearn wrapper)
│       │   ├── random_forest.py
│       │   ├── lgbm.py
│       │   ├── xgboost.py
│       │   └── benchmark.py       # Benchmarking utilities
│       │
│       ├── accelerated/           # High-performance implementations
│       │   ├── __init__.py
│       │   ├── polars_gnb.py      # Polars GNB
│       │   ├── vaex_gnb.py        # Vaex GNB
│       │   └── cpp/               # C++/CUDA (external)
│       │       ├── GNB.cpp
│       │       ├── GNB_cuda.cu
│       │       └── CMakeLists.txt
│       │
│       ├── training/              # Training pipelines
│       │   ├── __init__.py
│       │   ├── single_model.py    # Train individual model
│       │   └── multi_model.py     # Train all models (multimodel_training.py)
│       │
│       ├── alignment/             # Alignment algorithms
│       │   ├── __init__.py
│       │   └── predict_seeds.py   # Seed prediction & alignment
│       │
│       └── cli/                   # Command-line interface
│           ├── __init__.py
│           ├── main.py            # Main CLI entry point
│           ├── train.py           # Training CLI
│           ├── predict.py         # Prediction CLI
│           └── benchmark.py       # Benchmarking CLI
│
├── scripts/                       # HPC/SLURM scripts
│   ├── hpc/
│   │   ├── run_all.sh
│   │   ├── data_prep.sh
│   │   ├── gaussianNaiveBayes_classifier.sh
│   │   ├── random_forest.sh
│   │   └── (other SLURM scripts)
│   └── utils/
│       └── convert_to_hdf5.py
│
├── data/                          # Data directory (gitignored)
│   ├── reference_samples/
│   ├── reference_genomes/
│   └── README.md                  # Data description
│
├── models/                        # Model artifacts (gitignored)
│   ├── model_weights/
│   └── README.md                  # Model versioning info
│
├── tests/                         # Unit tests
│   ├── __init__.py
│   ├── test_encode_decode.py
│   ├── test_indexing.py
│   ├── test_models.py
│   └── test_alignment.py
│
├── research/                      # Research notebooks
│   ├── notebooks/
│   │   ├── multi_model_training.ipynb
│   │   ├── training_data_generation_v2.ipynb
│   │   └── model_analysis.ipynb
│   └── results/                   # Experimental results
│
├── docs/                          # Documentation
│   ├── architecture.md
│   ├── gnb_implementations.md
│   ├── overfitting_analysis.md
│   ├── benchmarking.md
│   ├── hpc_deployment.md
│   └── api/                       # API documentation
│
├── Validation/                    # Validation datasets (gitignored)
│   └── S*_reference.csv
│
└── logs/                          # Log files (gitignored)
```

## Migration Steps

### Phase 1: Create Package Structure (IMMEDIATE)

```bash
# Create src directory structure
mkdir -p src/pypair/{core,data,models,accelerated,training,alignment,cli}

# Create __init__.py files
touch src/pypair/__init__.py
touch src/pypair/{core,data,models,accelerated,training,alignment,cli}/__init__.py

# Create version file
echo "__version__ = '0.1.0'" > src/pypair/__version__.py
```

### Phase 2: Move Core Modules

```bash
# Move encoding/decoding
mv encode_decode.py src/pypair/core/

# Move I/O processing
mv io_processing.py src/pypair/core/

# Move indexing
mv gen_index.py src/pypair/core/indexing.py

# Move data processing
mv gen_training_data.py src/pypair/data/generator.py
mv input_encoder.py src/pypair/data/encoder.py
```

### Phase 3: Move Model Implementations

```bash
# Move classifiers to models/
mv gaussianNaiveBayes_classifier.py src/pypair/models/gaussian_nb.py
mv random_forest.py src/pypair/models/
mv LGBM_classifier.py src/pypair/models/lgbm.py
mv XGB_classifier.py src/pypair/models/xgboost.py
mv kNeighbours_classifier.py src/pypair/models/knn.py
mv nearest_centroid.py src/pypair/models/
mv decisiontrees_classifier.py src/pypair/models/decision_tree.py
mv extratrees_classifier.py src/pypair/models/extra_trees.py

# Move accelerated implementations
mv accelerated_GNB/GNB_polars/GNB_polars.py src/pypair/accelerated/polars_gnb.py
mv accelerated_GNB/GNB_vaex/GNB_HDF5_optimised.py src/pypair/accelerated/vaex_gnb.py
mv accelerated_GNB/Cpp/ src/pypair/accelerated/cpp/
```

### Phase 4: Move Training & Alignment

```bash
# Move training pipeline
mv multimodel_training.py src/pypair/training/multi_model.py

# Move prediction/alignment
mv predict_seeds.py src/pypair/alignment/
```

### Phase 5: Create CLI Structure

Create **src/pypair/cli/main.py**:

```python
#!/usr/bin/env python3
"""PyPair command-line interface."""
import argparse
from pypair import __version__


def main():
    parser = argparse.ArgumentParser(
        description="PyPair: Genomic Sequence Alignment",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"PyPair {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Train command
    train_parser = subparsers.add_parser("train", help="Train models")
    train_parser.add_argument("--model", required=True, help="Model type")
    train_parser.add_argument("--dataset", required=True, help="Dataset name")

    # Predict command
    predict_parser = subparsers.add_parser("predict", help="Predict alignments")
    predict_parser.add_argument("--input", required=True, help="Input FASTQ file")
    predict_parser.add_argument("--output", required=True, help="Output BAM file")

    # Benchmark command
    bench_parser = subparsers.add_parser("benchmark", help="Run benchmarks")
    bench_parser.add_argument("--implementations", nargs="+", help="Implementations to test")

    args = parser.parse_args()

    if args.command == "train":
        from pypair.cli.train import train_command
        train_command(args)
    elif args.command == "predict":
        from pypair.cli.predict import predict_command
        predict_command(args)
    elif args.command == "benchmark":
        from pypair.cli.benchmark import benchmark_command
        benchmark_command(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
```

### Phase 6: Organize HPC Scripts

```bash
# Create scripts directory
mkdir -p scripts/hpc scripts/utils

# Move SLURM scripts
mv *.sh scripts/hpc/

# Keep run_all.sh in root for convenience
cp scripts/hpc/run_all.sh ./
```

### Phase 7: Update Imports

**Critical**: All files must update their imports from:

```python
# OLD (flat structure)
from encode_decode import encode_sequence
from gen_index import GenIndex
from gaussianNaiveBayes_classifier import train_model

# NEW (package structure)
from pypair.core.encode_decode import encode_sequence
from pypair.core.indexing import GenIndex
from pypair.models.gaussian_nb import train_model
```

**Automated refactoring script**:

```python
# scripts/utils/refactor_imports.py
import os
import re
from pathlib import Path


def update_imports(file_path):
    """Update imports in a Python file."""
    with open(file_path, "r") as f:
        content = f.read()

    # Replace old imports with new ones
    replacements = {
        r"from encode_decode import": "from pypair.core.encode_decode import",
        r"from gen_index import": "from pypair.core.indexing import",
        r"from io_processing import": "from pypair.core.io_processing import",
        r"from gen_training_data import": "from pypair.data.generator import",
        r"from gaussianNaiveBayes_classifier import": "from pypair.models.gaussian_nb import",
        r"from predict_seeds import": "from pypair.alignment.predict_seeds import",
        # Add more as needed
    }

    for old, new in replacements.items():
        content = re.sub(old, new, content)

    with open(file_path, "w") as f:
        f.write(content)
    print(f"Updated: {file_path}")


# Run on all Python files
for py_file in Path("src/pypair").rglob("*.py"):
    update_imports(py_file)
```

### Phase 8: Update pyproject.toml

Ensure **pyproject.toml** specifies the correct package location:

```toml
[tool.setuptools.packages.find]
where = ["src"]
include = ["pypair*"]

[project.scripts]
pypair = "pypair.cli.main:main"
pypair-train = "pypair.cli.train:main"
pypair-predict = "pypair.cli.predict:main"
pypair-benchmark = "pypair.cli.benchmark:main"
```

### Phase 9: Reinstall Package

```bash
# Uninstall old version
pip uninstall pypair -y

# Install new structure
uv pip install -e .

# Verify CLI works
pypair --version
pypair train --help
pypair predict --help
```

### Phase 10: Update HPC Scripts

Update SLURM scripts to use new CLI:

```bash
# OLD
python gaussianNaiveBayes_classifier.py --dataset S1

# NEW
pypair train --model gaussian_nb --dataset S1
```

## Testing the Migration

### 1. Import Tests

```python
# tests/test_imports.py
def test_core_imports():
    from pypair.core import encode_decode, io_processing, indexing
    assert True

def test_model_imports():
    from pypair.models import gaussian_nb, random_forest
    assert True

def test_accelerated_imports():
    from pypair.accelerated import polars_gnb, vaex_gnb
    assert True

def test_cli_imports():
    from pypair.cli import main
    assert True
```

### 2. Functional Tests

```bash
# Test training
pypair train --model gaussian_nb --dataset S1

# Test prediction
pypair predict --input test.fastq --output test.bam

# Test benchmarking
pypair benchmark --implementations polars vaex
```

### 3. HPC Integration Test

```bash
# Submit test job
sbatch scripts/hpc/gaussianNaiveBayes_classifier.sh

# Verify output
cat logs/gnb_vaex_*.out
```

## Benefits of New Structure

✅ **Standard Python Package** - Follows PEP 518/621 conventions
✅ **Importable Modules** - `from pypair.core import encode_decode`
✅ **CLI Entry Points** - `pypair train`, `pypair predict`, etc.
✅ **Clean Separation** - Source code in `src/`, data in `data/`, scripts in `scripts/`
✅ **Testable** - Clear module boundaries for unit tests
✅ **Publishable** - Can be uploaded to PyPI
✅ **Maintainable** - Logical organization by functionality

## Rollback Plan

If migration causes issues:

```bash
# 1. Restore from git
git checkout HEAD -- .

# 2. Or keep backup
cp -r PyPair PyPair_backup_20250124
```

## Timeline

| Phase | Duration | Status |
|-------|----------|--------|
| 1-3: Create structure & move core | 1 hour | Pending |
| 4-5: Move models & training | 1 hour | Pending |
| 6-7: Update imports & scripts | 2 hours | Pending |
| 8-9: Package config & reinstall | 30 min | Pending |
| 10: Update HPC scripts | 1 hour | Pending |
| Testing & validation | 2 hours | Pending |
| **Total** | **~8 hours** | **Pending** |

## Recommendation

**When to Migrate?**

✅ **Do it soon** - Before adding more features
✅ **After current work** - Complete current research first
⚠️ **Not during active HPC jobs** - Wait for cluster jobs to finish

**Priority**: Medium-High

Current flat structure is functional but will become unmaintainable as project grows. Recommend completing migration within next 2 weeks.

---

**Document Version**: 1.0
**Last Updated**: 2025-01-24
**Status**: PENDING - Awaiting execution
**Author**: Tony Robinson (t.robinson@ulster.ac.uk)
