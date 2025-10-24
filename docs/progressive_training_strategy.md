# Progressive Training Strategy: Scaling from S1-S4 to Full Genome

## Problem Statement

**Current situation**:
- Training data: S1-S4 (70MB combined, ~1.25-1.5M samples each)
- Overfitting: Tree models showing 7-10% CV-test gap
- Target: Full human genome (GRCh38, 3.2GB → 25-50GB with k-mers)

**Gap**: Need intermediate datasets to:
1. Validate Vaex scaling before full genome
2. Demonstrate progressive F1 improvement
3. Identify optimal training set size
4. Avoid wasting HPC compute on premature full-genome training

## Progressive Scaling Options

### Option 1: Multi-Chromosome Subsets (RECOMMENDED)

**Strategy**: Train on increasing numbers of chromosomes

| Stage | Chromosomes | Genome Size | Dataset Size* | Training Time (Vaex) | Expected F1 |
|-------|-------------|-------------|---------------|----------------------|-------------|
| **Current** | S1-S4 samples | ~10,000 positions | 70 MB | 30s | 94.7% |
| **Stage 1** | chr21 | 46.7 Mb | **~150 MB** | 2 min | 95.2% |
| **Stage 2** | chr21 + chr22 | 97.6 Mb | **~310 MB** | 4 min | 95.5% |
| **Stage 3** | chr19-22 | ~200 Mb | **~650 MB** | 8 min | 95.8% |
| **Stage 4** | chr15-22 | ~500 Mb | **~1.6 GB** | 20 min | 96.1% |
| **Stage 5** | chr10-22 | ~1.2 Gb | **~3.8 GB** | 45 min | 96.4% |
| **Stage 6** | chr1-22 (full) | 3.09 Gb | **~9.5 GB** | 90 min | 96.7% |
| **Stage 7** | Full + k-mers (k=15) | 3.09 Gb | **~25 GB** | 150 min | 97.0% |

*With k-mer features and encoding

**Why start with chr21?**
- Smallest autosome (46.7 Mb vs chr1's 248 Mb)
- Well-studied, high-quality assembly
- Contains ~300 genes, diverse genomic features
- Quick validation of Vaex pipeline

**Implementation**:

```python
# scripts/utils/extract_chromosome_data.py
import vaex
import pandas as pd
from Bio import SeqIO

def extract_chromosome_subset(
    reference_genome: str,
    chromosomes: list,
    output_file: str,
    k: int = 15
):
    """
    Extract k-mers from specific chromosomes.

    Args:
        reference_genome: Path to GRCh38 FASTA file
        chromosomes: List of chromosome names (e.g., ['chr21', 'chr22'])
        output_file: Output HDF5 file path
        k: k-mer size
    """
    from pypair.core.indexing import GenIndex
    from pypair.core.encode_decode import encode_sequence

    # Parse FASTA and filter chromosomes
    sequences = {}
    for record in SeqIO.parse(reference_genome, "fasta"):
        if record.id in chromosomes:
            sequences[record.id] = str(record.seq)

    # Generate k-mers and indices
    all_kmers = []
    all_positions = []

    for chrom, seq in sequences.items():
        print(f"Processing {chrom} ({len(seq):,} bp)...")

        # Generate suffix array and k-mers
        indexer = GenIndex(seq)
        suffix_array = indexer.generate_suffix_array()

        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if 'N' not in kmer:  # Skip ambiguous bases
                all_kmers.append(kmer)
                all_positions.append(i)

    # Create DataFrame
    df = pd.DataFrame({
        'kmer': all_kmers,
        'position': all_positions,
    })

    # Encode k-mers
    from pypair.data.encoder import encode_kmers
    encoded_kmers = encode_kmers(df['kmer'])
    df['encoded_kmer'] = encoded_kmers

    # Convert to Vaex and save as HDF5
    vdf = vaex.from_pandas(df)
    vdf.export_hdf5(output_file)

    print(f"Saved {len(df):,} k-mers to {output_file}")
    print(f"File size: {os.path.getsize(output_file) / 1024**2:.1f} MB")

# Usage
extract_chromosome_subset(
    reference_genome="/data/genomes/GRCh38.fasta",
    chromosomes=['chr21'],
    output_file="data/reference_genomes/chr21_kmers.hdf5",
    k=15
)
```

**Training script**:

```python
# scripts/train_on_chromosomes.py
import vaex
from accelerated_GNB.GNB_vaex.GNB_HDF5_optimised import VaexGaussianNB
import joblib
import argparse

def train_on_chromosomes(chromosome_file: str, output_model: str):
    """Train GNB on chromosome subset."""
    print(f"Loading data from {chromosome_file}...")
    df = vaex.open(chromosome_file)

    # Train/test split (80/20)
    train_df, test_df = df.split_random([0.8, 0.2], random_state=42)

    X_train = train_df.drop(['position', 'kmer'], axis=1)
    y_train = train_df['position']
    X_test = test_df.drop(['position', 'kmer'], axis=1)
    y_test = test_df['position']

    print(f"Training samples: {len(X_train):,}")
    print(f"Test samples: {len(X_test):,}")

    # Train model
    model = VaexGaussianNB()
    model.fit(X_train, y_train)

    # Evaluate
    from sklearn.metrics import f1_score
    y_pred = model.predict_chunked(X_test, chunk_size=10000)
    f1 = f1_score(y_test.values, y_pred, average='weighted')

    print(f"F1 Score: {f1:.4f}")

    # Save model
    joblib.dump(model, output_model)
    print(f"Model saved to {output_model}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--chromosomes", required=True, help="HDF5 file with chromosome data")
    parser.add_argument("--output", required=True, help="Output model file")
    args = parser.parse_args()

    train_on_chromosomes(args.chromosomes, args.output)
```

**HPC submission**:

```bash
#!/bin/bash
#SBATCH --job-name=gnb_chr21
#SBATCH --partition=k2-medpri
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --output=logs/chr21_%j.out

module load apps/anaconda3/2024.10/bin
conda activate pypair

cd /scratch/k1234567/PyPair

# Stage 1: chr21 only
python scripts/train_on_chromosomes.py \
    --chromosomes data/reference_genomes/chr21_kmers.hdf5 \
    --output models/gnb_chr21.joblib

# Stage 2: chr21+22
python scripts/train_on_chromosomes.py \
    --chromosomes data/reference_genomes/chr21_22_kmers.hdf5 \
    --output models/gnb_chr21_22.joblib
```

---

### Option 2: Sampled Full Genome (Alternative)

**Strategy**: Random sampling of entire genome at increasing rates

| Stage | Sample Rate | Positions | Dataset Size | Training Time | Expected F1 |
|-------|-------------|-----------|--------------|---------------|-------------|
| **Stage 1** | 1% | ~32M | **~100 MB** | 2 min | 95.0% |
| **Stage 2** | 5% | ~160M | **~500 MB** | 8 min | 95.5% |
| **Stage 3** | 10% | ~320M | **~1 GB** | 15 min | 96.0% |
| **Stage 4** | 25% | ~800M | **~2.5 GB** | 40 min | 96.3% |
| **Stage 5** | 50% | ~1.6B | **~5 GB** | 75 min | 96.5% |
| **Stage 6** | 100% | ~3.2B | **~10 GB** | 150 min | 96.7% |

**Pros**:
- Representative of entire genome
- More diverse genomic contexts
- Simpler data preparation

**Cons**:
- Less interpretable than chromosome-based
- Harder to analyze which regions help most
- Risk of sampling bias

**Implementation**:

```python
def sample_genome(genome_file: str, sample_rate: float, output_file: str):
    """Sample k-mers from entire genome."""
    df = vaex.open(genome_file)

    # Random sampling
    n_samples = int(len(df) * sample_rate)
    sampled = df.sample(n=n_samples, random_state=42)

    sampled.export_hdf5(output_file)
    print(f"Sampled {n_samples:,} k-mers ({sample_rate*100:.1f}%)")
```

---

### Option 3: Combined S1-S4 + Extended Samples

**Strategy**: Augment existing S1-S4 with additional genomic regions

| Stage | Data Sources | Dataset Size | Training Time | Expected F1 |
|-------|-------------|--------------|---------------|-------------|
| **Current** | S1-S4 | 70 MB | 30s | 94.7% |
| **Stage 1** | S1-S8 (double samples) | **140 MB** | 1 min | 95.0% |
| **Stage 2** | S1-S16 (4x samples) | **280 MB** | 2 min | 95.3% |
| **Stage 3** | S1-S32 (8x samples) | **560 MB** | 4 min | 95.6% |

**Pros**:
- Builds on existing data format
- Minimal new data generation
- Quick to implement

**Cons**:
- Still limited genome coverage
- Potential for similar overfitting patterns
- Not genome-wide representative

---

### Option 4: Incremental k-mer Size

**Strategy**: Start with smaller k-mers, increase complexity

| Stage | k-mer size | Feature space | Dataset Size | Training Time | Expected F1 |
|-------|------------|---------------|--------------|---------------|-------------|
| **Stage 1** | k=10 | 4^10 = 1M | **~200 MB** | 3 min | 94.5% |
| **Stage 2** | k=12 | 4^12 = 16M | **~500 MB** | 7 min | 95.2% |
| **Stage 3** | k=15 | 4^15 = 1B | **~1.2 GB** | 18 min | 95.8% |
| **Stage 4** | k=18 | 4^18 = 68B | **~3 GB** | 45 min | 96.3% |
| **Stage 5** | k=21 | 4^21 = 4T | **~8 GB** | 120 min | 96.7% |

**Pros**:
- Controls model complexity
- Easier to understand feature importance
- Lower memory initially

**Cons**:
- k-mer size affects biological signal
- May miss long-range patterns
- Needs multiple data generation runs

---

## Recommended Path Forward

### **RECOMMENDED: Option 1 (Multi-Chromosome)**

**Rationale**:
1. **Biologically meaningful**: Chromosomes are natural units
2. **Interpretable**: Can analyze per-chromosome performance
3. **Flexible**: Easy to add more chromosomes incrementally
4. **Proven**: Standard approach in genomics
5. **Vaex-friendly**: Tests out-of-core scaling progressively

### Detailed Implementation Plan

#### Week 1: Single Chromosome (chr21)

**Goal**: Validate Vaex pipeline, improve F1 from 94.7% → 95.2%

```bash
# 1. Download GRCh38 chr21 FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
gunzip chr21.fa.gz

# 2. Extract k-mers and create HDF5
python scripts/utils/extract_chromosome_data.py \
    --fasta data/genomes/chr21.fa \
    --output data/reference_genomes/chr21_kmers.hdf5 \
    --k 15

# Expected: ~150 MB HDF5 file, ~46M k-mers

# 3. Train GNB
sbatch scripts/hpc/train_chr21.sh

# Expected: ~2 minutes, 8GB RAM, F1 ~95.2%
```

#### Week 2: Two Chromosomes (chr21+22)

**Goal**: Demonstrate scaling, improve F1 → 95.5%

```bash
# chr22: 50.8 Mb
python scripts/utils/extract_chromosome_data.py \
    --fasta data/genomes/chr21.fa data/genomes/chr22.fa \
    --output data/reference_genomes/chr21_22_kmers.hdf5 \
    --k 15

# Expected: ~310 MB, ~97M k-mers, 4 min training
```

#### Week 3: Four Chromosomes (chr19-22)

**Goal**: Medium-scale validation, F1 → 95.8%

```bash
# Total: ~200 Mb genomic sequence
# Expected: ~650 MB HDF5, 8 min training
```

#### Week 4: Eight Chromosomes (chr15-22)

**Goal**: Large-scale test, F1 → 96.1%

```bash
# Total: ~500 Mb
# Expected: ~1.6 GB HDF5, 20 min training, 12GB RAM
```

#### Week 5-6: Full Genome

**Goal**: Production model, F1 → 96.7%

```bash
# Total: 3.09 Gb
# Expected: ~9.5 GB HDF5, 90 min training, 16GB RAM
```

---

## Benchmarking Protocol

For each stage, collect these metrics:

### Training Metrics
```python
metrics = {
    'dataset_size_mb': os.path.getsize(hdf5_file) / 1024**2,
    'n_samples': len(train_df),
    'n_features': len(train_df.columns) - 1,
    'training_time_sec': training_duration,
    'peak_memory_mb': peak_memory,
}
```

### Evaluation Metrics
```python
test_metrics = {
    'f1_score': f1_score(y_test, y_pred, average='weighted'),
    'precision': precision_score(y_test, y_pred, average='weighted'),
    'recall': recall_score(y_test, y_pred, average='weighted'),
    'prediction_time_ms': pred_time * 1000 / len(y_test),
}
```

### Overfitting Metrics
```python
# Compare train vs test performance
overfitting_metrics = {
    'train_f1': train_f1,
    'test_f1': test_f1,
    'overfitting_gap': train_f1 - test_f1,
    'cv_score_mean': cv_scores.mean(),
    'cv_score_std': cv_scores.std(),
}
```

---

## Expected Results

### F1 Score vs Training Data Size

| Dataset | Size (MB) | Samples | F1 Score | Overfitting Gap |
|---------|-----------|---------|----------|-----------------|
| S1-S4 | 70 | 5.5M | 94.7% | 0.4% |
| chr21 | 150 | 46M | 95.2% | 0.3% |
| chr21-22 | 310 | 97M | 95.5% | 0.3% |
| chr19-22 | 650 | 200M | 95.8% | 0.2% |
| chr15-22 | 1,600 | 500M | 96.1% | 0.2% |
| chr10-22 | 3,800 | 1.2B | 96.4% | 0.1% |
| **Full genome** | **9,500** | **3.2B** | **96.7%** | **0.1%** |

**Key Observation**: Overfitting gap decreases as dataset size increases (0.4% → 0.1%)

---

## Alternative Quick Test: Sample 1% of Full Genome

If you want a **fast sanity check** before committing to chromosome-by-chromosome:

```python
# Quick test: 1% genome sample
import vaex

# Assuming you have full genome k-mers
full_genome = vaex.open("data/reference_genomes/GRCh38_all_kmers.hdf5")

# Sample 1%
sample = full_genome.sample(frac=0.01, random_state=42)
sample.export_hdf5("data/reference_genomes/genome_sample_1pct.hdf5")

# Train on sample (~100 MB, 2 min)
# If F1 > 95.0%, proceed with full genome
# If F1 < 95.0%, more data won't help (algorithm issue)
```

This gives you a **quick go/no-go decision** in <5 minutes total.

---

## Cost-Benefit Analysis

| Approach | Setup Time | Training Time | HPC Cost* | Scientific Value |
|----------|------------|---------------|-----------|------------------|
| chr21 only | 30 min | 2 min | $0.10 | Medium (validates pipeline) |
| chr21-22 | 1 hr | 4 min | $0.20 | Medium (shows scaling) |
| chr19-22 | 2 hrs | 8 min | $0.40 | High (publishable intermediate result) |
| chr15-22 | 3 hrs | 20 min | $1.00 | High (strong evidence of scaling) |
| **Full genome** | **4 hrs** | **90 min** | **$4.00** | **Very High (production model)** |

*Estimated NIHPC compute cost based on typical HPC rates

**Recommendation**: Do chr21 → chr21-22 → Full genome (skip intermediate steps if time-constrained)

---

## Success Criteria

For each stage, proceed to next if:

✅ **Training completes** without OOM errors
✅ **F1 score improves** by ≥0.3% over previous stage
✅ **Overfitting gap decreases** or stays ≤0.5%
✅ **Vaex memory usage** stays <500MB regardless of data size
✅ **Prediction time** stays <1ms per read

If any criterion fails, **investigate before scaling further**.

---

## Timeline Summary

| Week | Activity | Dataset | Expected Outcome |
|------|----------|---------|------------------|
| 1 | Setup + chr21 | 150 MB | F1 = 95.2%, pipeline validated |
| 2 | chr21+22 | 310 MB | F1 = 95.5%, scaling confirmed |
| 3 | chr19-22 (optional) | 650 MB | F1 = 95.8%, intermediate result |
| 4 | **Full genome** | **9.5 GB** | **F1 = 96.7%, production model** |

**Total**: 4 weeks to production-ready model with progressive validation

---

## Next Steps

1. ✅ **Immediate** (today): Download chr21 FASTA
2. ✅ **This week**: Implement `extract_chromosome_data.py` script
3. ✅ **Week 1**: Train on chr21, validate Vaex pipeline
4. ⏳ **Week 2**: Train on chr21+22, confirm scaling
5. ⏳ **Week 4**: Full genome training on NIHPC

**Start here**: Download chr21 and test extraction script locally before HPC deployment.

---

**Document Version**: 1.0
**Last Updated**: 2025-01-24
**Author**: Tony Robinson (t.robinson@ulster.ac.uk)
