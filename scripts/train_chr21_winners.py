#!/usr/bin/env python3
"""
Train winning models on chr21 genome-scale dataset (80M k-mers)

Based on S1-S4 validation results, train the top 3 models:
1. Nearest Centroid (99.88% F1, 0.12% gap)
2. KNN (99.88% F1, 0.12% gap)
3. Gaussian NB (95.51% F1, 1.95% gap)

This provides realistic batch profiling (50k-500k) and non-inferiority testing.

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import os
import sys
import time
import joblib
import vaex
import pandas as pd
import numpy as np
from pathlib import Path
from memory_profiler import memory_usage
from sklearn.metrics import f1_score, accuracy_score
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier

# Add project root
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from scripts.model_fixes_and_limitations import BatchedNearestCentroid
from scripts.noninferiority_testing import (
    comprehensive_noninferiority_analysis,
    generate_noninferiority_report
)

# Configuration
CHR21_DATA = Path("data/reference_genomes/chr21_kmers.hdf5")
BWA_SEC_PER_SEED = 0.003
BATCH_SIZES = [50000, 100000, 200000, 500000]


def load_chr21_data(train_fraction=0.80, sample_size=None):
    """
    Load chr21 data with train/test split.

    Args:
        train_fraction: Fraction for training (default 0.80)
        sample_size: Optional - limit dataset size for testing

    Returns:
        X_train, y_train, X_test, y_test as pandas DataFrames
    """
    print("="*70)
    print("Loading chr21 k-mer data...")
    print("="*70)

    # Load HDF5
    df = vaex.open(str(CHR21_DATA))
    total_rows = len(df)

    print(f"Total k-mers: {total_rows:,}")
    print(f"Columns: {list(df.columns)}")

    # Sample if requested
    if sample_size and sample_size < total_rows:
        print(f"\nSampling {sample_size:,} k-mers for faster training...")
        df = df.sample(n=sample_size, random_state=42)

    # Feature columns (exclude metadata)
    feature_cols = [col for col in df.columns
                    if col not in ['chromosome', 'position', 'kmer', 'strand']]

    print(f"\nFeature columns: {feature_cols}")
    print(f"Target: position (genomic location)")

    # Convert to pandas for sklearn compatibility
    print("\nConverting to pandas (may take a moment)...")
    df_pandas = df.to_pandas_df()

    # Prepare features and target
    X = df_pandas[feature_cols]
    y = df_pandas['position']

    # Train/test split
    split_idx = int(len(X) * train_fraction)
    X_train = X.iloc[:split_idx]
    y_train = y.iloc[:split_idx]
    X_test = X.iloc[split_idx:]
    y_test = y.iloc[split_idx:]

    print(f"\nTrain set: {len(X_train):,} samples")
    print(f"Test set: {len(X_test):,} samples")
    print(f"Unique positions: {y.nunique():,}")
    print(f"Features: {X.shape[1]}")

    return X_train, y_train, X_test, y_test


def profile_batch_inference(model, X_test, batch_sizes=BATCH_SIZES):
    """
    Profile inference at multiple batch sizes with memory tracking.

    Returns:
        DataFrame with batch profiling results
    """
    results = []

    for batch_size in batch_sizes:
        if batch_size > len(X_test):
            print(f"  Batch {batch_size:>7,}: SKIPPED (exceeds test set)")
            continue

        X_batch = X_test.iloc[:batch_size]

        # Warm-up
        try:
            _ = model.predict(X_batch.iloc[:min(100, len(X_batch))])
        except:
            pass

        # Measure memory
        def predict_wrapper():
            return model.predict(X_batch)

        try:
            mem_usage = memory_usage((predict_wrapper,), interval=0.1, max_iterations=1)
            mem_mb = max(mem_usage) if mem_usage else 0
        except Exception as e:
            print(f"  Memory profiling failed: {e}")
            mem_mb = 0

        # Timed runs (3 replicates)
        times = []
        for rep in range(3):
            start = time.time()
            try:
                _ = model.predict(X_batch)
                elapsed = time.time() - start
                times.append(elapsed)
            except Exception as e:
                print(f"  Batch {batch_size}: ERROR - {e}")
                break

        if times:
            avg_time = np.mean(times)
            std_time = np.std(times)
            time_per_seed = avg_time / batch_size
            speedup = BWA_SEC_PER_SEED / time_per_seed if time_per_seed > 0 else 0

            results.append({
                'batch_size': batch_size,
                'total_time_sec': avg_time,
                'time_std_sec': std_time,
                'time_per_seed_sec': time_per_seed,
                'time_per_seed_ms': time_per_seed * 1000,
                'seeds_per_sec': batch_size / avg_time if avg_time > 0 else 0,
                'memory_mb': mem_mb,
                'memory_per_seed_kb': (mem_mb * 1024) / batch_size if batch_size > 0 else 0,
                'vs_bwa_speedup': speedup,
                'vs_bwa_status': 'FASTER' if speedup > 1 else 'SLOWER'
            })

            print(f"  Batch {batch_size:>7,}: {time_per_seed*1000:.6f} ms/seed "
                  f"({speedup:>6.2f}x vs BWA) | Mem: {mem_mb:.2f} MB")

    return pd.DataFrame(results)


def train_and_profile_model(model_name, model, X_train, y_train, X_test, y_test):
    """Train model and profile inference performance."""
    print(f"\n{'='*70}")
    print(f"{model_name} on chr21")
    print(f"{'='*70}")

    results = {'model_name': model_name}

    # Train
    print(f"\n[1/4] Training on {len(X_train):,} samples...")
    start = time.time()
    try:
        model.fit(X_train, y_train)
        train_time = time.time() - start
        results['train_time_sec'] = train_time
        results['train_status'] = 'success'
        print(f"  [OK] Training completed in {train_time:.2f}s")
    except Exception as e:
        print(f"  [ERROR] Training failed: {e}")
        results['train_status'] = f'failed: {e}'
        return results, None

    # Evaluate on train
    print(f"\n[2/4] Evaluating on training set...")
    try:
        y_train_pred = model.predict(X_train[:10000])  # Sample for speed
        train_f1 = f1_score(y_train[:10000], y_train_pred, average='weighted', zero_division=0)
        results['train_f1'] = train_f1
        print(f"  Train F1: {train_f1:.4f} (sampled 10k)")
    except Exception as e:
        print(f"  WARNING: {e}")
        results['train_f1'] = None

    # Evaluate on test
    print(f"\n[3/4] Evaluating on test set ({len(X_test):,} samples)...")
    try:
        y_test_pred = model.predict(X_test)
        test_f1 = f1_score(y_test, y_test_pred, average='weighted', zero_division=0)
        test_acc = accuracy_score(y_test, y_test_pred)
        results['test_f1'] = test_f1
        results['test_accuracy'] = test_acc

        if results.get('train_f1'):
            results['overfitting_gap'] = results['train_f1'] - test_f1
            results['overfitting_gap_pct'] = results['overfitting_gap'] * 100

        print(f"  Test F1: {test_f1:.4f} ({test_f1*100:.2f}%)")
        print(f"  Test Accuracy: {test_acc:.4f}")
        if results.get('overfitting_gap'):
            print(f"  Overfitting Gap: {results['overfitting_gap_pct']:.2f}%")
    except Exception as e:
        print(f"  ERROR: {e}")
        results['test_f1'] = None
        return results, model

    # Batch profiling
    print(f"\n[4/4] Profiling batch inference (50k-500k)...")
    print(f"  BWA benchmark: {BWA_SEC_PER_SEED*1000:.3f} ms/seed")
    batch_results = profile_batch_inference(model, X_test, BATCH_SIZES)
    results['batch_profile'] = batch_results

    # Save model
    model_dir = Path("models/model_weights")
    model_dir.mkdir(parents=True, exist_ok=True)
    model_file = model_dir / f"{model_name.lower().replace(' ', '_')}_chr21.joblib"

    try:
        joblib.dump(model, model_file)
        results['model_path'] = str(model_file)
        print(f"\n  [OK] Model saved: {model_file}")
    except Exception as e:
        print(f"  WARNING: Failed to save model: {e}")

    return results, model


def main():
    print("="*70)
    print("PyPair: chr21 Genome-Scale Training & Profiling")
    print("="*70)
    print("\nWinning Models (from S1-S4 validation):")
    print("  1. Nearest Centroid: 99.88% F1, 0.12% gap")
    print("  2. KNN: 99.88% F1, 0.12% gap")
    print("  3. Gaussian NB: 95.51% F1, 1.95% gap")
    print("\nObjective:")
    print("  - Train on chr21 (80M k-mers)")
    print("  - Realistic batch profiling (50k-500k)")
    print("  - Non-inferiority testing vs BWA")
    print("="*70)

    # Load data (use full dataset or sample for testing)
    SAMPLE_SIZE = None  # Set to e.g., 1000000 for faster testing
    X_train, y_train, X_test, y_test = load_chr21_data(sample_size=SAMPLE_SIZE)

    # Define models
    models = {
        'Nearest Centroid': BatchedNearestCentroid(batch_size=1000),
        'KNN': KNeighborsClassifier(n_neighbors=5, weights='distance', n_jobs=-1),
        'Gaussian NB': GaussianNB()
    }

    # Train and profile each model
    all_results = []
    batch_profiles = {}

    for model_name, model in models.items():
        results, trained_model = train_and_profile_model(
            model_name, model, X_train, y_train, X_test, y_test
        )
        all_results.append(results)

        if 'batch_profile' in results and isinstance(results['batch_profile'], pd.DataFrame):
            batch_profiles[model_name] = results['batch_profile']

    # Save results
    print(f"\n{'='*70}")
    print("SAVING RESULTS")
    print(f"{'='*70}")

    # Remove batch_profile for CSV export
    results_for_csv = []
    for r in all_results:
        r_copy = r.copy()
        r_copy.pop('batch_profile', None)
        results_for_csv.append(r_copy)

    results_df = pd.DataFrame(results_for_csv)
    results_df.to_csv("results/chr21_model_results.csv", index=False)
    print(f"Main results saved: results/chr21_model_results.csv")

    # Save batch profiles
    if batch_profiles:
        batch_data = []
        for model_name, batch_df in batch_profiles.items():
            batch_df['model_name'] = model_name
            batch_data.append(batch_df)

        combined_batch = pd.concat(batch_data, ignore_index=True)
        combined_batch.to_csv("results/chr21_batch_profiling.csv", index=False)
        print(f"Batch profiling saved: results/chr21_batch_profiling.csv")

    # Summary
    print(f"\n{'='*70}")
    print("CHR21 VALIDATION SUMMARY")
    print(f"{'='*70}")

    print(f"\n{'Model':<20} {'Test F1':<10} {'Gap %':<10} {'Status'}")
    print("-"*60)
    for _, row in results_df.iterrows():
        gap_pct = row.get('overfitting_gap_pct', 0)
        status = row.get('train_status', 'unknown')
        test_f1 = row.get('test_f1', 0)
        print(f"{row['model_name']:<20} {test_f1:<10.4f} {gap_pct:<10.2f} {status}")

    # Non-inferiority testing
    print(f"\n{'='*70}")
    print("NON-INFERIORITY TESTING vs BWA")
    print(f"{'='*70}")

    for model_name, batch_df in batch_profiles.items():
        print(f"\n[{model_name}]")
        try:
            ni_analysis = comprehensive_noninferiority_analysis(
                batch_df, model_name, bwa_time=BWA_SEC_PER_SEED
            )

            # Generate report
            report_file = f"results/noninferiority_{model_name.lower().replace(' ', '_')}_chr21.txt"
            report = generate_noninferiority_report(ni_analysis, report_file)

            decision = ni_analysis['overall_decision']
            print(f"  Non-Inferior: {'YES' if decision['non_inferior'] else 'NO'}")
            print(f"  Superior: {'YES' if decision['superior'] else 'NO'}")
            print(f"  Clinical Readiness: {decision['clinical_readiness']}")
            print(f"  Report: {report_file}")
        except Exception as e:
            print(f"  ERROR: {e}")

    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*70}")
    print("\nNext steps:")
    print("  1. Review chr21 results vs S1-S4 (consistency check)")
    print("  2. Examine batch scaling patterns (50k-500k)")
    print("  3. Review non-inferiority reports for clinical deployment")
    print("  4. Compare inference times vs BWA benchmark")


if __name__ == "__main__":
    main()
