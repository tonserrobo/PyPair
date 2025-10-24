"""
PyPair: Fast Seed Mapping Profiler (No Memory Profiler)
=========================================================

Profiles trained models for seed mapping performance using lightweight profiling:
- Execution time (sec per seed) - HIGH PRIORITY for viva
- Throughput (seeds per second)
- Comparison vs BWA benchmark (0.003 sec/seed)

Usage:
    python scripts/profile_seed_mapping_fast.py

Output:
    results/seed_mapping_profile_fast.csv - Complete profiling metrics
    results/seed_mapping_summary_fast.txt - Human-readable summary
"""

import joblib
import pandas as pd
import numpy as np
import time
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
BWA_SEC_PER_SEED = 0.003  # BWA benchmark (20-year gold standard)
BATCH_SIZES = [1000, 5000, 10000, 14000]  # Include full test set
DATASETS = ['S1', 'S2', 'S3', 'S4']
RESULTS_DIR = Path("results")
MODEL_DIR = Path("models/model_weights")

# Models to profile (only successful models from validation)
MODELS_TO_PROFILE = {
    'KNN': 'knn',
    'Nearest Centroid': 'nearest_centroid',
    'Gaussian NB': 'gaussian_nb'
}


def load_dataset(dataset_name: str):
    """Load test data for a given dataset."""
    base_path = Path("data/reference_samples/cleaned")

    X_test = pd.read_csv(base_path / f"{dataset_name}_extended_encoded_X_test.csv")
    y_test = pd.read_csv(base_path / f"{dataset_name}_extended_encoded_y_test.csv")

    return X_test.values, y_test['org_location'].values


def load_model(model_name: str, dataset: str):
    """Load trained model from model_weights directory."""
    model_filename = f"{model_name}_{dataset.lower()}.joblib"
    model_path = MODEL_DIR / model_filename

    if not model_path.exists():
        raise FileNotFoundError(f"Model not found: {model_path}")

    return joblib.load(model_path)


def profile_batch_prediction(model, X_batch):
    """Profile time for a single batch prediction (5 runs for stability)."""

    # Time profiling (5 runs for better stability, take median)
    times = []
    for _ in range(5):
        start_time = time.perf_counter()  # Higher precision timer
        _ = model.predict(X_batch)
        elapsed = time.perf_counter() - start_time
        times.append(elapsed)

    # Use median to reduce impact of outliers
    median_time_sec = np.median(times)
    mean_time_sec = np.mean(times)
    std_time_sec = np.std(times)
    min_time_sec = np.min(times)
    max_time_sec = np.max(times)

    return {
        'median_time': median_time_sec,
        'mean_time': mean_time_sec,
        'std_time': std_time_sec,
        'min_time': min_time_sec,
        'max_time': max_time_sec
    }


def profile_model_on_dataset(model_name: str, dataset: str, verbose: bool = True):
    """
    Profile a single model on a dataset across multiple batch sizes.

    Returns:
        List of result dictionaries with profiling metrics
    """
    results = []

    if verbose:
        print(f"\n{'='*70}")
        print(f"Profiling: {model_name} on {dataset}")
        print(f"{'='*70}")

    # Load model and data
    try:
        model_file = MODELS_TO_PROFILE[model_name]
        model = load_model(model_file, dataset)
        X_test, y_test = load_dataset(dataset)

        if verbose:
            print(f"  Model loaded: {model_file}_{dataset.lower()}.joblib")
            print(f"  Test samples: {len(X_test):,}")
            print(f"  Features: {X_test.shape[1]}")

    except Exception as e:
        print(f"  ERROR loading model/data: {e}")
        return results

    # Profile across batch sizes
    for batch_size in BATCH_SIZES:
        if batch_size > len(X_test):
            batch_size = len(X_test)  # Use full test set
            batch_label = f"{batch_size} (full)"
        else:
            batch_label = f"{batch_size}"

        # Extract batch
        X_batch = X_test[:batch_size]

        try:
            # Profile prediction
            timing_stats = profile_batch_prediction(model, X_batch)

            # Calculate metrics (using median time for robustness)
            time_per_seed_sec = timing_stats['median_time'] / batch_size
            time_per_seed_ms = time_per_seed_sec * 1000
            time_per_seed_us = time_per_seed_sec * 1_000_000
            seeds_per_sec = batch_size / timing_stats['median_time'] if timing_stats['median_time'] > 0 else 0

            # Compare to BWA
            speedup_vs_bwa = BWA_SEC_PER_SEED / time_per_seed_sec if time_per_seed_sec > 0 else 0
            faster_than_bwa = speedup_vs_bwa > 1

            # Store results
            result = {
                'model_name': model_name,
                'dataset': dataset,
                'batch_size': batch_size,
                'median_time_sec': timing_stats['median_time'],
                'mean_time_sec': timing_stats['mean_time'],
                'std_time_sec': timing_stats['std_time'],
                'min_time_sec': timing_stats['min_time'],
                'max_time_sec': timing_stats['max_time'],
                'time_per_seed_sec': time_per_seed_sec,
                'time_per_seed_ms': time_per_seed_ms,
                'time_per_seed_us': time_per_seed_us,
                'seeds_per_sec': seeds_per_sec,
                'bwa_benchmark_sec': BWA_SEC_PER_SEED,
                'bwa_benchmark_ms': BWA_SEC_PER_SEED * 1000,
                'speedup_vs_bwa': speedup_vs_bwa,
                'faster_than_bwa': faster_than_bwa,
                'status': 'FASTER' if faster_than_bwa else 'SLOWER'
            }

            results.append(result)

            if verbose:
                print(f"  Batch {batch_label:>12}: {time_per_seed_ms:.6f} ms/seed "
                      f"({speedup_vs_bwa:>6.2f}x vs BWA) | {seeds_per_sec:,.0f} seeds/sec")

        except Exception as e:
            print(f"  Batch {batch_label:>12}: ERROR - {e}")
            results.append({
                'model_name': model_name,
                'dataset': dataset,
                'batch_size': batch_size,
                'status': f'ERROR: {str(e)}'
            })

    return results


def generate_summary_report(results_df: pd.DataFrame, output_file: Path):
    """Generate human-readable summary report."""

    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("PyPair: Seed Mapping Performance Profile (Fast)\n")
        f.write("="*80 + "\n\n")

        # Overall statistics
        successful = results_df[results_df['status'].isin(['FASTER', 'SLOWER'])]

        f.write("SUMMARY STATISTICS\n")
        f.write("-"*80 + "\n")
        f.write(f"Total profiles: {len(results_df)}\n")
        f.write(f"Successful: {len(successful)}\n")
        f.write(f"Failed: {len(results_df) - len(successful)}\n\n")

        if len(successful) > 0:
            # Best performers (by batch size = 10k for consistency)
            batch_10k = successful[successful['batch_size'] == 10000]

            if len(batch_10k) > 0:
                f.write("TOP PERFORMERS (10k seed batch)\n")
                f.write("-"*80 + "\n")

                top_performers = batch_10k.sort_values('time_per_seed_ms').head(5)

                for idx, row in top_performers.iterrows():
                    f.write(f"{row['model_name']} on {row['dataset']}:\n")
                    f.write(f"  Time/seed: {row['time_per_seed_ms']:.6f} ms ({row['time_per_seed_us']:.2f} us)\n")
                    f.write(f"  Throughput: {row['seeds_per_sec']:,.2f} seeds/sec\n")
                    f.write(f"  vs BWA: {row['speedup_vs_bwa']:.2f}x ({row['status']})\n")
                    f.write(f"  Stability: median={row['median_time_sec']:.6f}s, std={row['std_time_sec']:.6f}s\n\n")

            # BWA comparison summary
            f.write("\nBWA COMPARISON SUMMARY\n")
            f.write("-"*80 + "\n")
            f.write(f"BWA benchmark: {BWA_SEC_PER_SEED*1000:.3f} ms/seed\n\n")

            faster = successful[successful['faster_than_bwa'] == True]
            slower = successful[successful['faster_than_bwa'] == False]

            f.write(f"Faster than BWA: {len(faster)}/{len(successful)} ({len(faster)/len(successful)*100:.1f}%)\n")
            f.write(f"Slower than BWA: {len(slower)}/{len(successful)} ({len(slower)/len(successful)*100:.1f}%)\n\n")

            if len(faster) > 0:
                avg_speedup = faster['speedup_vs_bwa'].mean()
                max_speedup = faster['speedup_vs_bwa'].max()
                f.write(f"Average speedup (for faster models): {avg_speedup:.2f}x\n")
                f.write(f"Maximum speedup: {max_speedup:.2f}x\n\n")

            # Model-wise summary
            f.write("\nMODEL-WISE AVERAGE PERFORMANCE\n")
            f.write("-"*80 + "\n")

            model_summary = successful.groupby('model_name').agg({
                'time_per_seed_ms': 'mean',
                'time_per_seed_us': 'mean',
                'seeds_per_sec': 'mean',
                'speedup_vs_bwa': 'mean'
            }).sort_values('time_per_seed_ms')

            f.write(f"{'Model':<20} {'ms/seed':<12} {'us/seed':<12} {'seeds/sec':<15} {'BWA speedup'}\n")
            f.write("-"*80 + "\n")

            for model, row in model_summary.iterrows():
                f.write(f"{model:<20} {row['time_per_seed_ms']:<12.6f} "
                       f"{row['time_per_seed_us']:<12.2f} {row['seeds_per_sec']:<15,.0f} "
                       f"{row['speedup_vs_bwa']:.2f}x\n")

            # Batch size comparison
            f.write("\n\nBATCH SIZE IMPACT (Average across all models)\n")
            f.write("-"*80 + "\n")

            batch_summary = successful.groupby('batch_size').agg({
                'time_per_seed_ms': 'mean',
                'seeds_per_sec': 'mean',
                'speedup_vs_bwa': 'mean'
            }).sort_values('batch_size')

            f.write(f"{'Batch Size':<15} {'ms/seed':<12} {'seeds/sec':<15} {'BWA speedup'}\n")
            f.write("-"*80 + "\n")

            for batch, row in batch_summary.iterrows():
                f.write(f"{batch:<15,} {row['time_per_seed_ms']:<12.6f} "
                       f"{row['seeds_per_sec']:<15,.0f} {row['speedup_vs_bwa']:.2f}x\n")

    print(f"\n[OK] Summary report saved: {output_file}")


def main():
    """Main profiling workflow."""
    print("="*70)
    print("PyPair: Fast Seed Mapping Profiler")
    print("="*70)
    print(f"BWA benchmark: {BWA_SEC_PER_SEED*1000:.3f} ms/seed ({BWA_SEC_PER_SEED*1_000_000:.1f} us/seed)")
    print(f"Batch sizes: {', '.join([f'{b:,}' for b in BATCH_SIZES])}")
    print(f"Datasets: {', '.join(DATASETS)}")
    print(f"Models: {', '.join(MODELS_TO_PROFILE.keys())}")
    print(f"Timing runs: 5 per batch (median reported)")
    print("="*70)

    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Profile all models
    all_results = []

    for dataset in DATASETS:
        print(f"\n{'#'*70}")
        print(f"# Dataset: {dataset}")
        print(f"{'#'*70}")

        for model_name in MODELS_TO_PROFILE.keys():
            results = profile_model_on_dataset(model_name, dataset, verbose=True)
            all_results.extend(results)

    # Save results
    results_df = pd.DataFrame(all_results)
    output_csv = RESULTS_DIR / "seed_mapping_profile_fast.csv"
    results_df.to_csv(output_csv, index=False)
    print(f"\n{'='*70}")
    print(f"[OK] Results saved: {output_csv}")
    print(f"     Rows: {len(results_df)}")

    # Generate summary report
    output_summary = RESULTS_DIR / "seed_mapping_summary_fast.txt"
    generate_summary_report(results_df, output_summary)

    print("="*70)
    print("PROFILING COMPLETE!")
    print("="*70)

    # Display quick summary
    successful = results_df[results_df['status'].isin(['FASTER', 'SLOWER'])]
    if len(successful) > 0:
        print("\nQUICK SUMMARY:")
        avg_time_per_seed_ms = successful['time_per_seed_ms'].mean()
        avg_time_per_seed_us = successful['time_per_seed_us'].mean()
        avg_speedup = successful['speedup_vs_bwa'].mean()
        faster_count = len(successful[successful['faster_than_bwa'] == True])

        print(f"  Average time/seed: {avg_time_per_seed_ms:.6f} ms ({avg_time_per_seed_us:.2f} us)")
        print(f"  Average speedup vs BWA: {avg_speedup:.2f}x")
        print(f"  Faster than BWA: {faster_count}/{len(successful)} profiles")
        print(f"  Average throughput: {successful['seeds_per_sec'].mean():,.0f} seeds/sec")


if __name__ == "__main__":
    main()
