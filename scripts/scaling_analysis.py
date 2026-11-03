"""
PyPair: Scaling Analysis - KNN vs GNB Memory Comparison
========================================================

Tests the hypothesis that GNB will outperform KNN at genome scale due to
lower memory requirements.

Memory Complexity Analysis:
---------------------------
KNN:
    - Training: O(n × d) where n=samples, d=features
    - Must store ALL training samples
    - Prediction: O(n × d × k) for k-NN search
    - Memory: ~1GB per 1M samples (11 features × 8 bytes)

GNB:
    - Training: O(n × d)
    - Stores only class statistics (mean, variance per feature per class)
    - Prediction: O(c × d) where c=classes
    - Memory: O(c × d) = 10,000 classes × 11 features × 8 bytes = 0.88 MB

Expected Outcome:
    - Small data (S1-S4: ~45k samples): KNN competitive
    - Medium data (chr21: ~80M k-mers): KNN struggles, GNB maintains
    - Large data (genome: ~3B k-mers): GNB only viable option

Usage:
    python scripts/scaling_analysis.py

Output:
    results/scaling_analysis.csv - Memory and time metrics at each scale
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time
import psutil
import os
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
import joblib


# Configuration
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Data scale progression
SCALES = {
    'S1': {
        'train_path': 'data/reference_samples/cleaned/S1_extended_encoded_X_train.csv',
        'test_path': 'data/reference_samples/cleaned/S1_extended_encoded_X_test.csv',
        'y_train_path': 'data/reference_samples/cleaned/S1_extended_encoded_y_train.csv',
        'y_test_path': 'data/reference_samples/cleaned/S1_extended_encoded_y_test.csv',
        'expected_samples': 44_800,
        'description': 'Baseline (44.8k samples)'
    },
    'S1+S2': {
        'description': 'Combined S1+S2 (~90k samples)'
    },
    'S1+S2+S3': {
        'description': 'Combined S1+S2+S3 (~135k samples)'
    },
    'S1+S2+S3+S4': {
        'description': 'All baseline datasets (~180k samples)'
    }
}


def get_memory_usage_mb():
    """Get current process memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 ** 2)


def load_dataset(scale_name):
    """Load dataset for given scale."""
    print(f"Loading dataset: {scale_name}")

    if scale_name == 'S1':
        config = SCALES['S1']
        X_train = pd.read_csv(config['train_path']).values
        y_train = pd.read_csv(config['y_train_path'])['org_location'].values
        X_test = pd.read_csv(config['test_path']).values
        y_test = pd.read_csv(config['y_test_path'])['org_location'].values

    elif scale_name == 'S1+S2':
        # Load and concatenate S1 and S2
        X_train_s1 = pd.read_csv('data/reference_samples/cleaned/S1_extended_encoded_X_train.csv').values
        y_train_s1 = pd.read_csv('data/reference_samples/cleaned/S1_extended_encoded_y_train.csv')['org_location'].values
        X_test_s1 = pd.read_csv('data/reference_samples/cleaned/S1_extended_encoded_X_test.csv').values
        y_test_s1 = pd.read_csv('data/reference_samples/cleaned/S1_extended_encoded_y_test.csv')['org_location'].values

        X_train_s2 = pd.read_csv('data/reference_samples/cleaned/S2_extended_encoded_X_train.csv').values
        y_train_s2 = pd.read_csv('data/reference_samples/cleaned/S2_extended_encoded_y_train.csv')['org_location'].values
        X_test_s2 = pd.read_csv('data/reference_samples/cleaned/S2_extended_encoded_X_test.csv').values
        y_test_s2 = pd.read_csv('data/reference_samples/cleaned/S2_extended_encoded_y_test.csv')['org_location'].values

        X_train = np.vstack([X_train_s1, X_train_s2])
        y_train = np.concatenate([y_train_s1, y_train_s2])
        X_test = np.vstack([X_test_s1, X_test_s2])
        y_test = np.concatenate([y_test_s1, y_test_s2])

    elif scale_name == 'S1+S2+S3':
        # Load S1+S2+S3
        datasets = ['S1', 'S2', 'S3']
        X_train_list = []
        y_train_list = []
        X_test_list = []
        y_test_list = []

        for ds in datasets:
            X_train_list.append(pd.read_csv(f'data/reference_samples/cleaned/{ds}_extended_encoded_X_train.csv').values)
            y_train_list.append(pd.read_csv(f'data/reference_samples/cleaned/{ds}_extended_encoded_y_train.csv')['org_location'].values)
            X_test_list.append(pd.read_csv(f'data/reference_samples/cleaned/{ds}_extended_encoded_X_test.csv').values)
            y_test_list.append(pd.read_csv(f'data/reference_samples/cleaned/{ds}_extended_encoded_y_test.csv')['org_location'].values)

        X_train = np.vstack(X_train_list)
        y_train = np.concatenate(y_train_list)
        X_test = np.vstack(X_test_list)
        y_test = np.concatenate(y_test_list)

    elif scale_name == 'S1+S2+S3+S4':
        # Load all 4 datasets
        datasets = ['S1', 'S2', 'S3', 'S4']
        X_train_list = []
        y_train_list = []
        X_test_list = []
        y_test_list = []

        for ds in datasets:
            X_train_list.append(pd.read_csv(f'data/reference_samples/cleaned/{ds}_extended_encoded_X_train.csv').values)
            y_train_list.append(pd.read_csv(f'data/reference_samples/cleaned/{ds}_extended_encoded_y_train.csv')['org_location'].values)
            X_test_list.append(pd.read_csv(f'data/reference_samples/cleaned/{ds}_extended_encoded_X_test.csv').values)
            y_test_list.append(pd.read_csv(f'data/reference_samples/cleaned/{ds}_extended_encoded_y_test.csv')['org_location'].values)

        X_train = np.vstack(X_train_list)
        y_train = np.concatenate(y_train_list)
        X_test = np.vstack(X_test_list)
        y_test = np.concatenate(y_test_list)

    print(f"[OK] Loaded {len(X_train):,} training samples, {len(X_test):,} test samples")
    return X_train, y_train, X_test, y_test


def profile_model(model_name, model, X_train, y_train, X_test, y_test):
    """
    Profile model training and prediction with memory tracking.

    Returns dict with:
        - train_time_sec
        - train_memory_mb
        - predict_time_sec
        - predict_memory_mb
        - model_size_mb (serialized)
    """
    print(f"\nProfiling: {model_name}")
    print(f"  Training samples: {len(X_train):,}")

    results = {
        'model_name': model_name,
        'train_samples': len(X_train),
        'test_samples': len(X_test),
        'n_features': X_train.shape[1]
    }

    # Measure training
    print("  Training...")
    mem_before_train = get_memory_usage_mb()

    start_time = time.time()
    model.fit(X_train, y_train)
    train_time = time.time() - start_time

    mem_after_train = get_memory_usage_mb()
    train_memory_delta = mem_after_train - mem_before_train

    results['train_time_sec'] = train_time
    results['train_memory_mb'] = train_memory_delta
    results['train_memory_per_sample_kb'] = (train_memory_delta * 1024) / len(X_train)

    print(f"  [OK] Training: {train_time:.2f}s, Memory: {train_memory_delta:.2f} MB")

    # Measure prediction
    print("  Predicting...")
    mem_before_predict = get_memory_usage_mb()

    start_time = time.time()
    _ = model.predict(X_test)
    predict_time = time.time() - start_time

    mem_after_predict = get_memory_usage_mb()
    predict_memory_delta = mem_after_predict - mem_before_predict

    results['predict_time_sec'] = predict_time
    results['predict_memory_mb'] = predict_memory_delta
    results['predict_time_per_sample_ms'] = (predict_time / len(X_test)) * 1000

    print(f"  [OK] Prediction: {predict_time:.2f}s, Memory: {predict_memory_delta:.2f} MB")

    # Measure model size (serialized)
    print("  Measuring model size...")
    temp_path = Path("temp_model.joblib")
    joblib.dump(model, temp_path)
    model_size_mb = temp_path.stat().st_size / (1024 ** 2)
    temp_path.unlink()

    results['model_size_mb'] = model_size_mb
    results['model_size_per_sample_kb'] = (model_size_mb * 1024) / len(X_train)

    print(f"  [OK] Model size: {model_size_mb:.2f} MB")

    return results


def run_scaling_experiment():
    """Run scaling experiment across all data scales."""
    print("="*80)
    print("PyPair: Scaling Analysis - KNN vs GNB")
    print("="*80)
    print("Hypothesis: GNB scales better than KNN due to lower memory overhead")
    print("="*80)

    all_results = []

    # Test each scale
    for scale_name, config in SCALES.items():
        print(f"\n{'#'*80}")
        print(f"# Scale: {scale_name} - {config['description']}")
        print(f"{'#'*80}")

        # Load data
        try:
            X_train, y_train, X_test, y_test = load_dataset(scale_name)
        except FileNotFoundError as e:
            print(f"[SKIP] Dataset not available: {e}")
            continue

        # Profile KNN
        knn = KNeighborsClassifier(n_neighbors=5, weights='distance', n_jobs=1)
        knn_results = profile_model('KNN', knn, X_train, y_train, X_test, y_test)
        knn_results['scale'] = scale_name
        knn_results['scale_description'] = config['description']
        all_results.append(knn_results)

        # Profile GNB
        gnb = GaussianNB()
        gnb_results = profile_model('Gaussian NB', gnb, X_train, y_train, X_test, y_test)
        gnb_results['scale'] = scale_name
        gnb_results['scale_description'] = config['description']
        all_results.append(gnb_results)

        # Calculate efficiency metrics
        print(f"\n  COMPARISON:")
        print(f"  {'Metric':<30} {'KNN':<20} {'GNB':<20} {'Winner'}")
        print(f"  {'-'*75}")

        knn_train_mem = knn_results['train_memory_mb']
        gnb_train_mem = gnb_results['train_memory_mb']
        mem_winner = 'GNB' if gnb_train_mem < knn_train_mem else 'KNN'
        print(f"  {'Training Memory (MB)':<30} {knn_train_mem:<20.2f} {gnb_train_mem:<20.2f} {mem_winner}")

        knn_model_size = knn_results['model_size_mb']
        gnb_model_size = gnb_results['model_size_mb']
        size_winner = 'GNB' if gnb_model_size < knn_model_size else 'KNN'
        print(f"  {'Model Size (MB)':<30} {knn_model_size:<20.2f} {gnb_model_size:<20.2f} {size_winner}")

        knn_train_time = knn_results['train_time_sec']
        gnb_train_time = gnb_results['train_time_sec']
        time_winner = 'GNB' if gnb_train_time < knn_train_time else 'KNN'
        print(f"  {'Training Time (sec)':<30} {knn_train_time:<20.2f} {gnb_train_time:<20.2f} {time_winner}")

    # Save results
    results_df = pd.DataFrame(all_results)
    output_file = RESULTS_DIR / "scaling_analysis.csv"
    results_df.to_csv(output_file, index=False)

    print(f"\n{'='*80}")
    print("SCALING ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"Results saved: {output_file}")
    print(f"Total experiments: {len(all_results)}")

    # Summary analysis
    print(f"\n{'='*80}")
    print("HYPOTHESIS TEST RESULTS")
    print(f"{'='*80}")

    for scale in results_df['scale'].unique():
        scale_data = results_df[results_df['scale'] == scale]
        knn_data = scale_data[scale_data['model_name'] == 'KNN'].iloc[0]
        gnb_data = scale_data[scale_data['model_name'] == 'Gaussian NB'].iloc[0]

        print(f"\n{scale} ({knn_data['train_samples']:,} samples):")
        print(f"  Memory: KNN={knn_data['model_size_mb']:.2f}MB, GNB={gnb_data['model_size_mb']:.2f}MB")
        print(f"  Ratio: KNN is {knn_data['model_size_mb']/gnb_data['model_size_mb']:.1f}x larger than GNB")

    print(f"\n{'='*80}")
    print("CONCLUSION:")
    print(f"{'='*80}")

    # Calculate average scaling factor
    knn_avg_size = results_df[results_df['model_name'] == 'KNN']['model_size_mb'].mean()
    gnb_avg_size = results_df[results_df['model_name'] == 'Gaussian NB']['model_size_mb'].mean()
    scaling_factor = knn_avg_size / gnb_avg_size

    print(f"Average model size:")
    print(f"  KNN: {knn_avg_size:.2f} MB")
    print(f"  GNB: {gnb_avg_size:.2f} MB")
    print(f"  Scaling factor: {scaling_factor:.1f}x")
    print(f"\nHypothesis {'CONFIRMED' if scaling_factor > 10 else 'NEEDS MORE DATA'}:")
    print(f"  GNB scales better than KNN for genome-scale data")
    print(f"{'='*80}")


if __name__ == "__main__":
    run_scaling_experiment()
