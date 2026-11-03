#!/usr/bin/env python3
"""
Train Gaussian Naive Bayes on chromosome k-mer data.

This script trains a GNB model using Vaex for out-of-core processing,
enabling training on datasets larger than available RAM.

Usage:
    python train_on_chromosomes.py \\
        --data data/reference_genomes/chr21_kmers.hdf5 \\
        --output models/model_weights/gnb_chr21.joblib \\
        --test-size 0.2

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import argparse
import os
import sys
import time
from pathlib import Path
import joblib
import vaex
import numpy as np
from sklearn.metrics import (
    f1_score, precision_score, recall_score, accuracy_score,
    classification_report, confusion_matrix
)

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from accelerated_GNB.GNB_vaex.GNB_HDF5_optimised import VaexGaussianNB
from scripts.vaex_utils import to_numpy


def train_gnb_on_chromosomes(
    data_file: str,
    output_model: str,
    test_size: float = 0.2,
    random_state: int = 42
):
    """
    Train Gaussian Naive Bayes on chromosome data.

    Args:
        data_file: Path to HDF5 file with k-mer data
        output_model: Path to save trained model
        test_size: Fraction of data for testing
        random_state: Random seed for reproducibility
    """
    print("="*60)
    print("PyPair: Chromosome-Based GNB Training")
    print("="*60)
    print(f"Data file: {data_file}")
    print(f"Output model: {output_model}")
    print(f"Test size: {test_size*100:.0f}%")
    print(f"Random seed: {random_state}")
    print("="*60)

    # Load data
    print("\n[1/5] Loading data...")
    df = vaex.open(data_file)
    print(f"  Total samples: {len(df):,}")
    print(f"  Columns: {list(df.columns)}")

    # Identify feature columns (exclude metadata)
    exclude_cols = ['chromosome', 'position', 'kmer', 'strand']
    feature_cols = [col for col in df.columns if col not in exclude_cols]
    print(f"  Feature columns ({len(feature_cols)}): {feature_cols}")

    # Split data
    print(f"\n[2/5] Splitting data (train: {(1-test_size)*100:.0f}%, test: {test_size*100:.0f}%)...")
    train_df, test_df = df.split_random([1-test_size, test_size], random_state=random_state)

    X_train = train_df[feature_cols]
    y_train = train_df['position']  # Use genomic position as target
    X_test = test_df[feature_cols]
    y_test = test_df['position']

    print(f"  Training samples: {len(X_train):,}")
    print(f"  Test samples: {len(X_test):,}")

    # Train model
    print("\n[3/5] Training Gaussian Naive Bayes (Vaex)...")
    model = VaexGaussianNB()

    start_time = time.time()
    model.fit(X_train, y_train)
    train_time = time.time() - start_time

    print(f"  Training completed in {train_time:.2f} seconds ({train_time/60:.1f} minutes)")

    # Evaluate on test set
    print("\n[4/5] Evaluating model on test set...")
    start_time = time.time()
    y_pred = model.predict_chunked(X_test, chunk_size=10000)
    pred_time = time.time() - start_time

    print(f"  Prediction completed in {pred_time:.2f} seconds")
    print(f"  Prediction speed: {pred_time/len(X_test)*1000:.4f} ms per sample")

    # Compute metrics
    print("\n  Computing metrics...")

    # Convert to numpy if needed
    y_test_np = to_numpy(y_test)

    # For genomic position prediction, we need to discretize or use regression metrics
    # For now, let's compute classification metrics if positions are discrete
    # Otherwise, compute MAE/RMSE

    unique_positions = len(np.unique(y_test_np))
    print(f"  Unique positions in test set: {unique_positions:,}")

    if unique_positions < 1000:  # Treat as classification
        f1 = f1_score(y_test_np, y_pred, average='weighted')
        precision = precision_score(y_test_np, y_pred, average='weighted')
        recall = recall_score(y_test_np, y_pred, average='weighted')
        accuracy = accuracy_score(y_test_np, y_pred)

        print("\n" + "="*60)
        print("CLASSIFICATION METRICS")
        print("="*60)
        print(f"F1 Score:        {f1:.4f}")
        print(f"Precision:       {precision:.4f}")
        print(f"Recall:          {recall:.4f}")
        print(f"Accuracy:        {accuracy:.4f}")
        print("="*60)

        metrics = {
            'f1_score': f1,
            'precision': precision,
            'recall': recall,
            'accuracy': accuracy,
            'train_time_sec': train_time,
            'pred_time_sec': pred_time,
            'pred_time_ms_per_sample': pred_time / len(X_test) * 1000,
            'train_samples': len(X_train),
            'test_samples': len(X_test),
            'n_features': len(feature_cols)
        }

    else:  # Treat as regression
        from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

        mae = mean_absolute_error(y_test_np, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test_np, y_pred))
        r2 = r2_score(y_test_np, y_pred)

        print("\n" + "="*60)
        print("REGRESSION METRICS")
        print("="*60)
        print(f"MAE:             {mae:.2f} bp")
        print(f"RMSE:            {rmse:.2f} bp")
        print(f"R² Score:        {r2:.4f}")
        print("="*60)

        metrics = {
            'mae': mae,
            'rmse': rmse,
            'r2_score': r2,
            'train_time_sec': train_time,
            'pred_time_sec': pred_time,
            'pred_time_ms_per_sample': pred_time / len(X_test) * 1000,
            'train_samples': len(X_train),
            'test_samples': len(X_test),
            'n_features': len(feature_cols)
        }

    # Save model
    print(f"\n[5/5] Saving model to {output_model}...")
    os.makedirs(os.path.dirname(output_model), exist_ok=True)
    joblib.dump(model, output_model)

    model_size_mb = os.path.getsize(output_model) / (1024**2)
    print(f"  Model size: {model_size_mb:.2f} MB")

    # Save metrics
    metrics_file = output_model.replace('.joblib', '_metrics.json')
    import json
    with open(metrics_file, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"  Metrics saved to: {metrics_file}")

    print("\n" + "="*60)
    print("Training Complete!")
    print("="*60)

    return model, metrics


def main():
    parser = argparse.ArgumentParser(
        description="Train GNB on chromosome k-mer data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Train on chr21
    python train_on_chromosomes.py \\
        --data data/reference_genomes/chr21_kmers.hdf5 \\
        --output models/model_weights/gnb_chr21.joblib

    # Train on chr21+22 with custom test size
    python train_on_chromosomes.py \\
        --data data/reference_genomes/chr21_22_kmers.hdf5 \\
        --output models/model_weights/gnb_chr21_22.joblib \\
        --test-size 0.15
        """
    )

    parser.add_argument(
        "--data",
        required=True,
        help="Input HDF5 file with k-mer data"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output model file path (.joblib)"
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.2,
        help="Fraction of data for testing (default: 0.2)"
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)"
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.data):
        print(f"ERROR: Data file not found: {args.data}")
        sys.exit(1)

    if not args.output.endswith('.joblib'):
        print("ERROR: Output file must have .joblib extension")
        sys.exit(1)

    # Train model
    model, metrics = train_gnb_on_chromosomes(
        data_file=args.data,
        output_model=args.output,
        test_size=args.test_size,
        random_state=args.random_state
    )


if __name__ == "__main__":
    main()
