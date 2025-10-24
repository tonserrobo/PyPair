#!/usr/bin/env python3
"""
Train ALL 8 models on chr21 data for definitive comparison.

This script trains all classifier models on the same chr21 dataset
to prove which models overfit and which generalize.

Usage:
    python scripts/train_all_models_chr21.py \\
        --data data/reference_genomes/chr21_kmers.hdf5 \\
        --output results/chr21_model_comparison.csv

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import argparse
import os
import sys
import time
import json
from pathlib import Path
import joblib
import vaex
import numpy as np
import pandas as pd
from sklearn.metrics import (
    f1_score, precision_score, recall_score, accuracy_score,
    classification_report
)
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier

# Add project root
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


def train_and_evaluate_model(
    model_name: str,
    model,
    X_train, y_train,
    X_test, y_test,
    output_dir: str
):
    """Train and evaluate a single model."""
    print(f"\n{'='*60}")
    print(f"Training: {model_name}")
    print(f"{'='*60}")

    results = {
        'model_name': model_name,
        'train_samples': len(X_train),
        'test_samples': len(X_test),
        'n_features': X_train.shape[1] if hasattr(X_train, 'shape') else len(X_train.columns)
    }

    # Convert Vaex to numpy if needed
    if isinstance(X_train, vaex.dataframe.DataFrame):
        print("  Converting training data to numpy...")
        X_train_np = X_train.values
        y_train_np = y_train.values
    else:
        X_train_np = X_train
        y_train_np = y_train

    if isinstance(X_test, vaex.dataframe.DataFrame):
        print("  Converting test data to numpy...")
        X_test_np = X_test.values
        y_test_np = y_test.values
    else:
        X_test_np = X_test
        y_test_np = y_test

    # Train
    print("  Training model...")
    start_time = time.time()
    try:
        model.fit(X_train_np, y_train_np)
        train_time = time.time() - start_time
        results['train_time_sec'] = train_time
        results['train_status'] = 'success'
        print(f"  Training completed in {train_time:.2f}s")
    except Exception as e:
        print(f"  ERROR during training: {e}")
        results['train_status'] = f'failed: {str(e)}'
        results['train_time_sec'] = time.time() - start_time
        return results

    # Predict on training set (to measure overfitting)
    print("  Evaluating on training set...")
    try:
        start_time = time.time()
        y_train_pred = model.predict(X_train_np)
        train_pred_time = time.time() - start_time

        train_f1 = f1_score(y_train_np, y_train_pred, average='weighted')
        train_accuracy = accuracy_score(y_train_np, y_train_pred)

        results['train_f1'] = train_f1
        results['train_accuracy'] = train_accuracy
        results['train_pred_time_sec'] = train_pred_time
        print(f"    Train F1: {train_f1:.4f}")
        print(f"    Train Accuracy: {train_accuracy:.4f}")
    except Exception as e:
        print(f"    WARNING: Training prediction failed: {e}")
        results['train_f1'] = None
        results['train_accuracy'] = None

    # Predict on test set
    print("  Evaluating on test set...")
    try:
        start_time = time.time()
        y_pred = model.predict(X_test_np)
        pred_time = time.time() - start_time

        f1 = f1_score(y_test_np, y_pred, average='weighted')
        precision = precision_score(y_test_np, y_pred, average='weighted')
        recall = recall_score(y_test_np, y_pred, average='weighted')
        accuracy = accuracy_score(y_test_np, y_pred)

        results['test_f1'] = f1
        results['test_precision'] = precision
        results['test_recall'] = recall
        results['test_accuracy'] = accuracy
        results['pred_time_sec'] = pred_time
        results['pred_time_ms_per_sample'] = pred_time / len(X_test_np) * 1000

        # Calculate overfitting gap
        if results.get('train_f1'):
            results['overfitting_gap'] = results['train_f1'] - f1

        print(f"    Test F1: {f1:.4f}")
        print(f"    Test Precision: {precision:.4f}")
        print(f"    Test Recall: {recall:.4f}")
        print(f"    Test Accuracy: {accuracy:.4f}")
        print(f"    Prediction Time: {pred_time:.2f}s ({pred_time/len(X_test_np)*1000:.4f} ms/sample)")

        if results.get('overfitting_gap'):
            print(f"    Overfitting Gap: {results['overfitting_gap']*100:.2f}%")
            if results['overfitting_gap'] > 0.05:
                print(f"      ⚠️  HIGH OVERFITTING DETECTED!")
    except Exception as e:
        print(f"    ERROR during prediction: {e}")
        results['test_f1'] = None
        results['pred_time_sec'] = time.time() - start_time

    # Save model
    model_path = os.path.join(output_dir, f"{model_name.lower().replace(' ', '_')}_chr21.joblib")
    try:
        joblib.dump(model, model_path)
        model_size_mb = os.path.getsize(model_path) / (1024**2)
        results['model_size_mb'] = model_size_mb
        print(f"  Model saved: {model_path} ({model_size_mb:.2f} MB)")
    except Exception as e:
        print(f"    WARNING: Could not save model: {e}")
        results['model_size_mb'] = None

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Train all 8 models on chr21 for definitive comparison"
    )
    parser.add_argument(
        "--data",
        required=True,
        help="Input HDF5 file with chr21 k-mer data"
    )
    parser.add_argument(
        "--output",
        default="results/chr21_model_comparison.csv",
        help="Output CSV file for results"
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.2,
        help="Test set fraction (default: 0.2)"
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help="Specific models to train (default: all)"
    )

    args = parser.parse_args()

    print("="*60)
    print("PyPair: Multi-Model Chr21 Training")
    print("="*60)
    print(f"Data: {args.data}")
    print(f"Test size: {args.test_size*100:.0f}%")
    print("="*60)

    # Load data
    print("\nLoading data...")
    df = vaex.open(args.data)
    print(f"Total samples: {len(df):,}")

    # Prepare features
    exclude_cols = ['chromosome', 'position', 'kmer', 'strand']
    feature_cols = [col for col in df.columns if col not in exclude_cols]
    print(f"Features: {len(feature_cols)}")

    # Split data
    print(f"\nSplitting data...")
    train_df, test_df = df.split_random([1-args.test_size, args.test_size], random_state=42)

    X_train = train_df[feature_cols]
    y_train = train_df['position']
    X_test = test_df[feature_cols]
    y_test = test_df['position']

    print(f"Training samples: {len(X_train):,}")
    print(f"Test samples: {len(X_test):,}")

    # Define all models
    models = {
        'Gaussian NB': GaussianNB(),
        'Random Forest': RandomForestClassifier(
            n_estimators=100,
            max_depth=20,
            random_state=42,
            n_jobs=-1
        ),
        'Extra Trees': ExtraTreesClassifier(
            n_estimators=100,
            max_depth=20,
            random_state=42,
            n_jobs=-1
        ),
        'Decision Tree': DecisionTreeClassifier(
            max_depth=20,
            random_state=42
        ),
        'KNN': KNeighborsClassifier(
            n_neighbors=5,
            n_jobs=-1
        ),
        'Nearest Centroid': NearestCentroid(),
        'LGBM': LGBMClassifier(
            n_estimators=100,
            max_depth=20,
            random_state=42,
            verbose=-1
        ),
        'XGBoost': XGBClassifier(
            n_estimators=100,
            max_depth=20,
            random_state=42,
            n_jobs=-1,
            verbosity=0
        )
    }

    # Filter models if specified
    if args.models:
        models = {k: v for k, v in models.items() if k in args.models}

    print(f"\nTraining {len(models)} models...")

    # Create output directory
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    os.makedirs("models/model_weights", exist_ok=True)

    # Train all models
    all_results = []
    for model_name, model in models.items():
        try:
            results = train_and_evaluate_model(
                model_name, model,
                X_train, y_train,
                X_test, y_test,
                "models/model_weights"
            )
            all_results.append(results)
        except Exception as e:
            print(f"\nFATAL ERROR with {model_name}: {e}")
            all_results.append({
                'model_name': model_name,
                'train_status': f'fatal_error: {str(e)}'
            })

    # Save results
    print(f"\n{'='*60}")
    print("SAVING RESULTS")
    print(f"{'='*60}")

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(args.output, index=False)
    print(f"Results saved to: {args.output}")

    # Print comparison table
    print(f"\n{'='*60}")
    print("COMPARISON TABLE")
    print(f"{'='*60}")

    # Sort by overfitting gap (low to high)
    if 'overfitting_gap' in results_df.columns:
        results_df_sorted = results_df.sort_values('overfitting_gap')
    else:
        results_df_sorted = results_df

    print("\nModel Performance Summary:")
    print("-" * 80)
    for _, row in results_df_sorted.iterrows():
        print(f"\n{row['model_name']}:")
        if row.get('train_f1'):
            print(f"  Train F1:      {row['train_f1']:.4f}")
        if row.get('test_f1'):
            print(f"  Test F1:       {row['test_f1']:.4f}")
        if row.get('overfitting_gap') is not None:
            print(f"  Overfit Gap:   {row['overfitting_gap']*100:.2f}%", end="")
            if row['overfitting_gap'] > 0.07:
                print(" ❌ HIGH")
            elif row['overfitting_gap'] > 0.03:
                print(" ⚠️  MEDIUM")
            else:
                print(" ✅ LOW")
        if row.get('model_size_mb'):
            print(f"  Model Size:    {row['model_size_mb']:.2f} MB")
        if row.get('train_time_sec'):
            print(f"  Train Time:    {row['train_time_sec']:.2f}s")
        if row.get('pred_time_ms_per_sample'):
            print(f"  Pred Time:     {row['pred_time_ms_per_sample']:.4f} ms/sample")

    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*60}")

    # Identify best model
    if 'overfitting_gap' in results_df.columns:
        best_model = results_df.loc[results_df['overfitting_gap'].idxmin()]
        print(f"\n✅ Best Model (Lowest Overfitting): {best_model['model_name']}")
        print(f"   Overfitting Gap: {best_model['overfitting_gap']*100:.2f}%")
        print(f"   Test F1: {best_model['test_f1']:.4f}")


if __name__ == "__main__":
    main()
