#!/usr/bin/env python3
"""
Train ALL 8 models on S1 baseline dataset for definitive overfitting comparison.

This script trains all classifier models on the S1 dataset with proper PCA-validated
features to prove which models overfit and which generalize well.

Usage:
    python scripts/train_all_models_S1.py

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import os
import sys
import time
import pandas as pd
from pathlib import Path
from sklearn.metrics import (
    f1_score, precision_score, recall_score, accuracy_score
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
    X_test, y_test
):
    """Train and evaluate a single model."""
    print(f"\n{'='*60}")
    print(f"Training: {model_name}")
    print(f"{'='*60}")

    results = {
        'model_name': model_name,
        'train_samples': len(X_train),
        'test_samples': len(X_test),
        'n_features': X_train.shape[1]
    }

    # Train
    print(f"  Training on {len(X_train):,} samples with {X_train.shape[1]} features...")
    start_time = time.time()
    try:
        model.fit(X_train, y_train)
        train_time = time.time() - start_time
        results['train_time_sec'] = train_time
        results['train_status'] = 'success'
        print(f"  [OK] Training completed in {train_time:.2f}s")
    except Exception as e:
        print(f"  [ERROR] Training failed: {e}")
        results['train_status'] = f'failed: {str(e)}'
        results['train_time_sec'] = time.time() - start_time
        return results

    # Predict on training set (to measure overfitting)
    print("  Evaluating on training set...")
    try:
        start_time = time.time()
        y_train_pred = model.predict(X_train)
        train_pred_time = time.time() - start_time

        train_f1 = f1_score(y_train, y_train_pred, average='weighted', zero_division=0)
        train_accuracy = accuracy_score(y_train, y_train_pred)

        results['train_f1'] = train_f1
        results['train_accuracy'] = train_accuracy
        results['train_pred_time_sec'] = train_pred_time
        print(f"    Train F1: {train_f1:.4f} ({train_f1*100:.2f}%)")
        print(f"    Train Accuracy: {train_accuracy:.4f} ({train_accuracy*100:.2f}%)")
    except Exception as e:
        print(f"    WARNING: Training prediction failed: {e}")
        results['train_f1'] = None
        results['train_accuracy'] = None

    # Predict on test set
    print("  Evaluating on test set...")
    try:
        start_time = time.time()
        y_pred = model.predict(X_test)
        pred_time = time.time() - start_time

        test_f1 = f1_score(y_test, y_pred, average='weighted', zero_division=0)
        test_precision = precision_score(y_test, y_pred, average='weighted', zero_division=0)
        test_recall = recall_score(y_test, y_pred, average='weighted', zero_division=0)
        test_accuracy = accuracy_score(y_test, y_pred)

        results['test_f1'] = test_f1
        results['test_precision'] = test_precision
        results['test_recall'] = test_recall
        results['test_accuracy'] = test_accuracy
        results['pred_time_sec'] = pred_time
        results['pred_time_ms_per_sample'] = pred_time / len(X_test) * 1000

        # Calculate overfitting gap
        if results.get('train_f1') is not None:
            results['overfitting_gap'] = results['train_f1'] - test_f1
            results['overfitting_gap_pct'] = results['overfitting_gap'] * 100

        print(f"    Test F1: {test_f1:.4f} ({test_f1*100:.2f}%)")
        print(f"    Test Precision: {test_precision:.4f}")
        print(f"    Test Recall: {test_recall:.4f}")
        print(f"    Test Accuracy: {test_accuracy:.4f}")
        print(f"    Prediction Time: {pred_time:.2f}s ({pred_time/len(X_test)*1000:.4f} ms/sample)")

        if results.get('overfitting_gap'):
            gap_pct = results['overfitting_gap'] * 100
            print(f"\n    Overfitting Gap: {gap_pct:.2f}%")
            if results['overfitting_gap'] > 0.10:
                print(f"      [WARNING] SEVERE OVERFITTING DETECTED (>10%)")
            elif results['overfitting_gap'] > 0.05:
                print(f"      [WARNING] HIGH OVERFITTING DETECTED (>5%)")
            elif results['overfitting_gap'] > 0.01:
                print(f"      [MODERATE] OVERFITTING ({gap_pct:.2f}%)")
            else:
                print(f"      [EXCELLENT] GENERALIZATION ({gap_pct:.2f}%)")
    except Exception as e:
        print(f"    ERROR during prediction: {e}")
        results['test_f1'] = None
        results['pred_time_sec'] = time.time() - start_time

    return results


def main():
    print("="*60)
    print("PyPair: Multi-Model S1 Baseline Training")
    print("="*60)
    print("Dataset: S1 (PCA-validated features)")
    print("="*60)

    # Load S1 training and test data (cleaned)
    data_path = Path("data/reference_samples/cleaned")

    print("\nLoading S1 cleaned training data...")
    X_train = pd.read_csv(data_path / "S1_extended_encoded_x_train.csv")
    y_train = pd.read_csv(data_path / "S1_extended_encoded_y_train.csv")['org_location']

    print("Loading S1 cleaned test data...")
    X_test = pd.read_csv(data_path / "S1_extended_encoded_X_test.csv")
    y_test = pd.read_csv(data_path / "S1_extended_encoded_y_test.csv")['org_location']

    print(f"\nDataset loaded:")
    print(f"  Training samples: {len(X_train):,}")
    print(f"  Test samples: {len(X_test):,}")
    print(f"  Features: {X_train.shape[1]}")
    print(f"  Feature names: {list(X_train.columns)}")

    # Define all 8 models
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
            verbosity=0
        )
    }

    print(f"\nTraining {len(models)} models...")

    # Train all models
    all_results = []
    for model_name, model in models.items():
        try:
            results = train_and_evaluate_model(
                model_name, model,
                X_train, y_train,
                X_test, y_test
            )
            all_results.append(results)
        except Exception as e:
            print(f"\n[FATAL ERROR] with {model_name}: {e}")
            all_results.append({
                'model_name': model_name,
                'train_status': f'fatal: {str(e)}'
            })

    # Save results
    print(f"\n{'='*60}")
    print("SAVING RESULTS")
    print(f"{'='*60}")

    results_df = pd.DataFrame(all_results)
    output_file = "results/S1_model_comparison.csv"
    os.makedirs("results", exist_ok=True)
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")

    # Print comparison table
    print(f"\n{'='*60}")
    print("OVERFITTING COMPARISON TABLE")
    print(f"{'='*60}")
    print(f"\n{'Model':<20} {'Train F1':<10} {'Test F1':<10} {'Gap %':<10} {'Status'}")
    print("-" * 70)

    for _, row in results_df.iterrows():
        model_name = row['model_name']
        train_f1 = row.get('train_f1', None)
        test_f1 = row.get('test_f1', None)
        gap = row.get('overfitting_gap_pct', None)

        if train_f1 is not None and test_f1 is not None and gap is not None:
            if gap > 10:
                status = "[SEVERE]"
            elif gap > 5:
                status = "[HIGH]"
            elif gap > 1:
                status = "[MODERATE]"
            else:
                status = "[EXCELLENT]"
            print(f"{model_name:<20} {train_f1:.4f}     {test_f1:.4f}     {gap:>6.2f}%    {status}")
        else:
            print(f"{model_name:<20} FAILED")

    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*60}")
    print(f"\nKey Findings:")
    print(f"  - Models with gap < 1%: EXCELLENT generalization (suitable)")
    print(f"  - Models with gap > 5%: HIGH overfitting (likely unsuitable)")
    print(f"  - Models with gap > 10%: SEVERE overfitting (definitely unsuitable)")


if __name__ == "__main__":
    main()
