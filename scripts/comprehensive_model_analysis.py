#!/usr/bin/env python3
"""
Comprehensive model analysis for ML seed mapping.

This script:
1. Trains all 8 models on S1-S4 datasets with memory-safe configurations
2. Profiles inference performance at scale (50k-500k seed batches)
3. Compares against BWA benchmark (0.003 sec/seed)
4. Handles tree-based models with pruning to prevent memory overflow
5. Integrates memory profiling using existing benchmark infrastructure
6. Runs SHAP explainability analysis on leading models

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import os
import sys
import time
import joblib
import pandas as pd
import numpy as np
from pathlib import Path
from memory_profiler import memory_usage
from sklearn.metrics import f1_score, accuracy_score
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier

# Add project root
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import SHAP analysis, model fixes, and non-inferiority testing
from scripts.utils.shap_analysis import analyze_model_with_shap
from scripts.model_fixes_and_limitations import (
    get_fixed_models,
    create_label_encoder_wrapper,
    LITERATURE_NOTES
)
from scripts.noninferiority_testing import (
    comprehensive_noninferiority_analysis,
    generate_noninferiority_report
)

# BWA benchmark for comparison
BWA_SEC_PER_SEED = 0.003

# Batch sizes for incremental testing
BATCH_SIZES = [50000, 100000, 200000, 500000]


def get_memory_safe_models_old():
    """Get model configurations with memory constraints for tree-based models."""
    return {
        'Gaussian NB': GaussianNB(),

        'Random Forest': RandomForestClassifier(
            n_estimators=50,  # Reduced from 100
            max_depth=15,     # Reduced from 20
            max_leaf_nodes=500,  # Limit tree size
            max_samples=0.5,  # Use only 50% of data per tree
            random_state=42,
            n_jobs=-1
        ),

        'Extra Trees': ExtraTreesClassifier(
            n_estimators=50,
            max_depth=15,
            max_leaf_nodes=500,
            max_samples=0.5,
            random_state=42,
            n_jobs=-1
        ),

        'Decision Tree': DecisionTreeClassifier(
            max_depth=15,
            max_leaf_nodes=500,  # Limit leaves to prevent memory explosion
            random_state=42
        ),

        'KNN': KNeighborsClassifier(
            n_neighbors=5,
            n_jobs=-1
        ),

        'Nearest Centroid': NearestCentroid(),

        'LGBM': LGBMClassifier(
            n_estimators=50,  # Reduced from 100
            max_depth=15,
            num_leaves=200,   # Limit leaves
            random_state=42,
            verbose=-1
        ),

        'XGBoost': XGBClassifier(
            n_estimators=50,
            max_depth=15,
            max_leaves=200,
            random_state=42,
            verbosity=0
        )
    }


def profile_batch_inference(model, X_batch, batch_sizes=None):
    """
    Profile inference performance at different batch sizes with memory tracking.

    Args:
        model: Trained classifier model
        X_batch: Feature matrix for testing
        batch_sizes: List of batch sizes to test (default: BATCH_SIZES)

    Returns:
        dict with timing and memory stats for each batch size
    """
    if batch_sizes is None:
        batch_sizes = BATCH_SIZES

    results = {}

    for batch_size in batch_sizes:
        if batch_size > len(X_batch):
            print(f"      Batch {batch_size:>6,}: SKIPPED (exceeds dataset size)")
            continue

        X_subset = X_batch[:batch_size]

        # Warm-up run
        try:
            _ = model.predict(X_subset[:min(10, len(X_subset))])
        except:
            pass

        # Measure memory usage
        mem_usage_mb = 0
        try:
            def predict_wrapper():
                return model.predict(X_subset)

            mem_usage = memory_usage((predict_wrapper,), interval=0.1, max_iterations=1)
            mem_usage_mb = max(mem_usage) if mem_usage else 0
        except Exception as e:
            print(f"      Memory profiling failed for batch {batch_size}: {e}")

        # Timed runs (average of 3)
        times = []
        for _ in range(3):
            start = time.time()
            try:
                _ = model.predict(X_subset)
                elapsed = time.time() - start
                times.append(elapsed)
            except Exception as e:
                print(f"      Batch size {batch_size}: ERROR - {e}")
                break

        if times:
            avg_time = np.mean(times)
            std_time = np.std(times)
            time_per_seed = avg_time / batch_size

            # Compare to BWA
            speedup = BWA_SEC_PER_SEED / time_per_seed if time_per_seed > 0 else 0

            results[batch_size] = {
                'batch_size': batch_size,
                'total_time_sec': avg_time,
                'time_std_sec': std_time,
                'time_per_seed_sec': time_per_seed,
                'time_per_seed_ms': time_per_seed * 1000,
                'seeds_per_sec': batch_size / avg_time if avg_time > 0 else 0,
                'memory_mb': mem_usage_mb,
                'memory_per_seed_kb': (mem_usage_mb * 1024) / batch_size if batch_size > 0 else 0,
                'vs_bwa_speedup': speedup,
                'vs_bwa_status': 'FASTER' if speedup > 1 else 'SLOWER'
            }

            print(f"      Batch {batch_size:>6,}: {time_per_seed*1000:.6f} ms/seed "
                  f"({speedup:>6.2f}x vs BWA) | Mem: {mem_usage_mb:.2f} MB")

    return results


def train_and_evaluate_model(
    model_name: str,
    model,
    X_train, y_train,
    X_test, y_test,
    dataset_name: str,
    save_model: bool = True,
    limitations: dict = None
):
    """Train, evaluate, and profile a single model."""
    print(f"\n{'='*70}")
    print(f"Training: {model_name} on {dataset_name}")
    print(f"{'='*70}")

    results = {
        'dataset': dataset_name,
        'model_name': model_name,
        'train_samples': len(X_train),
        'test_samples': len(X_test),
        'n_features': X_train.shape[1]
    }

    # Add limitations to results for documentation
    if limitations:
        results['limitation'] = limitations['limitation']
        results['overfitting_risk'] = limitations['overfitting_risk']
        results['justification'] = limitations['justification']

    # Train
    print(f"  Training on {len(X_train):,} samples...")
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
        return results, None  # Return tuple with None for failed model

    # Evaluate on training set
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
    except Exception as e:
        print(f"    WARNING: Training prediction failed: {e}")
        results['train_f1'] = None
        results['train_accuracy'] = None

    # Evaluate on test set
    print("  Evaluating on test set...")
    try:
        start_time = time.time()
        y_pred = model.predict(X_test)
        pred_time = time.time() - start_time

        test_f1 = f1_score(y_test, y_pred, average='weighted', zero_division=0)
        test_accuracy = accuracy_score(y_test, y_pred)

        results['test_f1'] = test_f1
        results['test_accuracy'] = test_accuracy
        results['pred_time_sec'] = pred_time
        results['pred_time_ms_per_sample'] = pred_time / len(X_test) * 1000

        # Calculate overfitting gap
        if results.get('train_f1') is not None:
            results['overfitting_gap'] = results['train_f1'] - test_f1
            results['overfitting_gap_pct'] = results['overfitting_gap'] * 100

        print(f"    Test F1: {test_f1:.4f} ({test_f1*100:.2f}%)")
        print(f"    Test Accuracy: {test_accuracy:.4f}")

        if results.get('overfitting_gap'):
            gap_pct = results['overfitting_gap'] * 100
            print(f"    Overfitting Gap: {gap_pct:.2f}%")
    except Exception as e:
        print(f"    ERROR during prediction: {e}")
        results['test_f1'] = None

    # Profile batch inference
    print(f"  Profiling batch inference across {len(BATCH_SIZES)} batch sizes...")
    print(f"  BWA benchmark: {BWA_SEC_PER_SEED*1000:.3f} ms/seed")
    try:
        batch_profile = profile_batch_inference(model, X_test)

        # Store results for each batch size
        for batch_size, metrics in batch_profile.items():
            results[f'batch_{batch_size//1000}k_time_ms'] = metrics['time_per_seed_ms']
            results[f'batch_{batch_size//1000}k_speedup'] = metrics['vs_bwa_speedup']
            results[f'batch_{batch_size//1000}k_mem_mb'] = metrics['memory_mb']

        # Store all batch profile data
        results['batch_profile'] = batch_profile

    except Exception as e:
        print(f"    WARNING: Batch profiling failed: {e}")
        results['batch_profile'] = None

    # Save trained model
    if save_model and results.get('train_status') == 'success':
        model_dir = Path("models/model_weights")
        model_dir.mkdir(parents=True, exist_ok=True)
        model_file = model_dir / f"{model_name.lower().replace(' ', '_')}_{dataset_name.lower()}.joblib"

        try:
            joblib.dump(model, model_file)
            results['model_path'] = str(model_file)
            print(f"  [OK] Model saved: {model_file}")
        except Exception as e:
            print(f"    WARNING: Failed to save model: {e}")
            results['model_path'] = None

    return results, model


def main():
    print("="*70)
    print("PyPair: Comprehensive Multi-Model Analysis (S1-S4)")
    print("="*70)
    print("Features:")
    print("  - Memory-safe tree configurations (pruned branches)")
    print(f"  - Batch inference profiling ({', '.join([f'{b//1000}k' for b in BATCH_SIZES])} seeds)")
    print("  - Memory and execution time tracking")
    print("  - BWA comparison (0.003 sec/seed benchmark)")
    print("  - SHAP explainability for leading models")
    print("="*70)

    data_path = Path("data/reference_samples/cleaned")

    # Process each dataset
    datasets = ['S1', 'S2', 'S3', 'S4']
    all_results = []
    trained_models = {}  # Global dict for SHAP analysis later

    for dataset_name in datasets:
        print(f"\n{'#'*70}")
        print(f"# Processing Dataset: {dataset_name}")
        print(f"{'#'*70}")

        # Load data
        try:
            print(f"\nLoading {dataset_name} data...")
            X_train = pd.read_csv(data_path / f"{dataset_name}_extended_encoded_x_train.csv")
            y_train = pd.read_csv(data_path / f"{dataset_name}_extended_encoded_y_train.csv")['org_location']
            X_test = pd.read_csv(data_path / f"{dataset_name}_extended_encoded_X_test.csv")
            y_test = pd.read_csv(data_path / f"{dataset_name}_extended_encoded_y_test.csv")['org_location']

            print(f"  Training: {len(X_train):,} samples")
            print(f"  Test: {len(X_test):,} samples")
            print(f"  Features: {X_train.shape[1]}")
        except Exception as e:
            print(f"  ERROR loading {dataset_name}: {e}")
            continue

        # Train all models with fixes and limitations
        fixed_models = get_fixed_models()

        for model_name, (model, limitations) in fixed_models.items():
            # Log model limitations
            print(f"\n[LIMITATION] {model_name}:")
            print(f"  {limitations['limitation']}")
            print(f"  Expected: {limitations['performance_expectation']}")
            print(f"  Justification: {limitations['justification'][:80]}...")

            # Apply label encoding wrapper for XGBoost/LGBM
            if 'XGBoost' in model_name or 'LGBM' in model_name:
                model = create_label_encoder_wrapper(model, model_name)
                print(f"  [FIX] Applied label encoding wrapper (1-indexed -> 0-indexed)")
            try:
                results, trained_model = train_and_evaluate_model(
                    model_name, model,
                    X_train, y_train,
                    X_test, y_test,
                    dataset_name,
                    save_model=True,
                    limitations=limitations
                )
                all_results.append(results)

                # Store successful models for SHAP
                if results.get('train_status') == 'success' and results.get('test_f1'):
                    key = f"{dataset_name}_{model_name}"
                    trained_models[key] = {
                        'model': trained_model,
                        'model_path': results.get('model_path'),
                        'test_f1': results.get('test_f1'),
                        'overfitting_gap': results.get('overfitting_gap', 999),
                        'X_test': X_test,
                        'model_name': model_name,
                        'dataset': dataset_name
                    }

            except Exception as e:
                print(f"\n[FATAL ERROR] {model_name} on {dataset_name}: {e}")
                all_results.append({
                    'dataset': dataset_name,
                    'model_name': model_name,
                    'train_status': f'fatal: {str(e)}'
                })

    # Save comprehensive results
    print(f"\n{'='*70}")
    print("SAVING RESULTS")
    print(f"{'='*70}")

    # Save main results
    results_df = pd.DataFrame(all_results)
    output_file = "results/comprehensive_model_analysis.csv"
    os.makedirs("results", exist_ok=True)
    results_df.to_csv(output_file, index=False)
    print(f"Main results saved to: {output_file}")

    # Print summary tables
    print(f"\n{'='*70}")
    print("OVERFITTING SUMMARY (All Datasets)")
    print(f"{'='*70}")

    summary_df = results_df[results_df['train_f1'].notna()].copy()
    if not summary_df.empty:
        summary_df = summary_df[['dataset', 'model_name', 'train_f1', 'test_f1', 'overfitting_gap_pct']]
        summary_df = summary_df.sort_values(['dataset', 'overfitting_gap_pct'])

        for dataset in datasets:
            dataset_results = summary_df[summary_df['dataset'] == dataset]
            if not dataset_results.empty:
                print(f"\n{dataset}:")
                print(f"{'Model':<20} {'Train F1':<10} {'Test F1':<10} {'Gap %':<10}")
                print("-" * 50)
                for _, row in dataset_results.iterrows():
                    print(f"{row['model_name']:<20} {row['train_f1']:.4f}     "
                          f"{row['test_f1']:.4f}     {row['overfitting_gap_pct']:>6.2f}%")

    # Print inference performance summary
    print(f"\n{'='*70}")
    print("INFERENCE PERFORMANCE vs BWA")
    print(f"{'='*70}")
    print(f"BWA Benchmark: {BWA_SEC_PER_SEED*1000:.3f} ms/seed")
    print()

    # Check if batch profiling columns exist
    if 'batch_100k_time_per_seed_ms' in results_df.columns:
        perf_df = results_df[results_df['batch_100k_time_per_seed_ms'].notna()].copy()
        if not perf_df.empty:
            perf_df = perf_df[['dataset', 'model_name', 'batch_100k_time_per_seed_ms', 'batch_100k_vs_bwa_speedup']]
            perf_df = perf_df.sort_values('batch_100k_time_per_seed_ms')

            print(f"{'Dataset':<8} {'Model':<20} {'ms/seed':<12} {'vs BWA':<10}")
            print("-" * 50)
            for _, row in perf_df.iterrows():
                speedup = row['batch_100k_vs_bwa_speedup']
                status = "FASTER" if speedup > 1 else "SLOWER"
                print(f"{row['dataset']:<8} {row['model_name']:<20} {row['batch_100k_time_per_seed_ms']:>10.6f}  "
                      f"{speedup:>6.2f}x {status}")
        else:
            print("No batch profiling data available.")
    else:
        print("Batch profiling skipped - test sets too small (14k < 50k min batch)")
        print("For realistic batch profiling, use chr21 dataset (80M k-mers, ~20M test)")

    # Export batch size comparison for leading models
    print(f"\n{'='*70}")
    print("BATCH SIZE COMPARISON (50k-500k)")
    print(f"{'='*70}")

    # Extract batch profile data and create comparison CSV
    batch_comparison_data = []
    if 'batch_profile' in results_df.columns:
        for _, row in results_df.iterrows():
            if pd.notna(row.get('batch_profile')) and row.get('batch_profile'):
                for batch_size, metrics in row['batch_profile'].items():
                    batch_comparison_data.append({
                        'dataset': row['dataset'],
                        'model_name': row['model_name'],
                        'batch_size': batch_size,
                        'time_per_seed_ms': metrics['time_per_seed_ms'],
                        'seeds_per_sec': metrics['seeds_per_sec'],
                        'memory_mb': metrics['memory_mb'],
                        'memory_per_seed_kb': metrics['memory_per_seed_kb'],
                        'vs_bwa_speedup': metrics['vs_bwa_speedup'],
                        'vs_bwa_status': metrics['vs_bwa_status']
                    })

    if batch_comparison_data:
        batch_df = pd.DataFrame(batch_comparison_data)
        batch_output = "results/batch_size_comparison.csv"
        batch_df.to_csv(batch_output, index=False)
        print(f"Batch comparison saved to: {batch_output}")

        # Print comparison for top models
        top_models = results_df.nsmallest(3, 'overfitting_gap_pct')['model_name'].tolist()
        print(f"\nBatch scaling for top 3 models:")
        for model in top_models:
            model_batches = batch_df[batch_df['model_name'] == model]
            if not model_batches.empty:
                print(f"\n{model}:")
                print(f"{'Batch Size':<12} {'ms/seed':<12} {'Seeds/sec':<15} {'Mem (MB)':<10} {'vs BWA':<10}")
                print("-" * 60)
                for _, row in model_batches.iterrows():
                    print(f"{row['batch_size']:<12,} {row['time_per_seed_ms']:<12.6f} "
                          f"{row['seeds_per_sec']:<15,.0f} {row['memory_mb']:<10.2f} "
                          f"{row['vs_bwa_speedup']:<10.2f}x")
    else:
        print("\nNo batch data collected (test sets < 50k minimum batch size)")
        print("Batch profiling will be performed on chr21 dataset (80M k-mers)")

    # SHAP Explainability Analysis for Leading Models
    print(f"\n{'='*70}")
    print("SHAP EXPLAINABILITY ANALYSIS")
    print(f"{'='*70}")

    # Identify top 3 models (lowest overfitting gap with test_f1 > 0.9)
    print("\nIdentifying leading models for explainability analysis...")

    # Sort by overfitting gap and filter by performance
    shap_candidates = []
    for key, model_info in trained_models.items():
        if model_info['test_f1'] > 0.90 and model_info['overfitting_gap'] < 0.15:
            shap_candidates.append(model_info)

    shap_candidates.sort(key=lambda x: (x['overfitting_gap'], -x['test_f1']))

    if shap_candidates:
        print(f"\nFound {len(shap_candidates)} leading models:")
        for i, candidate in enumerate(shap_candidates[:3], 1):
            print(f"  {i}. {candidate['model_name']} ({candidate['dataset']}): "
                  f"F1={candidate['test_f1']:.4f}, Gap={candidate['overfitting_gap']*100:.2f}%")

        # Run SHAP on top 3
        print(f"\nRunning SHAP analysis on top 3 models...")
        for i, candidate in enumerate(shap_candidates[:3], 1):
            dataset = candidate['dataset']
            model_name = candidate['model_name']
            model = candidate['model']
            X_test = candidate['X_test']

            print(f"\n[{i}/3] SHAP: {model_name} on {dataset}")
            shap_output = Path(f"results/shap/{dataset.lower()}_{model_name.lower().replace(' ', '_')}")

            # Save test data temporarily for SHAP
            temp_data_file = shap_output.parent / f"temp_{dataset.lower()}_test.csv"
            temp_data_file.parent.mkdir(parents=True, exist_ok=True)
            X_test.to_csv(temp_data_file, index=False)

            # Determine model type for SHAP explainer
            if 'Tree' in model_name or 'Forest' in model_name or 'LGBM' in model_name or 'XGBoost' in model_name:
                model_type = 'tree'
            elif 'NB' in model_name or 'Centroid' in model_name or 'KNN' in model_name:
                model_type = 'kernel'
            else:
                model_type = 'auto'

            try:
                # Save model temporarily if not already saved
                if not candidate.get('model_path'):
                    temp_model_path = shap_output.parent / f"temp_{dataset.lower()}_{model_name.lower().replace(' ', '_')}.joblib"
                    joblib.dump(model, temp_model_path)
                    model_path = str(temp_model_path)
                else:
                    model_path = candidate['model_path']

                # Run SHAP analysis
                analyze_model_with_shap(
                    model_path=model_path,
                    data_path=str(temp_data_file),
                    output_dir=str(shap_output),
                    sample_size=1000,
                    model_type=model_type
                )
                print(f"  [OK] SHAP analysis complete: {shap_output}/")

                # Clean up temp data
                if temp_data_file.exists():
                    temp_data_file.unlink()

            except Exception as e:
                print(f"  [ERROR] SHAP analysis failed for {model_name} on {dataset}: {e}")

    else:
        print("\nNo models met criteria for SHAP analysis (test_f1 > 0.90, gap < 0.15)")

    # Non-Inferiority Testing vs BWA
    print(f"\n{'='*70}")
    print("NON-INFERIORITY TESTING vs BWA")
    print(f"{'='*70}")
    print("\nClinical Research Framework:")
    print(f"  BWA Reference: {BWA_SEC_PER_SEED*1000:.3f} ms/seed (20-year gold standard)")
    print(f"  Non-Inferiority Margin: 10% ({(BWA_SEC_PER_SEED*1.1)*1000:.3f} ms/seed)")
    print(f"  Significance Level: α = 0.05")

    # Run non-inferiority tests for models with complete batch data
    noninferiority_results = {}

    if batch_comparison_data:
        batch_df = pd.DataFrame(batch_comparison_data)

        # Test each model with sufficient data
        for model_name in batch_df['model_name'].unique():
            model_batches = batch_df[batch_df['model_name'] == model_name]

            # Only test models with all 4 batch sizes
            if len(model_batches) == len(BATCH_SIZES):
                print(f"\n[NON-INFERIORITY] {model_name}")

                try:
                    ni_analysis = comprehensive_noninferiority_analysis(
                        model_batches,
                        model_name,
                        bwa_time=BWA_SEC_PER_SEED
                    )

                    noninferiority_results[model_name] = ni_analysis

                    # Generate and save report
                    output_file = f"results/noninferiority_{model_name.lower().replace(' ', '_')}.txt"
                    report = generate_noninferiority_report(ni_analysis, output_file)

                    # Print summary
                    decision = ni_analysis['overall_decision']
                    print(f"  Non-Inferior: {'YES' if decision['non_inferior'] else 'NO'}")
                    print(f"  Superior: {'YES' if decision['superior'] else 'NO'}")
                    print(f"  Clinical Readiness: {decision['clinical_readiness']}")
                    print(f"  Report saved: {output_file}")

                except Exception as e:
                    print(f"  [ERROR] Non-inferiority test failed: {e}")

            else:
                print(f"\n[SKIP] {model_name} - Insufficient batch data ({len(model_batches)}/{len(BATCH_SIZES)})")

        # Summary table
        if noninferiority_results:
            print(f"\n{'='*70}")
            print("NON-INFERIORITY SUMMARY")
            print(f"{'='*70}")
            print(f"{'Model':<25} {'Non-Inferior':<15} {'Superior':<12} {'Status'}")
            print("-"*70)

            for model_name, analysis in noninferiority_results.items():
                decision = analysis['overall_decision']
                ni_status = 'YES' if decision['non_inferior'] else 'NO'
                sup_status = 'YES' if decision['superior'] else 'NO'
                clinical = decision['clinical_readiness']

                print(f"{model_name:<25} {ni_status:<15} {sup_status:<12} {clinical}")

    else:
        print("\nNo batch data available for non-inferiority testing.")

    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*70}")
    print("\nOutputs:")
    print(f"  - Main results: results/comprehensive_model_analysis.csv")
    print(f"  - Batch comparison: results/batch_size_comparison.csv")
    print(f"  - Non-inferiority reports: results/noninferiority_*.txt")
    print(f"  - Model weights: models/model_weights/")
    print(f"  - SHAP analysis: results/shap/")
    print(f"\nClinical Decision Framework:")
    print(f"  - Non-inferiority: Model ≤ BWA + 10% margin")
    print(f"  - Superiority: Model < BWA (statistically significant)")
    print(f"  - BWA gold standard: 3.000 ms/seed (20-year benchmark)")
    print(f"\nNext steps:")
    print(f"  1. Review non-inferiority reports for clinical deployment readiness")
    print(f"  2. Compare overfitting gaps across S1-S4 for model selection")
    print(f"  3. Examine batch scaling economies (50k-500k)")
    print(f"  4. Analyze SHAP feature importance for leading models")


if __name__ == "__main__":
    main()
