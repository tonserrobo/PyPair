#!/usr/bin/env python3
"""
SHAP (SHapley Additive exPlanations) Analysis for PyPair Models

This script performs comprehensive SHAP analysis on trained models to understand
feature importance and model explainability.

Usage:
    python shap_analysis.py \\
        --model models/model_weights/gnb_chr21.joblib \\
        --data data/reference_genomes/chr21_kmers.hdf5 \\
        --output results/shap/gnb_chr21 \\
        --samples 1000

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import argparse
import os
import sys
from pathlib import Path
import joblib
import vaex
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shap

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


def analyze_model_with_shap(
    model_path: str,
    data_path: str,
    output_dir: str,
    sample_size: int = 1000,
    model_type: str = "auto"
):
    """
    Perform comprehensive SHAP analysis on trained model.

    Args:
        model_path: Path to trained model (.joblib)
        data_path: Path to test data (.hdf5 or .csv)
        output_dir: Directory for output plots
        sample_size: Number of samples for SHAP analysis
        model_type: Model type (auto, tree, linear, kernel)
    """
    print("="*60)
    print("PyPair: SHAP Explainability Analysis")
    print("="*60)
    print(f"Model: {model_path}")
    print(f"Data: {data_path}")
    print(f"Output: {output_dir}")
    print(f"Sample size: {sample_size}")
    print("="*60)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Load model
    print("\n[1/6] Loading model...")
    model = joblib.load(model_path)
    print(f"  Model type: {type(model).__name__}")

    # Load data
    print("\n[2/6] Loading data...")
    if data_path.endswith('.hdf5'):
        df = vaex.open(data_path)
    elif data_path.endswith('.csv'):
        df = pd.read_csv(data_path)
    else:
        raise ValueError(f"Unsupported file type: {data_path}")

    # Prepare feature matrix
    print("\n[3/6] Preparing features...")
    # Drop non-feature columns
    feature_cols = [col for col in df.columns if col not in ['chromosome', 'position', 'kmer', 'strand']]

    if isinstance(df, vaex.dataframe.DataFrame):
        # Sample from Vaex
        X_sample = df[feature_cols].sample(n=min(sample_size, len(df))).to_pandas_df()
    else:
        # Sample from pandas
        X_sample = df[feature_cols].sample(n=min(sample_size, len(df)))

    print(f"  Sample shape: {X_sample.shape}")
    print(f"  Features: {list(X_sample.columns)}")

    # Create SHAP explainer
    print("\n[4/6] Creating SHAP explainer...")

    # Detect model type if auto
    if model_type == "auto":
        if hasattr(model, 'tree_') or 'Tree' in type(model).__name__ or 'Forest' in type(model).__name__:
            model_type = "tree"
        elif hasattr(model, 'coef_'):
            model_type = "linear"
        else:
            model_type = "kernel"

    print(f"  Explainer type: {model_type}")

    if model_type == "tree":
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_sample)
    elif model_type == "linear":
        explainer = shap.LinearExplainer(model, X_sample)
        shap_values = explainer.shap_values(X_sample)
    else:  # kernel explainer (model-agnostic, slower)
        background = shap.sample(X_sample, min(100, len(X_sample)))
        explainer = shap.KernelExplainer(model.predict, background)
        shap_values = explainer.shap_values(X_sample)

    print(f"  SHAP values shape: {np.array(shap_values).shape}")

    # Generate plots
    print("\n[5/6] Generating SHAP visualizations...")

    # Plot 1: Summary plot (global feature importance)
    print("  - Summary plot...")
    plt.figure(figsize=(12, 8))
    shap.summary_plot(shap_values, X_sample, show=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/shap_summary.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 2: Bar plot (mean absolute SHAP values)
    print("  - Feature importance bar plot...")
    plt.figure(figsize=(10, 8))
    shap.plots.bar(
        shap.Explanation(
            values=np.abs(shap_values).mean(0),
            feature_names=X_sample.columns
        ),
        show=False
    )
    plt.tight_layout()
    plt.savefig(f"{output_dir}/shap_importance.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 3: Beeswarm plot
    print("  - Beeswarm plot...")
    plt.figure(figsize=(12, 8))
    try:
        shap.plots.beeswarm(
            shap.Explanation(
                values=shap_values,
                data=X_sample.values,
                feature_names=X_sample.columns
            ),
            show=False
        )
        plt.tight_layout()
        plt.savefig(f"{output_dir}/shap_beeswarm.png", dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"    Warning: Beeswarm plot failed: {e}")

    # Plot 4: Waterfall plot (single prediction example)
    print("  - Waterfall plot (first sample)...")
    plt.figure(figsize=(10, 8))
    try:
        shap.plots.waterfall(
            shap.Explanation(
                values=shap_values[0],
                base_values=explainer.expected_value if hasattr(explainer, 'expected_value') else 0,
                data=X_sample.iloc[0],
                feature_names=X_sample.columns
            ),
            show=False
        )
        plt.tight_layout()
        plt.savefig(f"{output_dir}/shap_waterfall_sample0.png", dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"    Warning: Waterfall plot failed: {e}")

    # Plot 5: Dependence plots for top 3 features
    print("  - Dependence plots (top 3 features)...")
    feature_importance = np.abs(shap_values).mean(0)
    top_features_idx = np.argsort(feature_importance)[-3:][::-1]

    for idx in top_features_idx:
        feature_name = X_sample.columns[idx]
        plt.figure(figsize=(10, 6))
        try:
            shap.dependence_plot(
                idx,
                shap_values,
                X_sample,
                show=False
            )
            plt.tight_layout()
            plt.savefig(f"{output_dir}/shap_dependence_{feature_name}.png", dpi=300, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"    Warning: Dependence plot for {feature_name} failed: {e}")

    # Save data
    print("\n[6/6] Saving SHAP data...")
    np.save(f"{output_dir}/shap_values.npy", shap_values)
    X_sample.to_csv(f"{output_dir}/shap_sample_data.csv", index=False)

    # Compute and save feature importance
    feature_importance_df = pd.DataFrame({
        'feature': X_sample.columns,
        'mean_abs_shap': np.abs(shap_values).mean(0),
        'mean_shap': shap_values.mean(0),
        'std_shap': shap_values.std(0)
    }).sort_values('mean_abs_shap', ascending=False)

    feature_importance_df.to_csv(f"{output_dir}/feature_importance.csv", index=False)

    print("\n" + "="*60)
    print("SHAP Analysis Complete!")
    print("="*60)
    print(f"\nTop 10 Most Important Features:")
    print(feature_importance_df.head(10).to_string(index=False))
    print(f"\nAll outputs saved to: {output_dir}/")
    print("="*60)

    return feature_importance_df


def main():
    parser = argparse.ArgumentParser(
        description="SHAP explainability analysis for PyPair models",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Analyze GNB model on chr21
    python shap_analysis.py \\
        --model models/model_weights/gnb_chr21.joblib \\
        --data data/reference_genomes/chr21_kmers.hdf5 \\
        --output results/shap/gnb_chr21 \\
        --samples 1000

    # Analyze Random Forest (tree-based)
    python shap_analysis.py \\
        --model models/model_weights/random_forest_s1_s4.joblib \\
        --data data/reference_samples/S1_extended_encoded_X_test.csv \\
        --output results/shap/rf_baseline \\
        --model-type tree
        """
    )

    parser.add_argument(
        "--model",
        required=True,
        help="Path to trained model file (.joblib)"
    )
    parser.add_argument(
        "--data",
        required=True,
        help="Path to test data (.hdf5 or .csv)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for SHAP plots and data"
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=1000,
        help="Number of samples to use for SHAP analysis (default: 1000)"
    )
    parser.add_argument(
        "--model-type",
        choices=["auto", "tree", "linear", "kernel"],
        default="auto",
        help="Type of SHAP explainer (default: auto)"
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.model):
        print(f"ERROR: Model file not found: {args.model}")
        sys.exit(1)

    if not os.path.exists(args.data):
        print(f"ERROR: Data file not found: {args.data}")
        sys.exit(1)

    # Run analysis
    feature_importance = analyze_model_with_shap(
        model_path=args.model,
        data_path=args.data,
        output_dir=args.output,
        sample_size=args.samples,
        model_type=args.model_type
    )


if __name__ == "__main__":
    main()
