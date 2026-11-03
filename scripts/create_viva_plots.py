"""
PyPair: Viva Defense Visualization Generator
=============================================

Creates publication-quality plots for viva defense:
1. Model Performance Comparison (F1 scores)
2. Overfitting Gap Analysis
3. BWA Speedup Comparison
4. Sec/Seed Performance across batch sizes
5. Cross-Dataset Consistency

Usage:
    python scripts/create_viva_plots.py

Output:
    docs/figures/ - All visualization plots
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

# Set style for publication-quality plots
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.5)
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.family'] = 'DejaVu Sans'

# Create output directory
OUTPUT_DIR = Path("docs/figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    """Load validation and profiling results."""
    print("Loading results data...")

    # Model validation results
    validation_df = pd.read_csv("results/comprehensive_model_analysis.csv")

    # Seed mapping profiling results
    profiling_df = pd.read_csv("results/seed_mapping_profile_fast.csv")

    print(f"[OK] Loaded {len(validation_df)} validation results")
    print(f"[OK] Loaded {len(profiling_df)} profiling results")

    return validation_df, profiling_df


def plot_model_performance(validation_df):
    """Plot 1: Model Performance Comparison (F1 Scores)."""
    print("Creating Plot 1: Model Performance Comparison...")

    # Filter successful models
    successful = validation_df[validation_df['test_f1'].notna()].copy()

    # Group by model and average across datasets
    model_summary = successful.groupby('model_name').agg({
        'test_f1': ['mean', 'std'],
        'train_f1': 'mean'
    }).reset_index()

    model_summary.columns = ['model_name', 'test_f1_mean', 'test_f1_std', 'train_f1_mean']
    model_summary = model_summary.sort_values('test_f1_mean', ascending=False)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Bar plot
    x_pos = np.arange(len(model_summary))
    bars = ax.bar(x_pos, model_summary['test_f1_mean'] * 100,
                   yerr=model_summary['test_f1_std'] * 100,
                   capsize=5, alpha=0.8, color='steelblue', edgecolor='black')

    # Highlight top performers
    for i, (idx, row) in enumerate(model_summary.iterrows()):
        if row['test_f1_mean'] > 0.99:
            bars[i].set_color('forestgreen')

    ax.set_xlabel('Model', fontsize=14, fontweight='bold')
    ax.set_ylabel('Test F1 Score (%)', fontsize=14, fontweight='bold')
    ax.set_title('PyPair Model Performance Comparison (S1-S4 Average)',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(model_summary['model_name'], rotation=45, ha='right')
    ax.set_ylim([0, 105])
    ax.axhline(y=99, color='red', linestyle='--', linewidth=2, label='99% Threshold')
    ax.grid(axis='y', alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "01_model_performance.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / '01_model_performance.png'}")
    plt.close()


def plot_overfitting_analysis(validation_df):
    """Plot 2: Overfitting Gap Analysis."""
    print("Creating Plot 2: Overfitting Gap Analysis...")

    # Filter successful models with both train and test F1
    successful = validation_df[
        (validation_df['test_f1'].notna()) &
        (validation_df['train_f1'].notna())
    ].copy()

    successful['overfitting_pct'] = successful['overfitting_gap'] * 100

    # Group by model
    gap_summary = successful.groupby('model_name').agg({
        'overfitting_pct': ['mean', 'std']
    }).reset_index()

    gap_summary.columns = ['model_name', 'gap_mean', 'gap_std']
    gap_summary = gap_summary.sort_values('gap_mean')

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Bar plot
    x_pos = np.arange(len(gap_summary))
    bars = ax.bar(x_pos, gap_summary['gap_mean'],
                   yerr=gap_summary['gap_std'],
                   capsize=5, alpha=0.8, edgecolor='black')

    # Color code: green < 1%, yellow 1-5%, red > 5%
    for i, (idx, row) in enumerate(gap_summary.iterrows()):
        if row['gap_mean'] < 1:
            bars[i].set_color('forestgreen')
        elif row['gap_mean'] < 5:
            bars[i].set_color('gold')
        else:
            bars[i].set_color('tomato')

    ax.set_xlabel('Model', fontsize=14, fontweight='bold')
    ax.set_ylabel('Overfitting Gap (%)', fontsize=14, fontweight='bold')
    ax.set_title('Overfitting Gap: Train F1 - Test F1 (Lower is Better)',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(gap_summary['model_name'], rotation=45, ha='right')
    ax.axhline(y=1, color='green', linestyle='--', linewidth=2, label='< 1% (Minimal)')
    ax.axhline(y=5, color='orange', linestyle='--', linewidth=2, label='< 5% (Acceptable)')
    ax.grid(axis='y', alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "02_overfitting_gap.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / '02_overfitting_gap.png'}")
    plt.close()


def plot_bwa_speedup(profiling_df):
    """Plot 3: BWA Speedup Comparison."""
    print("Creating Plot 3: BWA Speedup Comparison...")

    # Filter for 10k batch size (standard)
    batch_10k = profiling_df[profiling_df['batch_size'] == 10000].copy()

    # Group by model
    speedup_summary = batch_10k.groupby('model_name').agg({
        'speedup_vs_bwa': ['mean', 'std']
    }).reset_index()

    speedup_summary.columns = ['model_name', 'speedup_mean', 'speedup_std']
    speedup_summary = speedup_summary.sort_values('speedup_mean', ascending=False)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Bar plot
    x_pos = np.arange(len(speedup_summary))
    bars = ax.bar(x_pos, speedup_summary['speedup_mean'],
                   yerr=speedup_summary['speedup_std'],
                   capsize=5, alpha=0.8, color='steelblue', edgecolor='black')

    # Highlight models faster than BWA
    for i, (idx, row) in enumerate(speedup_summary.iterrows()):
        if row['speedup_mean'] > 10:
            bars[i].set_color('forestgreen')
        elif row['speedup_mean'] > 1:
            bars[i].set_color('gold')
        else:
            bars[i].set_color('tomato')

    ax.set_xlabel('Model', fontsize=14, fontweight='bold')
    ax.set_ylabel('Speedup vs BWA (x)', fontsize=14, fontweight='bold')
    ax.set_title('PyPair Speedup vs BWA (10k seed batch)',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(speedup_summary['model_name'], rotation=45, ha='right')
    ax.axhline(y=1, color='red', linestyle='--', linewidth=2, label='BWA Baseline (1x)')
    ax.grid(axis='y', alpha=0.3)
    ax.legend()

    # Add value labels on bars
    for i, (idx, row) in enumerate(speedup_summary.iterrows()):
        ax.text(i, row['speedup_mean'] + row['speedup_std'] + 1,
                f"{row['speedup_mean']:.1f}x",
                ha='center', va='bottom', fontweight='bold', fontsize=12)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "03_bwa_speedup.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / '03_bwa_speedup.png'}")
    plt.close()


def plot_batch_performance(profiling_df):
    """Plot 4: Sec/Seed Performance across Batch Sizes."""
    print("Creating Plot 4: Batch Size Performance...")

    # Group by model and batch size
    batch_summary = profiling_df.groupby(['model_name', 'batch_size']).agg({
        'time_per_seed_ms': 'mean'
    }).reset_index()

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Line plot for each model
    for model in batch_summary['model_name'].unique():
        model_data = batch_summary[batch_summary['model_name'] == model]
        ax.plot(model_data['batch_size'], model_data['time_per_seed_ms'],
                marker='o', linewidth=2, markersize=8, label=model)

    # BWA baseline
    ax.axhline(y=3.0, color='red', linestyle='--', linewidth=2, label='BWA (3.0 ms/seed)')

    ax.set_xlabel('Batch Size (seeds)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Time per Seed (ms)', fontsize=14, fontweight='bold')
    ax.set_title('Seed Mapping Performance vs Batch Size',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "04_batch_performance.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / '04_batch_performance.png'}")
    plt.close()


def plot_cross_dataset_consistency(validation_df):
    """Plot 5: Cross-Dataset Consistency (S1-S4)."""
    print("Creating Plot 5: Cross-Dataset Consistency...")

    # Filter successful models
    successful = validation_df[validation_df['test_f1'].notna()].copy()

    # Pivot to show F1 scores across datasets
    consistency_df = successful.pivot_table(
        index='model_name',
        columns='dataset',
        values='test_f1',
        aggfunc='mean'
    )

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Heatmap
    sns.heatmap(consistency_df * 100, annot=True, fmt='.2f', cmap='RdYlGn',
                vmin=0, vmax=100, cbar_kws={'label': 'Test F1 Score (%)'},
                linewidths=0.5, ax=ax)

    ax.set_xlabel('Dataset', fontsize=14, fontweight='bold')
    ax.set_ylabel('Model', fontsize=14, fontweight='bold')
    ax.set_title('Cross-Dataset Consistency: Test F1 Scores (%)',
                 fontsize=16, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "05_cross_dataset_consistency.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / '05_cross_dataset_consistency.png'}")
    plt.close()


def create_summary_infographic():
    """Create summary infographic with key metrics."""
    print("Creating summary infographic...")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('PyPair: Performance Summary', fontsize=20, fontweight='bold', y=0.98)

    # Metric 1: Accuracy (KNN)
    ax1.text(0.5, 0.6, '99.88%', fontsize=80, ha='center', va='center',
             fontweight='bold', color='forestgreen')
    ax1.text(0.5, 0.3, 'Test F1 Score (KNN)', fontsize=20, ha='center', va='center')
    ax1.text(0.5, 0.15, 'Consistent across S1-S4', fontsize=14, ha='center', va='center',
             style='italic')
    ax1.axis('off')

    # Metric 2: Speedup
    ax2.text(0.5, 0.6, '35.7x', fontsize=80, ha='center', va='center',
             fontweight='bold', color='steelblue')
    ax2.text(0.5, 0.3, 'Faster than BWA', fontsize=20, ha='center', va='center')
    ax2.text(0.5, 0.15, '0.084 ms/seed vs 3.0 ms/seed', fontsize=14, ha='center', va='center',
             style='italic')
    ax2.axis('off')

    # Metric 3: Overfitting Gap
    ax3.text(0.5, 0.6, '0.12%', fontsize=80, ha='center', va='center',
             fontweight='bold', color='forestgreen')
    ax3.text(0.5, 0.3, 'Overfitting Gap', fontsize=20, ha='center', va='center')
    ax3.text(0.5, 0.15, 'Minimal generalization gap', fontsize=14, ha='center', va='center',
             style='italic')
    ax3.axis('off')

    # Metric 4: Throughput
    ax4.text(0.5, 0.6, '24k', fontsize=80, ha='center', va='center',
             fontweight='bold', color='steelblue')
    ax4.text(0.5, 0.3, 'Seeds/Second', fontsize=20, ha='center', va='center')
    ax4.text(0.5, 0.15, 'Peak throughput (1k batch)', fontsize=14, ha='center', va='center',
             style='italic')
    ax4.axis('off')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "00_summary_infographic.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / '00_summary_infographic.png'}")
    plt.close()


def main():
    """Generate all viva defense plots."""
    print("="*70)
    print("PyPair: Viva Defense Visualization Generator")
    print("="*70)

    # Load data
    validation_df, profiling_df = load_data()

    # Create plots
    create_summary_infographic()
    plot_model_performance(validation_df)
    plot_overfitting_analysis(validation_df)
    plot_bwa_speedup(profiling_df)
    plot_batch_performance(profiling_df)
    plot_cross_dataset_consistency(validation_df)

    print("="*70)
    print("ALL PLOTS CREATED")
    print("="*70)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Total plots: 6")
    print("="*70)


if __name__ == "__main__":
    main()
