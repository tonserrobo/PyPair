"""
PyPair: Scaling Analysis Visualization
=======================================

Creates plots showing KNN vs GNB performance at different data scales.

Usage:
    python scripts/create_scaling_plots.py

Output:
    docs/figures/scaling_*.png - Scaling analysis plots
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

# Set style
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.5)
plt.rcParams['figure.figsize'] = (12, 8)

OUTPUT_DIR = Path("docs/figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_scaling_data():
    """Load scaling analysis results."""
    df = pd.read_csv("results/scaling_analysis.csv")
    print(f"[OK] Loaded {len(df)} scaling results")
    return df


def plot_memory_scaling(df):
    """Plot memory usage vs data scale."""
    fig, ax = plt.subplots(figsize=(12, 8))

    for model in df['model_name'].unique():
        model_data = df[df['model_name'] == model]
        ax.plot(model_data['train_samples'], model_data['model_size_mb'],
                marker='o', linewidth=2, markersize=10, label=model)

    ax.set_xlabel('Training Samples', fontsize=14, fontweight='bold')
    ax.set_ylabel('Model Size (MB)', fontsize=14, fontweight='bold')
    ax.set_title('Model Size Scaling: KNN vs GNB',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "scaling_01_memory.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / 'scaling_01_memory.png'}")
    plt.close()


def plot_memory_per_sample(df):
    """Plot memory per sample (shows KNN's O(n) vs GNB's O(1))."""
    fig, ax = plt.subplots(figsize=(12, 8))

    for model in df['model_name'].unique():
        model_data = df[df['model_name'] == model]
        ax.plot(model_data['train_samples'], model_data['model_size_per_sample_kb'],
                marker='o', linewidth=2, markersize=10, label=model)

    ax.set_xlabel('Training Samples', fontsize=14, fontweight='bold')
    ax.set_ylabel('Model Size per Sample (KB)', fontsize=14, fontweight='bold')
    ax.set_title('Memory Efficiency: KNN (O(n)) vs GNB (O(1))',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xscale('log')
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "scaling_02_memory_per_sample.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / 'scaling_02_memory_per_sample.png'}")
    plt.close()


def plot_training_time_scaling(df):
    """Plot training time vs data scale."""
    fig, ax = plt.subplots(figsize=(12, 8))

    for model in df['model_name'].unique():
        model_data = df[df['model_name'] == model]
        ax.plot(model_data['train_samples'], model_data['train_time_sec'],
                marker='o', linewidth=2, markersize=10, label=model)

    ax.set_xlabel('Training Samples', fontsize=14, fontweight='bold')
    ax.set_ylabel('Training Time (seconds)', fontsize=14, fontweight='bold')
    ax.set_title('Training Time Scaling: KNN vs GNB',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "scaling_03_training_time.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / 'scaling_03_training_time.png'}")
    plt.close()


def plot_efficiency_comparison(df):
    """Side-by-side comparison at each scale."""
    scales = df['scale'].unique()
    n_scales = len(scales)

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    # Memory comparison
    ax1 = axes[0]
    width = 0.35
    x = np.arange(n_scales)

    knn_sizes = []
    gnb_sizes = []

    for scale in scales:
        scale_data = df[df['scale'] == scale]
        knn_sizes.append(scale_data[scale_data['model_name'] == 'KNN']['model_size_mb'].values[0])
        gnb_sizes.append(scale_data[scale_data['model_name'] == 'Gaussian NB']['model_size_mb'].values[0])

    ax1.bar(x - width/2, knn_sizes, width, label='KNN', color='steelblue')
    ax1.bar(x + width/2, gnb_sizes, width, label='GNB', color='forestgreen')

    ax1.set_xlabel('Data Scale', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Model Size (MB)', fontsize=12, fontweight='bold')
    ax1.set_title('Model Size Comparison', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(scales, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # Time comparison
    ax2 = axes[1]
    knn_times = []
    gnb_times = []

    for scale in scales:
        scale_data = df[df['scale'] == scale]
        knn_times.append(scale_data[scale_data['model_name'] == 'KNN']['train_time_sec'].values[0])
        gnb_times.append(scale_data[scale_data['model_name'] == 'Gaussian NB']['train_time_sec'].values[0])

    ax2.bar(x - width/2, knn_times, width, label='KNN', color='steelblue')
    ax2.bar(x + width/2, gnb_times, width, label='GNB', color='forestgreen')

    ax2.set_xlabel('Data Scale', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Training Time (sec)', fontsize=12, fontweight='bold')
    ax2.set_title('Training Time Comparison', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(scales, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "scaling_04_comparison.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / 'scaling_04_comparison.png'}")
    plt.close()


def create_projection_plot(df):
    """Extrapolate to genome scale (3B k-mers)."""
    fig, ax = plt.subplots(figsize=(12, 8))

    # Fit power law to existing data
    for model in df['model_name'].unique():
        model_data = df[df['model_name'] == model].copy()

        # Fit log-log linear regression (power law)
        log_samples = np.log10(model_data['train_samples'].values)
        log_size = np.log10(model_data['model_size_mb'].values)

        # Fit
        coeffs = np.polyfit(log_samples, log_size, 1)
        slope, intercept = coeffs

        # Extrapolate to 3B
        projection_samples = np.logspace(4, 9.5, 100)  # 10k to 3B
        projection_size = 10 ** (slope * np.log10(projection_samples) + intercept)

        # Plot actual data
        ax.scatter(model_data['train_samples'], model_data['model_size_mb'],
                   s=100, label=f'{model} (actual)', zorder=3)

        # Plot projection
        ax.plot(projection_samples, projection_size,
                linestyle='--', linewidth=2, label=f'{model} (projected)')

    # Highlight genome scale
    ax.axvline(x=3e9, color='red', linestyle=':', linewidth=2, label='Human Genome (3B k-mers)')

    # Highlight memory limits
    ax.axhline(y=1024, color='orange', linestyle=':', linewidth=2, label='1TB Memory Limit')
    ax.axhline(y=256, color='yellow', linestyle=':', linewidth=2, label='256GB Memory (typical HPC node)')

    ax.set_xlabel('Training Samples (k-mers)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Model Size (MB)', fontsize=14, fontweight='bold')
    ax.set_title('Projected Memory Requirements: Genome Scale (3B k-mers)',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "scaling_05_genome_projection.png", dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {OUTPUT_DIR / 'scaling_05_genome_projection.png'}")
    plt.close()


def main():
    """Generate all scaling plots."""
    print("="*70)
    print("PyPair: Scaling Analysis Visualization")
    print("="*70)

    df = load_scaling_data()

    plot_memory_scaling(df)
    plot_memory_per_sample(df)
    plot_training_time_scaling(df)
    plot_efficiency_comparison(df)
    create_projection_plot(df)

    print("="*70)
    print("ALL SCALING PLOTS CREATED")
    print("="*70)
    print(f"Output: {OUTPUT_DIR}")
    print("="*70)


if __name__ == "__main__":
    main()
