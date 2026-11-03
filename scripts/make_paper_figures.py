#!/usr/bin/env python3
"""Generate the figures used in the condensed PyPair paper from canonical CSVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    }
)

BLUE = "#1f4e79"
LIGHT_BLUE = "#9dc3e6"
ORANGE = "#c55a11"
GREY = "#7f7f7f"
GREEN = "#548235"
RED = "#c00000"


def fig_model_performance(results: Path, out: Path) -> None:
    df = pd.read_csv(results / "canonical_model_comparison.csv")
    df = df.sort_values("test_f1_pct", ascending=False)
    x = range(len(df))
    width = 0.38

    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    ax.bar([i - width / 2 for i in x], df["cv_score_pct"], width, label="CV F1", color=LIGHT_BLUE)
    bars = ax.bar([i + width / 2 for i in x], df["test_f1_pct"], width, label="Held-out F1", color=BLUE)
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3, f"{bar.get_height():.1f}",
                ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["model"], rotation=25, ha="right")
    ax.set_ylabel("Weighted F1 (%)")
    ax.set_ylim(85, 102)
    ax.axhline(95, color=GREY, linestyle=":", linewidth=0.8)
    ax.legend(frameon=False, loc="lower left")
    ax.set_title("Cross-validation and held-out F1 by model")
    fig.savefig(out / "fig1_model_performance.png")
    plt.close(fig)


def fig_overfitting_gap(results: Path, out: Path) -> None:
    df = pd.read_csv(results / "canonical_model_comparison.csv")
    df = df.sort_values("overfitting_gap_pp", ascending=True)
    colors = [GREEN if gap < 1 else ORANGE if gap < 5 else RED for gap in df["overfitting_gap_pp"]]

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    bars = ax.barh(df["model"], df["overfitting_gap_pp"], color=colors)
    for bar in bars:
        ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2, f"{bar.get_width():.1f} pp",
                va="center", fontsize=8)
    ax.set_xlabel("CV–test gap (percentage points)")
    ax.set_xlim(0, max(df["overfitting_gap_pp"]) + 1.8)
    ax.set_title("Generalisation gap: lower is better")
    fig.savefig(out / "fig2_overfitting_gap.png")
    plt.close(fig)


def fig_scaling(results: Path, out: Path) -> None:
    df = pd.read_csv(results / "scaling_analysis.csv")
    df = df[df["model_name"].isin(["KNN", "Gaussian NB"])].copy()
    order = ["S1", "S1+S2", "S1+S2+S3", "S1+S2+S3+S4"]
    labels = ["S1", "S1+S2", "S1–S3", "S1–S4"]
    df["scale"] = pd.Categorical(df["scale"], categories=order, ordered=True)
    df = df.sort_values("scale")

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
    for model, color, marker in [("KNN", ORANGE, "o"), ("Gaussian NB", BLUE, "s")]:
        subset = df[df["model_name"] == model]
        axes[0].plot(range(len(subset)), subset["model_size_mb"], marker=marker, color=color, label=model)
        axes[1].plot(range(len(subset)), subset["train_memory_mb"], marker=marker, color=color, label=model)
    for ax, title, ylabel in [
        (axes[0], "Model size", "Serialized model (MB)"),
        (axes[1], "Training memory", "Peak training memory (MB)"),
    ]:
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False)
    fig.suptitle("Measured scaling from S1 to S1–S4", y=1.02)
    fig.savefig(out / "fig3_scaling.png")
    plt.close(fig)


def fig_bwa_baseline(results: Path, out: Path) -> None:
    bwa = pd.read_csv(results / "bwa_wsl_summary.csv")
    profile = pd.read_csv(results / "seed_mapping_profile_fast.csv")
    ten_k = profile[profile["batch_size"] == 10000].groupby("model_name")["time_per_seed_ms"].mean()

    rows = [
        ("BWA 100 kb (4t)", float(bwa.loc[bwa["group"] == "s1", "ms_per_read_mean"].iloc[0]), BLUE),
        ("BWA chr21 (4t)", float(bwa.loc[bwa["group"] == "chr21_exact", "ms_per_read_mean"].iloc[0]), BLUE),
        ("BWA chr21 1% err (4t)", float(bwa.loc[bwa["group"] == "chr21_err1", "ms_per_read_mean"].iloc[0]), BLUE),
        ("BWA chr21 (1t)", float(bwa.loc[bwa["group"] == "chr21_exact_t1", "ms_per_read_mean"].iloc[0]), BLUE),
        ("KNN (10k batch)", float(ten_k.get("KNN", float("nan"))), ORANGE),
        ("GNB (10k batch)", float(ten_k.get("Gaussian NB", float("nan"))), ORANGE),
        ("Declared whole-genome\nBWA baseline", 3.0, GREY),
    ]
    rows = [row for row in rows if row[1] == row[1]]
    labels = [row[0] for row in rows]
    values = [row[1] for row in rows]
    colors = [row[2] for row in rows]

    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    bars = ax.barh(labels, values, color=colors)
    ax.set_xscale("log")
    ax.set_xlabel("Milliseconds per read/seed (log scale)")
    ax.set_xlim(0.003, 6)
    for bar, value in zip(bars, values):
        ax.text(value * 1.08, bar.get_y() + bar.get_height() / 2, f"{value:.4g}", va="center", fontsize=8)
    ax.set_title("Measured per-read cost: BWA-MEM vs ML seed mapping")
    ax.invert_yaxis()
    fig.savefig(out / "fig4_bwa_baseline.png")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate paper figures from canonical results")
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    fig_model_performance(args.results_dir, args.out_dir)
    fig_overfitting_gap(args.results_dir, args.out_dir)
    fig_scaling(args.results_dir, args.out_dir)
    fig_bwa_baseline(args.results_dir, args.out_dir)
    print(f"wrote 4 figures to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
