#!/usr/bin/env python3
"""
Non-Inferiority Testing for ML-based Seed Mapping vs BWA

Clinical research framework for validating that PyPair models are
non-inferior to BWA (the 20-year gold standard).

Non-inferiority criterion: ML model inference time ≤ BWA + margin
Margin: 10% (industry standard for clinical non-inferiority)

Statistical tests:
1. One-sided t-test (H0: ML >= BWA + margin, H1: ML < BWA + margin)
2. 95% confidence interval (upper bound must be < BWA + margin)
3. Effect size (Cohen's d)

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import argparse
import math
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from scipy import stats


# BWA Benchmark (from literature and user specification)
BWA_SEC_PER_SEED = 0.003  # 3 milliseconds per seed
BWA_THROUGHPUT = 1 / BWA_SEC_PER_SEED  # ~333 seeds/sec

# Non-inferiority margin (10% slower is acceptable)
NONINFERIORITY_MARGIN = 0.10  # 10%
NONINFERIORITY_THRESHOLD = BWA_SEC_PER_SEED * (1 + NONINFERIORITY_MARGIN)  # 3.3 ms

# Clinical significance threshold
ALPHA = 0.05  # 95% confidence
POWER = 0.80  # 80% power


def noninferiority_ttest(
    ml_times: np.ndarray,
    bwa_time: float = BWA_SEC_PER_SEED,
    margin: float = NONINFERIORITY_MARGIN,
    alpha: float = ALPHA
) -> Dict:
    """
    One-sided t-test for non-inferiority.

    H0: μ_ML >= μ_BWA + margin  (ML is inferior)
    H1: μ_ML < μ_BWA + margin   (ML is non-inferior)

    Args:
        ml_times: Array of ML model inference times (sec/seed)
        bwa_time: BWA reference time (sec/seed)
        margin: Non-inferiority margin (proportion)
        alpha: Significance level

    Returns:
        dict with test statistics and decision
    """
    threshold = bwa_time * (1 + margin)

    # One-sample t-test against threshold
    t_stat, p_value = stats.ttest_1samp(
        ml_times,
        threshold,
        alternative='less'  # Testing if ML < threshold
    )

    # Confidence interval (one-sided upper bound)
    mean_ml = np.mean(ml_times)
    se_ml = stats.sem(ml_times)
    df = len(ml_times) - 1
    t_critical = stats.t.ppf(1 - alpha, df)
    ci_upper = mean_ml + t_critical * se_ml

    # Effect size (Cohen's d)
    std_ml = np.std(ml_times, ddof=1)
    cohens_d = (threshold - mean_ml) / std_ml if std_ml > 0 else np.inf

    # Decision: Non-inferior if p < alpha AND CI upper bound < threshold
    non_inferior = (p_value < alpha) and (ci_upper < threshold)

    return {
        'ml_mean_sec': mean_ml,
        'ml_std_sec': std_ml,
        'bwa_threshold_sec': threshold,
        't_statistic': t_stat,
        'p_value': p_value,
        'ci_upper_bound': ci_upper,
        'cohens_d': cohens_d,
        'non_inferior': non_inferior,
        'decision': 'NON-INFERIOR' if non_inferior else 'INCONCLUSIVE',
        'speedup_vs_bwa': bwa_time / mean_ml if mean_ml > 0 else np.inf,
        'percent_faster': ((bwa_time - mean_ml) / bwa_time * 100) if mean_ml > 0 else np.inf
    }


def noninferiority_ttest_from_summary(
    mean: float,
    std: float,
    n: int,
    bwa_time: float = BWA_SEC_PER_SEED,
    margin: float = NONINFERIORITY_MARGIN,
    alpha: float = ALPHA,
) -> Dict:
    """One-sample one-sided t-test from summary statistics.

    This is the exact test for a mean against a non-inferiority threshold; it
    does not simulate replicates. ``std`` is the sample standard deviation
    (ddof=1).
    """
    if n < 2:
        raise ValueError("n must be at least 2 for a t-test")

    threshold = bwa_time * (1 + margin)
    df = n - 1
    se = std / math.sqrt(n)
    t_stat = (mean - threshold) / se if se > 0 else (-math.inf if mean < threshold else math.inf)
    p_value = float(stats.t.cdf(t_stat, df))
    t_critical = float(stats.t.ppf(1 - alpha, df))
    ci_upper = mean + t_critical * se
    cohens_d = (threshold - mean) / std if std > 0 else math.inf
    non_inferior = bool(p_value < alpha and ci_upper < threshold)
    superior = bool(p_value < alpha and ci_upper < bwa_time)

    return {
        'ml_mean_sec': mean,
        'ml_std_sec': std,
        'bwa_threshold_sec': threshold,
        'threshold': threshold,
        'bwa_time_sec': bwa_time,
        'n': n,
        'standard_error': se,
        't_statistic': t_stat,
        'p_value': p_value,
        'ci_upper_bound': ci_upper,
        'cohens_d': cohens_d,
        'non_inferior': non_inferior,
        'superior': superior,
        'decision': 'NON-INFERIOR' if non_inferior else 'INCONCLUSIVE',
        'speedup_vs_bwa': bwa_time / mean if mean > 0 else math.inf,
        'percent_faster': ((bwa_time - mean) / bwa_time * 100) if mean > 0 else math.inf,
    }


def superiority_ttest(
    ml_times: np.ndarray,
    bwa_time: float = BWA_SEC_PER_SEED,
    alpha: float = ALPHA
) -> Dict:
    """
    One-sided t-test for superiority (ML faster than BWA).

    H0: μ_ML >= μ_BWA  (ML is not faster)
    H1: μ_ML < μ_BWA   (ML is faster)

    Args:
        ml_times: Array of ML model inference times (sec/seed)
        bwa_time: BWA reference time (sec/seed)
        alpha: Significance level

    Returns:
        dict with test statistics and decision
    """
    t_stat, p_value = stats.ttest_1samp(
        ml_times,
        bwa_time,
        alternative='less'
    )

    mean_ml = np.mean(ml_times)
    se_ml = stats.sem(ml_times)
    df = len(ml_times) - 1
    t_critical = stats.t.ppf(1 - alpha, df)
    ci_upper = mean_ml + t_critical * se_ml

    superior = (p_value < alpha) and (ci_upper < bwa_time)

    return {
        'ml_mean_sec': mean_ml,
        'bwa_time_sec': bwa_time,
        't_statistic': t_stat,
        'p_value': p_value,
        'ci_upper_bound': ci_upper,
        'superior': superior,
        'decision': 'SUPERIOR' if superior else 'NOT SUPERIOR',
        'speedup': bwa_time / mean_ml if mean_ml > 0 else np.inf
    }


def comprehensive_noninferiority_analysis(
    batch_results: pd.DataFrame,
    model_name: str,
    bwa_time: float = BWA_SEC_PER_SEED,
    margin: float = NONINFERIORITY_MARGIN
) -> Dict:
    """
    Comprehensive non-inferiority analysis for a model across batch sizes.

    Args:
        batch_results: DataFrame with columns [batch_size, time_per_seed_sec, ...]
        model_name: Name of the model
        bwa_time: BWA reference time
        margin: Non-inferiority margin

    Returns:
        dict with comprehensive statistics
    """
    results = {
        'model_name': model_name,
        'bwa_reference_sec': bwa_time,
        'noninferiority_threshold_sec': bwa_time * (1 + margin),
        'margin_percent': margin * 100,
        'batch_analyses': {}
    }

    # Analyze each batch size
    for _, row in batch_results.iterrows():
        batch_size = int(row['batch_size'])
        n_runs = int(row.get('n_runs', 5))

        if 'mean_time_sec' in row and 'std_time_sec' in row:
            mean_time = float(row['mean_time_sec']) / batch_size
            std_population = float(row['std_time_sec']) / batch_size
        else:
            mean_time = float(row['time_per_seed_sec'])
            std_population = float(row.get('time_std_sec', mean_time * 0.05))

        std_time = (
            std_population * math.sqrt(n_runs / (n_runs - 1))
            if n_runs > 1 else 0.0
        )

        ni_test = noninferiority_ttest_from_summary(
            mean=mean_time,
            std=std_time,
            n=n_runs,
            bwa_time=bwa_time,
            margin=margin,
        )

        results['batch_analyses'][batch_size] = {
            'mean_time_sec': mean_time,
            'std_time_sec': std_time,
            'noninferiority': ni_test,
            'superiority': {
                'p_value': ni_test['p_value'],
                'superior': ni_test['superior'],
                'decision': 'SUPERIOR' if ni_test['superior'] else 'NOT SUPERIOR',
            },
            'memory_mb': row.get('memory_mb', None),
            'speedup': ni_test['speedup_vs_bwa'],
        }

    # Overall decision across all batch sizes
    all_non_inferior = all(
        analysis['noninferiority']['non_inferior']
        for analysis in results['batch_analyses'].values()
    )

    all_superior = all(
        analysis['superiority']['superior']
        for analysis in results['batch_analyses'].values()
    )

    results['overall_decision'] = {
        'non_inferior': all_non_inferior,
        'superior': all_superior,
        'clinical_readiness': 'READY' if all_non_inferior else 'NOT READY',
        'recommendation': (
            'DEPLOY - Superior to BWA' if all_superior else
            'DEPLOY - Non-inferior to BWA' if all_non_inferior else
            'DO NOT DEPLOY - Inferior to BWA'
        )
    }

    return results


def generate_noninferiority_report(
    analysis_results: Dict,
    output_file: str = None
) -> str:
    """
    Generate a clinical-style non-inferiority report.

    Args:
        analysis_results: Output from comprehensive_noninferiority_analysis
        output_file: Optional file path to save report

    Returns:
        Formatted report string
    """
    model = analysis_results['model_name']
    bwa = analysis_results['bwa_reference_sec'] * 1000  # Convert to ms
    threshold = analysis_results['noninferiority_threshold_sec'] * 1000
    margin = analysis_results['margin_percent']

    report = []
    report.append("="*80)
    report.append("NON-INFERIORITY ANALYSIS REPORT")
    report.append("="*80)
    report.append(f"\nModel: {model}")
    report.append(f"Reference Standard: BWA ({bwa:.3f} ms/seed)")
    report.append(f"Non-Inferiority Threshold: {threshold:.3f} ms/seed")
    report.append(f"Margin: {margin:.1f}% (Clinical Standard)")
    report.append(f"Significance Level: α = {ALPHA}")
    report.append("\n" + "-"*80)
    report.append("BATCH SIZE ANALYSIS")
    report.append("-"*80)

    for batch_size, analysis in analysis_results['batch_analyses'].items():
        ni = analysis['noninferiority']
        sup = analysis['superiority']

        report.append(f"\nBatch Size: {batch_size:,} seeds")
        report.append(f"  Mean Time: {ni['ml_mean_sec']*1000:.6f} ms/seed")
        report.append(f"  Std Dev: {ni['ml_std_sec']*1000:.6f} ms/seed")
        report.append(f"  Speedup: {ni['speedup_vs_bwa']:.2f}x")
        report.append(f"  Percent Faster: {ni['percent_faster']:.1f}%")
        report.append(f"\n  Non-Inferiority Test:")
        report.append(f"    t-statistic: {ni['t_statistic']:.4f}")
        report.append(f"    p-value: {ni['p_value']:.6f}")
        report.append(f"    95% CI Upper: {ni['ci_upper_bound']*1000:.6f} ms/seed")
        report.append(f"    Cohen's d: {ni['cohens_d']:.4f}")
        report.append(f"    Decision: {ni['decision']}")
        report.append(f"\n  Superiority Test:")
        report.append(f"    p-value: {sup['p_value']:.6f}")
        report.append(f"    Decision: {sup['decision']}")
        if analysis.get('memory_mb'):
            report.append(f"  Memory: {analysis['memory_mb']:.2f} MB")

    report.append("\n" + "="*80)
    report.append("OVERALL CLINICAL DECISION")
    report.append("="*80)

    overall = analysis_results['overall_decision']
    report.append(f"\nNon-Inferior to BWA: {'YES' if overall['non_inferior'] else 'NO'}")
    report.append(f"Superior to BWA: {'YES' if overall['superior'] else 'NO'}")
    report.append(f"Clinical Readiness: {overall['clinical_readiness']}")
    report.append(f"\nRECOMMENDATION: {overall['recommendation']}")

    report.append("\n" + "="*80)
    report.append("CLINICAL INTERPRETATION")
    report.append("="*80)

    if overall['superior']:
        report.append("\nThis model demonstrates SUPERIORITY over BWA across all batch sizes.")
        report.append("Statistical significance (p < 0.05) and clinical significance (effect size)")
        report.append("confirm that this model is meaningfully faster than the gold standard.")
        report.append("\nCONCLUSION: SUITABLE FOR CLINICAL DEPLOYMENT")

    elif overall['non_inferior']:
        report.append("\nThis model demonstrates NON-INFERIORITY to BWA across all batch sizes.")
        report.append("Performance is within the 10% margin, confirming that this model can")
        report.append("replace BWA without clinically meaningful performance degradation.")
        report.append("\nCONCLUSION: SUITABLE FOR CLINICAL DEPLOYMENT")

    else:
        report.append("\nThis model FAILS to demonstrate non-inferiority to BWA.")
        report.append("Performance exceeds the 10% margin, indicating clinically meaningful")
        report.append("degradation compared to the gold standard.")
        report.append("\nCONCLUSION: NOT SUITABLE FOR CLINICAL DEPLOYMENT")

    report.append("\n" + "="*80)
    report.append(f"Generated: 2025-01-24")
    report.append(f"Author: Tony Robinson (t.robinson@ulster.ac.uk)")
    report.append("="*80)

    report_text = "\n".join(report)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_text)

    return report_text


def _profile_tests(
    profile: pd.DataFrame,
    bwa_time: float,
    margin: float,
    alpha: float,
    n_runs: int,
) -> pd.DataFrame:
    """Run the summary-statistics test for every profiled configuration."""
    rows = []
    for _, row in profile.iterrows():
        batch_size = int(row["batch_size"])
        mean = float(row["mean_time_sec"]) / batch_size
        std_population = float(row["std_time_sec"]) / batch_size
        std = (
            std_population * math.sqrt(n_runs / (n_runs - 1))
            if n_runs > 1 else 0.0
        )
        result = noninferiority_ttest_from_summary(
            mean=mean,
            std=std,
            n=n_runs,
            bwa_time=bwa_time,
            margin=margin,
            alpha=alpha,
        )
        rows.append({
            "model": row["model_name"],
            "dataset": row["dataset"],
            "batch_size": batch_size,
            "n_runs": n_runs,
            "mean_ms_per_seed": mean * 1000,
            "std_ms_per_seed": std * 1000,
            "threshold_ms": result["threshold"] * 1000,
            "t_statistic": result["t_statistic"],
            "p_value": result["p_value"],
            "ci_upper_ms": result["ci_upper_bound"] * 1000,
            "cohens_d": result["cohens_d"],
            "speedup_vs_bwa": result["speedup_vs_bwa"],
            "non_inferior": result["non_inferior"],
            "superior": result["superior"],
            "decision": result["decision"],
        })
    return pd.DataFrame(rows)


def _accuracy_test(
    accuracy: pd.DataFrame,
    bwa_accuracy: float,
    margin_pp: float,
    alpha: float,
) -> Dict:
    """Test GNB accuracy non-inferiority by transforming F1 to a loss."""
    if "F1 Score" not in accuracy.columns:
        raise ValueError("accuracy file must contain an 'F1 Score' column")

    f1 = accuracy["F1 Score"].astype(float) * 100
    mean_f1 = float(f1.mean())
    std_f1 = float(f1.std(ddof=1))
    n = int(len(f1))
    if n < 2:
        raise ValueError("accuracy file needs at least two rows")

    bwa_loss = 100.0 - bwa_accuracy
    mean_loss = 100.0 - mean_f1
    margin = margin_pp / bwa_loss
    result = noninferiority_ttest_from_summary(
        mean=mean_loss,
        std=std_f1,
        n=n,
        bwa_time=bwa_loss,
        margin=margin,
        alpha=alpha,
    )
    return {
        "n_datasets": n,
        "mean_f1": mean_f1,
        "std_f1": std_f1,
        "bwa_accuracy": bwa_accuracy,
        "margin_pp": margin_pp,
        "threshold_f1": 100.0 - result["threshold"],
        "mean_loss": mean_loss,
        "p_value": result["p_value"],
        "ci_upper_loss": result["ci_upper_bound"],
        "ci_lower_f1": 100.0 - result["ci_upper_bound"],
        "cohens_d": result["cohens_d"],
        "non_inferior": result["non_inferior"],
        "decision": result["decision"],
    }


def _write_report(
    speed: pd.DataFrame,
    accuracy: Dict | None,
    bwa_ms: float,
    margin_ms: float,
    alpha: float,
    n_runs: int,
    path: Path,
) -> None:
    lines = [
        "# Non-inferiority results",
        "",
        "## Method",
        "",
        f"- BWA-MEM baseline: {bwa_ms:.3f} ms/seed.",
        f"- Non-inferiority margin: {margin_ms:.3f} ms (threshold {bwa_ms + margin_ms:.3f} ms/seed).",
        f"- n = {n_runs} replicate timings per configuration, from `results/seed_mapping_profile_fast.csv`.",
        "- One-sided one-sample t-test computed exactly from the summary mean and sample SD; no simulated replicates.",
        f"- Decision: non-inferior if p < {alpha:.2f} and the one-sided 95% upper CI is below the threshold.",
        "",
        "## Speed",
        "",
        "| Model | Dataset | Batch | Mean ms/seed | 95% CI upper | Threshold | p | Speed-up | Decision |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for _, row in speed.iterrows():
        p_display = "< 0.001" if row['p_value'] < 0.001 else f"{row['p_value']:.3g}"
        lines.append(
            f"| {row['model']} | {row['dataset']} | {row['batch_size']:,} | "
            f"{row['mean_ms_per_seed']:.6f} | {row['ci_upper_ms']:.6f} | "
            f"{row['threshold_ms']:.3f} | {p_display} | "
            f"{row['speedup_vs_bwa']:.2f}x | {row['decision']} |"
        )

    lines += ["", "## Model summary", "",
              "| Model | Configurations | Non-inferior | Superior | Mean speed-up |",
              "|---|---:|---:|---:|---:|"]
    for model, group in speed.groupby("model"):
        lines.append(
            f"| {model} | {len(group)} | {int(group['non_inferior'].sum())} | "
            f"{int(group['superior'].sum())} | {group['speedup_vs_bwa'].mean():.2f}x |"
        )

    if accuracy:
        lines += [
            "",
            "## Accuracy (Gaussian NB)",
            "",
            f"- Mean per-dataset F1: {accuracy['mean_f1']:.2f}% (SD {accuracy['std_f1']:.2f}, n = {accuracy['n_datasets']}).",
            f"- BWA-MEM reference accuracy: {accuracy['bwa_accuracy']:.2f}%.",
            f"- Non-inferiority threshold: {accuracy['threshold_f1']:.2f}% (margin {accuracy['margin_pp']:.2f} pp).",
            f"- One-sided 95% lower F1 bound: {accuracy['ci_lower_f1']:.2f}%.",
            f"- p-value (loss < threshold): {accuracy['p_value']:.3g}.",
            f"- Decision: **{accuracy['decision']}**.",
        ]

    lines += [
        "",
        "## Caveats",
        "",
        "- The profiler retains summary statistics over five runs, not the raw run-level timings. "
        "The t-test uses the summary mean and sample SD, which is exact for a one-sample test "
        "but cannot detect run-level outliers.",
        "- The accuracy test treats the four per-dataset F1 values as replicates and the BWA "
        "accuracy as a fixed reference. A paired test would require BWA accuracy on the same "
        "four samples.",
        "- The model results in the repository come from more than one protocol; this report uses "
        "`seed_mapping_profile_fast.csv` and `NaiveBayes_results.csv` only. See "
        "`docs/results_provenance.md`.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Non-inferiority tests for PyPair vs BWA-MEM")
    parser.add_argument("--profile", type=Path, default=Path("results/seed_mapping_profile_fast.csv"))
    parser.add_argument("--accuracy", type=Path, default=Path("results/NaiveBayes_results.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/noninferiority_results.csv"))
    parser.add_argument("--report", type=Path, default=Path("results/noninferiority_report.md"))
    parser.add_argument("--bwa-ms", type=float, default=3.0)
    parser.add_argument("--margin-ms", type=float, default=0.5)
    parser.add_argument("--accuracy-margin-pp", type=float, default=2.0)
    parser.add_argument("--bwa-accuracy", type=float, default=96.2)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--n-runs", type=int, default=5)
    args = parser.parse_args(argv)

    if not args.profile.is_file():
        raise FileNotFoundError(f"profile not found: {args.profile}")

    profile = pd.read_csv(args.profile)
    required = {"model_name", "dataset", "batch_size", "mean_time_sec", "std_time_sec"}
    missing = required - set(profile.columns)
    if missing:
        raise ValueError(f"profile missing columns: {sorted(missing)}")

    bwa_time = args.bwa_ms / 1000
    margin = (args.margin_ms / 1000) / bwa_time
    speed = _profile_tests(profile, bwa_time, margin, args.alpha, args.n_runs)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    speed.to_csv(args.out, index=False)

    accuracy = None
    if args.accuracy.is_file():
        accuracy = _accuracy_test(
            pd.read_csv(args.accuracy),
            args.bwa_accuracy,
            args.accuracy_margin_pp,
            args.alpha,
        )

    _write_report(speed, accuracy, args.bwa_ms, args.margin_ms, args.alpha, args.n_runs, args.report)

    print(f"speed rows: {len(speed)}")
    print(f"non-inferior: {int(speed['non_inferior'].sum())}/{len(speed)}")
    print(f"superior: {int(speed['superior'].sum())}/{len(speed)}")
    if accuracy:
        print(f"accuracy: {accuracy['decision']} (p={accuracy['p_value']:.3g})")
    print(f"wrote {args.out}")
    print(f"wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
