# Results provenance

This note records which result files are safe to cite and why the repository
contains several apparently conflicting numbers. It is the reference for the
values in Chapter 3 and for `results/canonical_model_comparison.csv`.

## Canonical results

| Value | Canonical source |
|---|---|
| Eight-model aggregate comparison (CV, F1, gap, size, latency) | `results/canonical_model_comparison.csv` (transcribed from Chapter 3 Table 3.6) |
| GNB per-dataset F1 and resources | `results/NaiveBayes_results.csv` |
| Nearest Centroid per-dataset F1 and resources | `results/nearest_centroid_results.csv` |
| S1→S4 scaling | `results/scaling_analysis.csv` |
| Throughput (32 configurations) | `results/seed_mapping_profile_fast.csv` |
| Non-inferiority tests | `results/noninferiority_results.csv`, `results/noninferiority_report.md` |
| Protocol B cleaned saved-model evaluation | `results/cleaned_model_evaluation.csv` |
| BWA-MEM baseline (WSL) | `results/bwa_wsl_summary.csv`, `results/bwa_wsl_report.md` |

The canonical files above are tracked in git. The large SHAP artefacts and the
contents of `results/legacy/` are ignored.

The canonical aggregate table is:

| Model | CV % | Test F1 % | Gap pp | Train s | Size MB | ms/seed |
|---|---:|---:|---:|---:|---:|---:|
| KNN | 99.7 | 94.4 | 5.3 | 12 | 2.0 | 0.084 |
| Nearest Centroid | 98.0 | 95.7 | 2.3 | 8 | 1.2 | 0.012 |
| Gaussian NB | 95.1 | 94.7 | 0.4 | 15 | 7.0 | 0.84 |
| Decision Tree | 96.2 | 94.8 | 1.4 | 45 | 0.025 | 0.003 |
| Extra Trees | 99.2 | 94.9 | 4.3 | 320 | 118 | 1.2 |
| Random Forest | 100.0 | 92.8 | 7.2 | 280 | 117 | 1.1 |
| XGBoost | 100.0 | 91.8 | 8.2 | 450 | 102 | 2.3 |
| LightGBM | 100.0 | 90.5 | 9.5 | 380 | 233 | 1.8 |

## Why the files diverge

The repository contains results from several protocol versions:

1. **Protocol A — aggregate eight-model run** (`docs/overfitting_analysis.md`, `results/*_results.csv`). This is the protocol-consistent baseline used in Chapter 3. It yields the table above with KNN 94.4.
2. **Protocol B — cleaned saved-model evaluation** (`results/cleaned_model_evaluation.csv`). The saved `*_s1..s4.joblib` models were evaluated on the cleaned test sets: KNN 99.87%, GNB 95.08–95.51%. This is a coherent protocol, but it is not the aggregate table. The 99.88% KNN figure in the draft belongs to Protocol B, not to the aggregate run; it is a protocol mismatch, not a transcription error.
3. **`S1_model_comparison.csv`** (11 features, 10,000 classes, 44,800/14,000 split). This run reports GNB 95.08%, KNN 84.75%, Nearest Centroid 99.87% on S1. It is a later refactor with a different feature matrix and failure handling. It is not directly comparable with the aggregate table.
4. **`comprehensive_model_analysis.csv`**. This run records mixed outcomes: Extra Trees fails with `too many values to unpack`, Random Forest has no test F1, and KNN and Nearest Centroid both report 99.87%. The KNN values duplicate Nearest Centroid, so the KNN column is not independent.
5. **`chr21_model_comparison.csv`**. Every model is a fatal error (`'DataFrameLocal' object has no attribute 'to_numpy'`). No chromosome-scale training completed.

The old `results/KNeighborsClassifier_results.csv` was a stale copy of the GNB metrics: its precision, recall and F1 were identical to `NaiveBayes_results.csv` to four decimal places. It has been moved to `results/legacy/`. The clean replacement is `results/cleaned_model_evaluation.csv`, which evaluates the saved KNN and GNB models directly and shows KNN 99.87% and GNB 95.08–95.51%.

The GNB model size also has two legitimate values: 7 MB in the aggregate table and 1.9 MB for the per-sample serialized model in `scaling_analysis.csv`. The difference is the number of classes and the serialization format. State which one is meant.

## BWA-MEM baseline (WSL)

BWA 0.7.19-r1273 was built and run in WSL. On the 47 Mb chr21 reference, BWA-MEM
averaged 0.0173 ms/read with four threads and 0.0627 ms/read with one thread
(10,000 exact 100 bp reads, three runs). On the 100 kb samples it averaged
0.0050 ms/read. The mapped rate was 100%; the exact-position rate was 92.8% on
chr21 (alternative loci in repeats) and 100% on the 100 kb samples.

The chapter's 3.0 ms/seed BWA-MEM baseline is a **whole-genome** figure. It is
not reproduced on chr21 or on the 100 kb samples. The formal non-inferiority
result (32/32 at a 3.5 ms/seed threshold) therefore holds against the declared
whole-genome baseline, not against the measured chr21 baseline. Against the
chr21 baseline the Python GNB prototype is slower than BWA-MEM. See
`results/bwa_wsl_report.md` and `results/bwa_wsl_summary.csv`.

## What Chapter 3 uses

- Table 3.6: `results/canonical_model_comparison.csv`.
- Table 3.7: `results/seed_mapping_profile_fast.csv` (10,000-seed batch and all-batch means).
- Table 3.8: `results/scaling_analysis.csv`.
- Table 3.10: `results/NaiveBayes_results.csv`.
- Table 3.11 and Section 4.8: `results/noninferiority_results.csv` and `results/noninferiority_report.md`.

Sections 4.2.4 and 5.5 of Chapter 3 explicitly flag the provenance problem. Do not replace those caveats with a single number until a canonical rerun exists.

## Non-inferiority method

`scripts/noninferiority_testing.py` now computes the one-sided one-sample
t-test from summary statistics (mean, sample SD, n = 5 profiler runs). It does
not simulate replicates. The default margins are 0.5 ms/seed for speed and
2 percentage points for accuracy, matching Chapter 3.

Outputs:

- `results/noninferiority_results.csv`: one row per model × dataset × batch.
- `results/noninferiority_report.md`: method, speed table, model summary, accuracy test and caveats.

Caveats:

- The profiler stores summary statistics over five runs, not raw run-level timings. The t-test is exact for the summary mean and sample SD, but cannot detect run-level outliers.
- The accuracy test treats the four per-dataset F1 values as replicates and the BWA accuracy (96.2%) as a fixed reference. A paired test needs BWA accuracy on the same four samples.
- The huge Cohen's d values in the CSV reflect the distance between the observed mean and the threshold, not a meaningful practical scale. Report direction and p-value, not the magnitude.

## Recommended reconciliation

1. Freeze one configuration: feature set, seed length, class construction, split, random seeds, scikit-learn version.
2. Re-run the eight models and write a single `results/canonical_model_comparison.csv` with a `commit` and `run_id` column.
3. Re-run `profile_seed_mapping_fast.py` with raw replicate timings saved (not only summary statistics).
4. Re-run the chromosome-21 Vaex path after fixing the `DataFrameLocal`/`to_numpy` error; keep the failure file only as a legacy artefact.
5. Update Chapter 3 and `results/README.md` in the same commit.
