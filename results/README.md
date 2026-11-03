# Results: provenance and canonical set

The `results/` directory contains files from several protocol versions. Use the
files listed under **Canonical** for any reported number. Do not mix values
across files.

## Canonical

| File | Produced by | Used for |
|---|---|---|
| `canonical_model_comparison.csv` | Chapter 3 Table 3.6 | Eight-model aggregate comparison (CV, F1, gap, size, latency) |
| `NaiveBayes_results.csv` | `models/classifiers/gaussian_nb.py` | GNB per-dataset F1, memory and timing |
| `nearest_centroid_results.csv` | `models/classifiers/nearest_centroid.py` | Nearest Centroid per-dataset F1, memory and timing |
| `scaling_analysis.csv` | `scripts/scaling_analysis.py` | Measured S1→S4 model size, memory and prediction time |
| `seed_mapping_profile_fast.csv` | `scripts/profile_seed_mapping_fast.py` | 32 timing profiles (KNN/GNB × S1–S4 × four batch sizes) |
| `seed_mapping_summary_fast.txt` | `scripts/profile_seed_mapping_fast.py` | Human-readable throughput summary |
| `noninferiority_results.csv` | `scripts/noninferiority_testing.py` | Per-configuration non-inferiority tests |
| `noninferiority_report.md` | `scripts/noninferiority_testing.py` | Method, speed table, accuracy test, caveats |
| `cleaned_model_evaluation.csv` | `scripts/evaluate_saved_models.py` | Protocol B: saved KNN/GNB evaluated on cleaned test sets |
| `bwa_wsl_summary.csv` | `scripts/parse_bwa_sam.py` | Measured BWA-MEM baseline in WSL |
| `bwa_wsl_report.md` | WSL BWA benchmark | Method, results and the whole-genome baseline caveat |

Protocol A is the aggregate eight-model table (`canonical_model_comparison.csv`).
Protocol B is the cleaned saved-model evaluation (`cleaned_model_evaluation.csv`).
The KNN 99.88% value belongs to Protocol B; the aggregate KNN value is 94.4%.
Do not mix them.

The BWA summary is a measured chr21/100 kb baseline. The chapter's 3.0 ms/seed
figure is a whole-genome baseline; see `bwa_wsl_report.md`.

## Legacy (moved to `results/legacy/`)

| File | Why it is not canonical |
|---|---|
| `KNeighborsClassifier_results.csv` | Stale copy of the GNB metrics; clean replacement is `cleaned_model_evaluation.csv` (KNN 99.87% on Protocol B) |
| `comprehensive_model_analysis.csv` | Mixed protocol: several models fail (`too many values to unpack`, `DataFrameLocal`), and KNN values duplicate Nearest Centroid |
| `S1_model_comparison.csv` | Different feature/protocol run (11 features, 10k classes); KNN 84.75% and Nearest Centroid 99.87% do not match the aggregate table |
| `chr21_model_comparison.csv` | Every model is a fatal error (`'DataFrameLocal' object has no attribute 'to_numpy'`); no training completed |

The canonical files in the table above are tracked in git; the SHAP artefacts,
the raw BWA SAM/time files (`results/bwa_wsl/`) and `results/legacy/` are ignored.
The Gaussian NB model weights are tracked with Git LFS (`git lfs pull`).

## Reproduction

```bash
# Throughput profile (canonical timing source)
python scripts/profile_seed_mapping_fast.py

# Scaling analysis
python scripts/scaling_analysis.py

# Non-inferiority tests (uses the profile and GNB accuracy above)
python scripts/noninferiority_testing.py \
  --profile results/seed_mapping_profile_fast.csv \
  --accuracy results/NaiveBayes_results.csv \
  --margin-ms 0.5 --accuracy-margin-pp 2.0 --bwa-accuracy 96.2
```

## Rule

Every reported value should name the CSV and the git commit that produced it.
If a result file is regenerated, update `docs/results_provenance.md` and the
chapter table in the same commit.
