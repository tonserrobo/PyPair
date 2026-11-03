# Non-inferiority results

## Method

- BWA-MEM baseline: 3.000 ms/seed.
- Non-inferiority margin: 0.500 ms (threshold 3.500 ms/seed).
- n = 5 replicate timings per configuration, from `results/seed_mapping_profile_fast.csv`.
- One-sided one-sample t-test computed exactly from the summary mean and sample SD; no simulated replicates.
- Decision: non-inferior if p < 0.05 and the one-sided 95% upper CI is below the threshold.

## Speed

| Model | Dataset | Batch | Mean ms/seed | 95% CI upper | Threshold | p | Speed-up | Decision |
|---|---|---:|---:|---:|---:|---:|---:|---|
| KNN | S1 | 1,000 | 0.042240 | 0.043464 | 3.500 | < 0.001 | 71.02x | NON-INFERIOR |
| KNN | S1 | 5,000 | 0.062592 | 0.064542 | 3.500 | < 0.001 | 47.93x | NON-INFERIOR |
| KNN | S1 | 10,000 | 0.083964 | 0.085519 | 3.500 | < 0.001 | 35.73x | NON-INFERIOR |
| KNN | S1 | 14,000 | 0.090726 | 0.091665 | 3.500 | < 0.001 | 33.07x | NON-INFERIOR |
| Gaussian NB | S1 | 1,000 | 0.383310 | 0.393392 | 3.500 | < 0.001 | 7.83x | NON-INFERIOR |
| Gaussian NB | S1 | 5,000 | 0.327593 | 0.329637 | 3.500 | < 0.001 | 9.16x | NON-INFERIOR |
| Gaussian NB | S1 | 10,000 | 0.360880 | 0.397932 | 3.500 | < 0.001 | 8.31x | NON-INFERIOR |
| Gaussian NB | S1 | 14,000 | 1.358873 | 1.542208 | 3.500 | < 0.001 | 2.21x | NON-INFERIOR |
| KNN | S2 | 1,000 | 0.056807 | 0.075975 | 3.500 | < 0.001 | 52.81x | NON-INFERIOR |
| KNN | S2 | 5,000 | 0.107938 | 0.147362 | 3.500 | < 0.001 | 27.79x | NON-INFERIOR |
| KNN | S2 | 10,000 | 0.154030 | 0.189698 | 3.500 | < 0.001 | 19.48x | NON-INFERIOR |
| KNN | S2 | 14,000 | 0.176895 | 0.200948 | 3.500 | < 0.001 | 16.96x | NON-INFERIOR |
| Gaussian NB | S2 | 1,000 | 0.990066 | 1.016711 | 3.500 | < 0.001 | 3.03x | NON-INFERIOR |
| Gaussian NB | S2 | 5,000 | 0.724957 | 0.735121 | 3.500 | < 0.001 | 4.14x | NON-INFERIOR |
| Gaussian NB | S2 | 10,000 | 0.722663 | 0.777236 | 3.500 | < 0.001 | 4.15x | NON-INFERIOR |
| Gaussian NB | S2 | 14,000 | 1.400645 | 1.430555 | 3.500 | < 0.001 | 2.14x | NON-INFERIOR |
| KNN | S3 | 1,000 | 0.097929 | 0.125779 | 3.500 | < 0.001 | 30.63x | NON-INFERIOR |
| KNN | S3 | 5,000 | 0.117409 | 0.157611 | 3.500 | < 0.001 | 25.55x | NON-INFERIOR |
| KNN | S3 | 10,000 | 0.131042 | 0.168232 | 3.500 | < 0.001 | 22.89x | NON-INFERIOR |
| KNN | S3 | 14,000 | 0.166299 | 0.204080 | 3.500 | < 0.001 | 18.04x | NON-INFERIOR |
| Gaussian NB | S3 | 1,000 | 0.995971 | 1.016397 | 3.500 | < 0.001 | 3.01x | NON-INFERIOR |
| Gaussian NB | S3 | 5,000 | 0.529128 | 0.734647 | 3.500 | < 0.001 | 5.67x | NON-INFERIOR |
| Gaussian NB | S3 | 10,000 | 0.373638 | 0.461567 | 3.500 | < 0.001 | 8.03x | NON-INFERIOR |
| Gaussian NB | S3 | 14,000 | 1.431286 | 1.439701 | 3.500 | < 0.001 | 2.10x | NON-INFERIOR |
| KNN | S4 | 1,000 | 0.061516 | 0.079383 | 3.500 | < 0.001 | 48.77x | NON-INFERIOR |
| KNN | S4 | 5,000 | 0.135835 | 0.170639 | 3.500 | < 0.001 | 22.09x | NON-INFERIOR |
| KNN | S4 | 10,000 | 0.138884 | 0.179994 | 3.500 | < 0.001 | 21.60x | NON-INFERIOR |
| KNN | S4 | 14,000 | 0.190968 | 0.197743 | 3.500 | < 0.001 | 15.71x | NON-INFERIOR |
| Gaussian NB | S4 | 1,000 | 0.966360 | 0.979510 | 3.500 | < 0.001 | 3.10x | NON-INFERIOR |
| Gaussian NB | S4 | 5,000 | 0.730104 | 0.763421 | 3.500 | < 0.001 | 4.11x | NON-INFERIOR |
| Gaussian NB | S4 | 10,000 | 0.702256 | 0.733992 | 3.500 | < 0.001 | 4.27x | NON-INFERIOR |
| Gaussian NB | S4 | 14,000 | 1.402245 | 1.408503 | 3.500 | < 0.001 | 2.14x | NON-INFERIOR |

## Model summary

| Model | Configurations | Non-inferior | Superior | Mean speed-up |
|---|---:|---:|---:|---:|
| Gaussian NB | 16 | 16 | 16 | 4.59x |
| KNN | 16 | 16 | 16 | 31.88x |

## Accuracy (Gaussian NB)

- Mean per-dataset F1: 94.72% (SD 0.07, n = 4).
- BWA-MEM reference accuracy: 96.20%.
- Non-inferiority threshold: 94.20% (margin 2.00 pp).
- One-sided 95% lower F1 bound: 94.63%.
- p-value (loss < threshold): 0.000392.
- Decision: **NON-INFERIOR**.

## Caveats

- The profiler retains summary statistics over five runs, not the raw run-level timings. The t-test uses the summary mean and sample SD, which is exact for a one-sample test but cannot detect run-level outliers.
- The accuracy test treats the four per-dataset F1 values as replicates and the BWA accuracy as a fixed reference. A paired test would require BWA accuracy on the same four samples.
- The model results in the repository come from more than one protocol; this report uses `seed_mapping_profile_fast.csv` and `NaiveBayes_results.csv` only. See `docs/results_provenance.md`.
