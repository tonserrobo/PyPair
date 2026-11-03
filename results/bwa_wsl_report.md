# BWA-MEM baseline measured in WSL

BWA 0.7.19-r1273 was built from source in WSL Ubuntu 22.04 and run against
synthetic 100 bp reads generated from `data/genomes/chr21.fa` and the four
100 kb samples `data/reference_samples/s1..s4.fasta`.

## Method

- Reads: 10,000 from chr21 (exact and 1% substitution error), 2,000 from each
  100 kb sample (exact). Read IDs and true positions are in `data/bwa_bench/`.
- Three timed runs per configuration with `/usr/bin/time -v`.
- chr21: `bwa mem -t 4` and `bwa mem -t 1` on the exact reads.
- 100 kb samples: `bwa mem -t 4` on exact reads.
- Accuracy: primary SAM alignment, `POS - 1` compared with the truth position;
  exact and within-10 bp rates.
- Summary: `results/bwa_wsl_summary.csv`; raw rows: `results/bwa_wsl_benchmark.csv`.

## Results

| Configuration | Reads | Mapped | Exact position | Within 10 bp | Mean ms/read | Mean max RSS |
|---|---:|---:|---:|---:|---:|---:|
| chr21 exact, 4 threads | 10,000 | 100% | 92.76% | 92.76% | 0.0173 | 91.5 MB |
| chr21 1% error, 4 threads | 10,000 | 100% | 92.57% | 92.65% | 0.0267 | 94.1 MB |
| chr21 exact, 1 thread | 10,000 | 100% | 92.76% | 92.76% | 0.0627 | 90.2 MB |
| S1–S4, 4 threads | 2,000 each | 100% | 100% | 100% | 0.0050 | ~4.5 MB |

The chr21 exact-position rate is below 100% because 100 bp reads in repetitive
or homologous regions can be placed at an alternative but equally valid locus;
the within-10 bp rate is the same as the exact rate, so the misses are not
small offsets.

## What this means for Chapter 3

The chapter's non-inferiority test uses a **3.0 ms/seed whole-genome BWA-MEM
baseline**. The WSL measurements do not reproduce that figure: on a 47 Mb chr21
reference BWA-MEM is ~0.017 ms/read with four threads and ~0.063 ms/read with
one thread, and on a 100 kb reference it is ~0.005 ms/read.

Therefore:

- The 3.0 ms/seed baseline is a whole-genome figure, not a measurement on the
  small references used for the ML samples.
- On chr21, the Python GNB prototype (0.5–0.8 ms/seed) is slower than BWA-MEM
  by roughly an order of magnitude. The KNN prototype (0.12–0.16 ms/seed) is
  within about 2–10× of BWA-MEM depending on thread count.
- The ML speed claim in the chapter depends on the whole-genome baseline. It is
  not supported by the chr21 test. The hardware kernel in Chapter 4 is the
  route to a like-for-like speed comparison on the same reference.

The non-inferiority test should be read with this in mind: the formal result
(32/32 non-inferior at 3.5 ms/seed) holds against the declared whole-genome
baseline, not against the measured chr21 baseline. A separate non-inferiority
run against the chr21 baseline would reject speed non-inferiority for GNB.

## Reproduction

```bash
# In WSL, after building BWA in ~/bwa
python scripts/make_synthetic_reads.py --reference data/genomes/chr21.fa \
  --out-prefix data/bwa_bench/chr21_exact --n-reads 10000 --seed 42
bash scripts/run_bwa_benchmark.sh          # copy the script and data to ~/bwa_bench first
bash scripts/run_bwa_single_thread.sh
```

The raw SAM and `/usr/bin/time` files are kept under `results/bwa_wsl/` and are
not tracked.
