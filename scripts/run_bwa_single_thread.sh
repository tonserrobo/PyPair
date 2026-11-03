#!/usr/bin/env bash
set -euo pipefail

cd "$HOME/bwa_bench"
for i in 1 2 3; do
  /usr/bin/time -v "$HOME/bwa/bwa" mem -t 1 chr21.fa chr21_exact_reads.fastq \
    > "chr21_exact_t1_run${i}.sam" 2> "chr21_exact_t1_time${i}.txt"
  echo "t1 run ${i} done"
done

cp chr21_exact_t1_run*.sam chr21_exact_t1_time*.txt \
  "/mnt/c/Users/Tony/OneDrive - Ulster University/Projects/PyPair/results/bwa_wsl/"
