#!/usr/bin/env bash
set -euo pipefail

cd "$HOME/bwa_bench"
BWA="$HOME/bwa/bwa"

for label in chr21_exact chr21_err1; do
  for i in 1 2 3; do
    /usr/bin/time -v "$BWA" mem -t 4 chr21.fa "${label}_reads.fastq" \
      > "${label}_run${i}.sam" 2> "${label}_time${i}.txt"
    echo "$label run $i done"
  done
done

for s in s1 s2 s3 s4; do
  "$BWA" index "${s}.fasta" > "${s}_index.log" 2>&1
  for i in 1 2 3; do
    /usr/bin/time -v "$BWA" mem -t 4 "${s}.fasta" "${s}_reads.fastq" \
      > "${s}_run${i}.sam" 2> "${s}_time${i}.txt"
    echo "$s run $i done"
  done
done
