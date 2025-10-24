#!/usr/bin/env python3
"""
Add extended features to chr21 k-mer dataset

Adds the same features used in S1-S4 training:
- Nucleotide counts (A, C, G, T, first)
- Reverse complement encoding (rc_seeds)
- Run-length encoding features (runlength_*)

This makes chr21 compatible with models trained on S1-S4.

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import sys
import time
import vaex
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter

# Configuration
INPUT_FILE = Path("data/reference_genomes/chr21_kmers.hdf5")
OUTPUT_FILE = Path("data/reference_genomes/chr21_kmers_extended.hdf5")
BATCH_SIZE = 1_000_000  # Process 1M rows at a time


def nucleotide_counts(kmer_series):
    """
    Compute nucleotide counts from k-mer string.

    Returns dict with: first, A, C, G, T
    """
    results = {
        'first': [],
        'A': [],
        'C': [],
        'G': [],
        'T': []
    }

    for kmer in kmer_series:
        if pd.isna(kmer) or len(kmer) == 0:
            results['first'].append(0)
            results['A'].append(0)
            results['C'].append(0)
            results['G'].append(0)
            results['T'].append(0)
        else:
            # First nucleotide (encoded: A=0, C=1, G=2, T=3)
            first_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
            results['first'].append(first_map.get(kmer[0], 0))

            # Count nucleotides
            counts = Counter(kmer)
            results['A'].append(counts.get('A', 0))
            results['C'].append(counts.get('C', 0))
            results['G'].append(counts.get('G', 0))
            results['T'].append(counts.get('T', 0))

    return results


def reverse_complement(kmer):
    """Get reverse complement of k-mer."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if pd.isna(kmer) or len(kmer) == 0:
        return ""
    return ''.join(complement.get(base, 'N') for base in reversed(kmer))


def encode_kmer_binary(kmer):
    """
    Encode k-mer as binary integer (2-bit encoding).
    A=00, C=01, G=10, T=11
    """
    if pd.isna(kmer) or len(kmer) == 0:
        return 0

    encoding = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    binary_str = ''.join(encoding.get(base, '00') for base in kmer)

    # Limit to 63 bits to fit in int64
    binary_str = binary_str[:63]

    try:
        return int(binary_str, 2)
    except:
        return 0


def run_length_encode(sequence):
    """
    Run-length encoding of sequence.

    Returns: output (position), value (count), character
    """
    if pd.isna(sequence) or len(sequence) == 0:
        return 0, 0, 0

    runs = []
    current_char = sequence[0]
    current_count = 1

    for i in range(1, len(sequence)):
        if sequence[i] == current_char:
            current_count += 1
        else:
            runs.append((current_char, current_count))
            current_char = sequence[i]
            current_count = 1

    runs.append((current_char, current_count))

    # Return first run info (most compressed)
    if runs:
        char_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        return 0, runs[0][1], char_map.get(runs[0][0], 0)
    else:
        return 0, 0, 0


def process_batch_features(batch_df):
    """
    Add extended features to a batch of k-mers.

    Input: DataFrame with 'kmer' column
    Returns: DataFrame with all extended features
    """
    print(f"  Processing batch of {len(batch_df):,} k-mers...")

    # Nucleotide counts
    print("    Computing nucleotide counts...")
    nuc_counts = nucleotide_counts(batch_df['kmer'])
    batch_df['first'] = nuc_counts['first']
    batch_df['A'] = nuc_counts['A']
    batch_df['C'] = nuc_counts['C']
    batch_df['G'] = nuc_counts['G']
    batch_df['T'] = nuc_counts['T']

    # Reverse complement
    print("    Computing reverse complement...")
    batch_df['rc_kmer'] = batch_df['kmer'].apply(reverse_complement)
    batch_df['rc_seeds'] = batch_df['rc_kmer'].apply(encode_kmer_binary)

    # Run-length encoding (forward)
    print("    Computing run-length encoding (forward)...")
    rle_forward = batch_df['kmer'].apply(run_length_encode)
    batch_df['runlength_output'] = [x[0] for x in rle_forward]
    batch_df['runlength_value'] = [x[1] for x in rle_forward]

    # Run-length encoding (reverse complement)
    print("    Computing run-length encoding (RC)...")
    rle_rc = batch_df['rc_kmer'].apply(run_length_encode)
    batch_df['runlength_rc_output'] = [x[0] for x in rle_rc]
    batch_df['runlength_rc_value'] = [x[1] for x in rle_rc]
    batch_df['runlength_rc_character'] = [x[2] for x in rle_rc]

    # Drop temporary columns
    batch_df = batch_df.drop(columns=['rc_kmer'])

    return batch_df


def main():
    print("="*80)
    print("chr21 Extended Feature Engineering")
    print("="*80)
    print(f"Input: {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    print(f"Batch size: {BATCH_SIZE:,}")
    print("="*80)

    # Check input exists
    if not INPUT_FILE.exists():
        print(f"\nERROR: Input file not found: {INPUT_FILE}")
        print("Please run chr21 k-mer extraction first.")
        sys.exit(1)

    # Load chr21 data
    print("\nLoading chr21 k-mers...")
    df_vaex = vaex.open(str(INPUT_FILE))
    total_rows = len(df_vaex)
    print(f"Total k-mers: {total_rows:,}")
    print(f"Current columns: {list(df_vaex.columns)}")

    # Process in batches
    print(f"\nProcessing in batches of {BATCH_SIZE:,}...")
    start_time = time.time()

    processed_batches = []
    num_batches = (total_rows + BATCH_SIZE - 1) // BATCH_SIZE

    for batch_idx in range(num_batches):
        batch_start = batch_idx * BATCH_SIZE
        batch_end = min(batch_start + BATCH_SIZE, total_rows)

        print(f"\nBatch {batch_idx + 1}/{num_batches}: rows {batch_start:,} - {batch_end:,}")

        # Extract batch as pandas
        batch_vaex = df_vaex[batch_start:batch_end]
        batch_df = batch_vaex.to_pandas_df()

        # Add extended features
        batch_df = process_batch_features(batch_df)

        # Convert back to vaex
        batch_vaex_extended = vaex.from_pandas(batch_df)
        processed_batches.append(batch_vaex_extended)

        # Progress update
        elapsed = time.time() - start_time
        rows_done = batch_end
        rows_per_sec = rows_done / elapsed if elapsed > 0 else 0
        eta_sec = (total_rows - rows_done) / rows_per_sec if rows_per_sec > 0 else 0

        print(f"  Progress: {rows_done:,}/{total_rows:,} ({rows_done/total_rows*100:.1f}%)")
        print(f"  Speed: {rows_per_sec:,.0f} rows/sec")
        print(f"  ETA: {eta_sec/60:.1f} minutes")

    # Concatenate all batches
    print(f"\nConcatenating {len(processed_batches)} batches...")
    df_extended = vaex.concat(processed_batches)

    print(f"\nExtended dataset:")
    print(f"  Rows: {len(df_extended):,}")
    print(f"  Columns: {list(df_extended.columns)}")

    # Save to HDF5
    print(f"\nSaving to {OUTPUT_FILE}...")
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    df_extended.export_hdf5(str(OUTPUT_FILE))

    total_time = time.time() - start_time

    print("\n" + "="*80)
    print("FEATURE ENGINEERING COMPLETE!")
    print("="*80)
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Output file: {OUTPUT_FILE}")
    print(f"File size: {OUTPUT_FILE.stat().st_size / (1024**3):.2f} GB")
    print("\nExtended features added:")
    print("  - first, A, C, G, T (nucleotide counts)")
    print("  - rc_seeds (reverse complement encoding)")
    print("  - runlength_* (RLE compression features)")
    print("\nReady for training with models trained on S1-S4!")
    print("="*80)


if __name__ == "__main__":
    main()
