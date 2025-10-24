#!/usr/bin/env python3
"""
Extract k-mers from specific chromosomes and create HDF5 training data.

This script parses FASTA files for specified chromosomes, generates k-mers,
and creates memory-mapped HDF5 files for Vaex-based out-of-core training.

Usage:
    python extract_chromosome_data.py \
        --fasta data/genomes/chr21.fa \
        --output data/reference_genomes/chr21_kmers.hdf5 \
        --k 15 \
        --chromosomes chr21

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List, Dict
import pandas as pd
import vaex
from Bio import SeqIO
from tqdm import tqdm

# Add project root to path for imports
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


def encode_kmer(kmer: str) -> int:
    """
    Encode DNA k-mer to integer using 2-bit encoding.

    A=00, C=01, G=10, T=11

    Args:
        kmer: DNA sequence string

    Returns:
        Integer encoding of k-mer
    """
    encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    result = 0
    for base in kmer:
        result = (result << 2) | encoding.get(base.upper(), 0)
    return result


def reverse_complement(seq: str) -> str:
    """Compute reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq.upper()))


def extract_kmers_from_chromosome_batched(
    sequence: str,
    chrom_name: str,
    k: int = 15,
    include_rc: bool = True,
    batch_size: int = 1_000_000
):
    """
    Extract all k-mers from a chromosome sequence in batches (generator).

    This avoids loading all k-mers into memory at once by yielding batches
    of DataFrames that can be written incrementally to HDF5.

    Args:
        sequence: DNA sequence string
        chrom_name: Chromosome name (e.g., 'chr21')
        k: k-mer size
        include_rc: Include reverse complement
        batch_size: Number of k-mers to process before yielding a batch

    Yields:
        DataFrame batches with columns: chromosome, position, kmer, encoded_kmer, gc_content, strand
    """
    data = []
    seq_len = len(sequence)
    total_kmers = 0

    print(f"  Extracting k-mers (k={k}) from {chrom_name} ({seq_len:,} bp)...")
    print(f"  Processing in batches of {batch_size:,} k-mers to manage memory...")

    for i in tqdm(range(seq_len - k + 1), desc=f"  {chrom_name}"):
        kmer = sequence[i:i+k]

        # Skip k-mers with ambiguous bases
        if 'N' in kmer.upper():
            continue

        # Calculate GC content
        gc_count = kmer.upper().count('G') + kmer.upper().count('C')
        gc_content = gc_count / k

        # Encode k-mer
        encoded = encode_kmer(kmer)

        # Add forward k-mer
        data.append({
            'chromosome': chrom_name,
            'position': i,
            'kmer': kmer,
            'encoded_kmer': encoded,
            'gc_content': gc_content,
            'strand': 'forward'
        })

        # Add reverse complement if requested
        if include_rc:
            rc_kmer = reverse_complement(kmer)
            rc_encoded = encode_kmer(rc_kmer)
            data.append({
                'chromosome': chrom_name,
                'position': i,
                'kmer': rc_kmer,
                'encoded_kmer': rc_encoded,
                'gc_content': gc_content,
                'strand': 'reverse'
            })

        # Yield batch when size is reached
        if len(data) >= batch_size:
            df_batch = pd.DataFrame(data)
            total_kmers += len(df_batch)
            yield df_batch
            data = []  # Clear for next batch

    # Yield remaining data
    if data:
        df_batch = pd.DataFrame(data)
        total_kmers += len(df_batch)
        yield df_batch

    print(f"  Extracted {total_kmers:,} k-mers from {chrom_name}")


def extract_chromosome_data(
    fasta_files: List[str],
    chromosomes: List[str],
    output_file: str,
    k: int = 15,
    include_rc: bool = True,
    batch_size: int = 1_000_000
):
    """
    Extract k-mers from specified chromosomes and save to HDF5.

    Uses batched processing to avoid loading all k-mers into memory at once.

    Args:
        fasta_files: List of FASTA file paths
        chromosomes: List of chromosome names to extract (e.g., ['chr21', 'chr22'])
        output_file: Output HDF5 file path
        k: k-mer size
        include_rc: Include reverse complement k-mers
        batch_size: Number of k-mers to process per batch (default: 1 million)
    """
    print("="*60)
    print("PyPair: Chromosome k-mer Extraction")
    print("="*60)
    print(f"Input FASTA files: {len(fasta_files)}")
    print(f"Target chromosomes: {chromosomes}")
    print(f"k-mer size: {k}")
    print(f"Include reverse complement: {include_rc}")
    print(f"Batch size: {batch_size:,} k-mers")
    print(f"Output: {output_file}")
    print("="*60)

    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Collect all batch DataFrames
    all_batches = []
    total_rows = 0

    # Parse each FASTA file
    for fasta_file in fasta_files:
        print(f"\nProcessing {fasta_file}...")

        if not os.path.exists(fasta_file):
            print(f"  ERROR: File not found: {fasta_file}")
            continue

        # Read FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            chrom_name = record.id

            # Check if this chromosome should be extracted
            if chromosomes and chrom_name not in chromosomes:
                print(f"  Skipping {chrom_name} (not in target list)")
                continue

            # Extract k-mers in batches
            sequence = str(record.seq)
            for batch_df in extract_kmers_from_chromosome_batched(
                sequence, chrom_name, k, include_rc, batch_size
            ):
                # Convert to Vaex and collect
                vdf_batch = vaex.from_pandas(batch_df)
                all_batches.append(vdf_batch)
                total_rows += len(vdf_batch)

                print(f"    Processed batch: {len(vdf_batch):,} rows (total: {total_rows:,})")

    if total_rows == 0:
        print("\nERROR: No k-mers extracted. Check chromosome names and FASTA files.")
        return

    # Combine all batches
    print("\nCombining all batches...")
    final_vdf = vaex.concat(all_batches)
    print(f"Total k-mers extracted: {len(final_vdf):,}")

    # Save final file (skip derived features to avoid memory issues)
    print(f"\nSaving final HDF5 file...")
    print(f"Note: Skipping derived features (kmer_rank, relative_position) to avoid memory issues.")
    print(f"Core features included: {list(final_vdf.columns)}")
    final_vdf.export_hdf5(output_file)

    file_size_mb = os.path.getsize(output_file) / (1024**2)
    print(f"\n{'='*60}")
    print(f"SUCCESS!")
    print(f"{'='*60}")
    print(f"Output file: {output_file}")
    print(f"File size: {file_size_mb:.1f} MB")
    print(f"Total rows: {len(final_vdf):,}")
    print(f"Columns: {list(final_vdf.columns)}")
    print(f"{'='*60}")

    # Print sample
    print("\nSample data (first 5 rows):")
    print(final_vdf.head(5).to_pandas_df())


def main():
    parser = argparse.ArgumentParser(
        description="Extract k-mers from chromosomes for PyPair training",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Extract chr21 only
    python extract_chromosome_data.py \\
        --fasta data/genomes/chr21.fa \\
        --output data/reference_genomes/chr21_kmers.hdf5 \\
        --chromosomes chr21

    # Extract chr21 and chr22
    python extract_chromosome_data.py \\
        --fasta data/genomes/chr21.fa data/genomes/chr22.fa \\
        --output data/reference_genomes/chr21_22_kmers.hdf5 \\
        --chromosomes chr21 chr22

    # Extract all chromosomes from file
    python extract_chromosome_data.py \\
        --fasta data/genomes/GRCh38.fa \\
        --output data/reference_genomes/full_genome_kmers.hdf5
        """
    )

    parser.add_argument(
        "--fasta",
        nargs="+",
        required=True,
        help="Input FASTA file(s) containing reference genome"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output HDF5 file path"
    )
    parser.add_argument(
        "--chromosomes",
        nargs="*",
        default=None,
        help="Specific chromosome names to extract (e.g., chr21 chr22). If not specified, extracts all."
    )
    parser.add_argument(
        "--k",
        type=int,
        default=15,
        help="k-mer size (default: 15)"
    )
    parser.add_argument(
        "--no-rc",
        action="store_true",
        help="Exclude reverse complement k-mers"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1_000_000,
        help="Number of k-mers to process per batch (default: 1,000,000)"
    )

    args = parser.parse_args()

    # Validate inputs
    for fasta_file in args.fasta:
        if not os.path.exists(fasta_file):
            print(f"ERROR: FASTA file not found: {fasta_file}")
            sys.exit(1)

    # Extract data
    extract_chromosome_data(
        fasta_files=args.fasta,
        chromosomes=args.chromosomes,
        output_file=args.output,
        k=args.k,
        include_rc=not args.no_rc,
        batch_size=args.batch_size
    )


if __name__ == "__main__":
    main()
