"""
PyPair: Genome Index Generation Module
=======================================

Generates suffix array and FM-index from FASTA reference genome.
Equivalent to BWA's 'bwa index <ref.fa>' command.

Usage:
    from models.index import GenIndex

    indexer = GenIndex()
    index_df = indexer.generate_index(reference_genome, seed_length=28)
"""

import pandas as pd
import os


def rotations(s, index):
    """
    Generate all rotations of input string s.

    Creates the Burrows-Wheeler Matrix by rotating the input string
    by one character at each step.

    Args:
        s: Input string (reference genome with '$' sentinel)
        index: DataFrame to store rotations

    Returns:
        DataFrame with 'pre_sorted' and 'org_location' columns
    """
    print("Making rotations")
    rotations_list = []

    for i in range(len(s)):
        row = s[i:] + s[:i]
        rotations_list.append(row)

    index['pre_sorted'] = rotations_list
    index['org_location'] = index['pre_sorted'].index

    return index


def suffix_array(index):
    """
    Build suffix array from rotations by lexicographic sorting.

    The suffix array is the key data structure for efficient
    string matching in genomic sequences.

    Args:
        index: DataFrame with 'pre_sorted' column

    Returns:
        DataFrame with 'suffix_array' column (sorted suffixes)
    """
    print("Building suffix array column")
    index['suffix_array'] = index['pre_sorted']
    index = index.sort_values(by=['suffix_array'])
    index = index.drop(['pre_sorted'], axis=1)
    index = index.reset_index(drop=True)

    return index


def first_last(index):
    """
    Build First and Last (BWT) columns of the FM-index.

    - First column: First character of each sorted suffix
    - Last column (BWT): Burrows-Wheeler Transform (last character of each rotation)

    Args:
        index: DataFrame with 'suffix_array' column

    Returns:
        DataFrame with 'first' and 'BWT' columns
    """
    print("Building BWT")
    first = []
    last = []

    for row in index['suffix_array']:
        first.append(row[:1])
        last.append(row[-1:])

    index['first'] = first
    index['BWT'] = last

    return index


def char_count(index):
    """
    Populate FM-index character counts from BWT.

    Creates cumulative count columns for each character ($, A, C, G, T).
    These counts enable O(1) character occurrence queries.

    Args:
        index: DataFrame with 'BWT' column

    Returns:
        DataFrame with character count columns
    """
    print("Populating FM index character counts")

    symbols = ['$', 'A', 'C', 'G', 'T']

    for symbol in symbols:
        count = 0
        count_list = []
        for char in index['BWT']:
            if char == symbol:
                count += 1
            count_list.append(count)
        index[symbol] = count_list

    return index


def k_length_seed(k, index):
    """
    Extract k-length seed strings from suffix array.

    Equivalent to BWA's seeding step. Each k-mer represents
    a potential alignment anchor point.

    Args:
        k: Seed length (default 28bp)
        index: DataFrame with 'suffix_array' column

    Returns:
        DataFrame with 'k-seed' column
    """
    print("Generating k-length seeds")
    k_length_seed_string = []

    try:
        for suffix in index['suffix_array']:
            k_seed = suffix[:k]
            k_length_seed_string.append(k_seed)
    except Exception as e:
        print(f"An error occurred generating the k-length seed from suffix array: {e}")

    index['k-seed'] = k_length_seed_string

    return index


def k_length_seed_extension(k, index):
    """
    Extract extended seed region for full read mapping.

    Extracts characters k to 2k for downstream pairwise alignment.

    Args:
        k: Seed length
        index: DataFrame with 'suffix_array' column

    Returns:
        DataFrame with 'k-seed-extend' column
    """
    k_length_seed_extension_list = []

    try:
        for suffix in index['suffix_array']:
            end_point = k + k
            k_extend = suffix[k:end_point]
            k_length_seed_extension_list.append(k_extend)
    except Exception as e:
        print(f"An error occurred creating an extended full read length seed: {e}")

    index['k-seed-extend'] = k_length_seed_extension_list

    return index


class GenIndex:
    """
    Genome index generation class.

    Generates suffix array, BWT, and FM-index from reference genome.
    Should be called first when model weights are not available for
    your genome reference.

    Methods:
        generate_index(reference_genome, seed_length): Build complete index

    Output Columns:
        - org_location: Original genomic position
        - suffix_array: Lexicographically sorted suffix
        - first: First character of suffix
        - BWT: Burrows-Wheeler Transform
        - $, A, C, G, T: FM-index cumulative counts
        - k-seed: First k characters of suffix
        - k-seed-extend: Characters k to 2k
    """

    def generate_index(self, reference_genome, seed_length):
        """
        Generate reference index from genome string.

        Equivalent to 'bwa index <ref.fa>' but outputs a single
        CSV file for training data exploration.

        Args:
            reference_genome: DNA sequence string (e.g., from FASTA)
            seed_length: k-mer length (default 28)

        Returns:
            DataFrame with complete index

        Output:
            data/reference_genomes/indexed_reference.csv
        """
        # Add sentinel character
        seed = "$" + reference_genome

        index = pd.DataFrame()

        # Build index components
        print(f"Indexing {len(reference_genome):,}bp reference genome...")
        index = rotations(seed, index)
        index = suffix_array(index)
        index = first_last(index)
        index = char_count(index)
        index = k_length_seed(seed_length, index)
        index = k_length_seed_extension(seed_length, index)

        # Save to file
        output_dir = "data/reference_genomes"
        try:
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, "indexed_reference.csv")
            index.to_csv(output_file, index=False)
            print(f"[OK] Index saved: {output_file}")
        except OSError as e:
            print(f"Cannot create file or directory: {e}")
            output_dir = "data_sets"
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, "indexed_reference.csv")
            index.to_csv(output_file, index=False)
            print(f"[OK] Index saved (fallback): {output_file}")

        print(f"Index shape: {index.shape}")
        return index
