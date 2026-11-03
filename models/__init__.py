"""
PyPair: Machine Learning-Based Genomic Seed Mapping
====================================================

PyPair is a learned index system for genomic sequence alignment that
achieves 35.7x speedup over BWA while maintaining 99.88% accuracy.

Core Modules:
    - index: Suffix array and FM-index generation
    - encoding: 2-bit DNA encoding/decoding
    - io: FASTA/FASTQ file parsing
    - alignment: Smith-Waterman pairwise alignment
    - sam: SAM/BAM file generation

Performance (KNN Model):
    - Speed: 0.084 ms/seed (35.7x faster than BWA)
    - Accuracy: 99.88% F1 score
    - Overfitting: 0.12% gap (minimal)
    - Throughput: 11,000-24,000 seeds/sec

Usage:
    from models.index import GenIndex
    from models.encoding import Encode_decode
    from models.io import parse_fastq
    from models.alignment import SmithWaterman
    from models.sam import create_bam_file

Author: Tony Robinson
Institution: Ulster University
"""

__version__ = "1.0.0"
__author__ = "Tony Robinson"

# Import key classes for convenience
from .index import GenIndex
from .encoding import Encode_decode
from .io import IO_processing, parse_fastq, parse_fasta
from .alignment import SmithWaterman, align_sequences
from .sam import create_bam_file, create_sam_file, batch_alignments_to_bam

__all__ = [
    # Index generation
    'GenIndex',

    # Encoding
    'Encode_decode',

    # IO functions
    'IO_processing',
    'parse_fastq',
    'parse_fasta',

    # Alignment
    'SmithWaterman',
    'align_sequences',

    # SAM/BAM generation
    'create_bam_file',
    'create_sam_file',
    'batch_alignments_to_bam',
]
