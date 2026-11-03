"""
PyPair: SAM/BAM File Generation Module
=======================================

Generates Sequence Alignment Map (SAM/BAM) files from alignment results.
SAM is the standard format for storing genomic alignments.

Usage:
    from models.sam import create_bam_file

    alignments = [
        {
            "query_name": "read1",
            "sequence": "ACGTACGT",
            "start": 12345,
            "cigar": "8M"
        }
    ]
    create_bam_file(alignments, "output.bam")
"""

import pysam


def create_bam_file(alignments, output_filename, reference_name="ref", reference_length=1000):
    """
    Create a BAM file from sequence alignments.

    BAM (Binary Alignment Map) is the compressed binary version of SAM.
    Used by downstream tools like SAMtools, IGV, GATK, etc.

    Args:
        alignments: List of alignment dictionaries with keys:
            - 'query_name': Read identifier
            - 'sequence': Read sequence
            - 'start': 0-based genomic position (predicted by ML model)
            - 'cigar': CIGAR string (from Smith-Waterman)
        output_filename: Path to output BAM file
        reference_name: Reference sequence name (default: "ref")
        reference_length: Reference sequence length (default: 1000)

    Output:
        Binary BAM file at output_filename

    SAM/BAM Fields (11 required):
        1. QNAME: Query name
        2. FLAG: Bitwise flags (0=unmapped, 16=reverse, etc.)
        3. RNAME: Reference sequence name
        4. POS: 1-based leftmost position
        5. MAPQ: Mapping quality (255=unavailable)
        6. CIGAR: Alignment operations
        7. RNEXT: Reference name of mate (-1=unavailable)
        8. PNEXT: Position of mate (-1=unavailable)
        9. TLEN: Template length (0=unavailable)
        10. SEQ: Query sequence
        11. QUAL: Base quality scores

    Example:
        >>> alignments = [
        ...     {
        ...         "query_name": "read_001",
        ...         "sequence": "ACGTACGTACGT",
        ...         "start": 1000,
        ...         "cigar": "12M"
        ...     }
        ... ]
        >>> create_bam_file(alignments, "output.bam")
    """
    # Create SAM header
    header = {
        "HD": {"VN": "1.0"},  # Header version
        "SQ": [{
            "SN": reference_name,
            "LN": reference_length
        }]
    }

    with pysam.AlignmentFile(output_filename, "wb", header=header) as bamfile:
        for aln in alignments:
            # Create AlignedSegment object
            a = pysam.AlignedSegment()

            # Required fields
            a.query_name = aln["query_name"]
            a.query_sequence = aln["sequence"]
            a.flag = 0                          # 0 = mapped, primary alignment
            a.reference_id = 0                  # First (and only) reference
            a.reference_start = aln["start"]    # 0-based position
            a.mapping_quality = 255             # 255 = unavailable
            a.cigarstring = aln["cigar"]

            # Mate pair fields (not used in single-end alignment)
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0

            # Quality scores (assume high quality Q30 for all bases)
            a.query_qualities = pysam.array.array('B', [30] * len(aln["sequence"]))

            # Write to BAM file
            bamfile.write(a)

    print(f"[OK] BAM file created: {output_filename}")
    print(f"     Alignments written: {len(alignments)}")


def create_sam_file(alignments, output_filename, reference_name="ref", reference_length=1000):
    """
    Create a SAM (text) file from sequence alignments.

    SAM is human-readable but larger than BAM. Use for debugging or
    when text format is required.

    Args:
        alignments: List of alignment dictionaries (same format as create_bam_file)
        output_filename: Path to output SAM file
        reference_name: Reference sequence name
        reference_length: Reference sequence length

    Output:
        Text SAM file at output_filename

    Example:
        >>> alignments = [...]
        >>> create_sam_file(alignments, "output.sam")
    """
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{
            "SN": reference_name,
            "LN": reference_length
        }]
    }

    with pysam.AlignmentFile(output_filename, "w", header=header) as samfile:
        for aln in alignments:
            a = pysam.AlignedSegment()
            a.query_name = aln["query_name"]
            a.query_sequence = aln["sequence"]
            a.flag = 0
            a.reference_id = 0
            a.reference_start = aln["start"]
            a.mapping_quality = 255
            a.cigarstring = aln["cigar"]
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            a.query_qualities = pysam.array.array('B', [30] * len(aln["sequence"]))

            samfile.write(a)

    print(f"[OK] SAM file created: {output_filename}")
    print(f"     Alignments written: {len(alignments)}")


def alignment_to_dict(query_name, sequence, predicted_position, cigar_string):
    """
    Convert alignment components to dictionary format for SAM/BAM generation.

    Args:
        query_name: Read identifier
        sequence: DNA sequence
        predicted_position: Genomic position (from ML model)
        cigar_string: CIGAR string (from Smith-Waterman)

    Returns:
        dict: Alignment dictionary ready for create_bam_file()

    Example:
        >>> aln_dict = alignment_to_dict("read1", "ACGT", 1000, "4M")
        >>> create_bam_file([aln_dict], "output.bam")
    """
    return {
        "query_name": query_name,
        "sequence": sequence,
        "start": predicted_position,
        "cigar": cigar_string
    }


def batch_alignments_to_bam(predictions, sequences, query_names, cigar_strings, output_file):
    """
    Batch process multiple predictions into BAM file.

    Convenience function for converting PyPair ML predictions + Smith-Waterman
    alignments into standard BAM format.

    Args:
        predictions: List of predicted genomic positions (from ML model)
        sequences: List of read sequences
        query_names: List of read identifiers
        cigar_strings: List of CIGAR strings (from Smith-Waterman)
        output_file: Path to output BAM file

    Example:
        >>> predictions = [1000, 2000, 3000]
        >>> sequences = ["ACGT", "TGCA", "GGCC"]
        >>> names = ["read1", "read2", "read3"]
        >>> cigars = ["4M", "4M", "4M"]
        >>> batch_alignments_to_bam(predictions, sequences, names, cigars, "out.bam")
    """
    alignments = []

    for pred, seq, name, cigar in zip(predictions, sequences, query_names, cigar_strings):
        alignments.append({
            "query_name": name,
            "sequence": seq,
            "start": pred,
            "cigar": cigar
        })

    create_bam_file(alignments, output_file)

    return alignments
