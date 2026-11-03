"""
PyPair: Pairwise Sequence Alignment Module
===========================================

Smith-Waterman local alignment algorithm for downstream confirmation
and final sequence alignment mapping.

Usage:
    from models.alignment import SmithWaterman

    sw = SmithWaterman()
    aligned_seq1, aligned_seq2, score = sw.compute("ACGTACGT", "ACGTACGT")
    cigar = sw.compute_cigar(aligned_seq1, aligned_seq2)
"""


class SmithWaterman:
    """
    Smith-Waterman local pairwise sequence alignment algorithm.

    Implements affine gap penalties for biologically realistic alignments.
    Used for downstream confirmation after ML-based seed mapping.

    Scoring:
        - Match: +2
        - Mismatch: -1
        - Gap open: -2
        - Gap extend: -1
    """

    @staticmethod
    def compute(
        seq1,
        seq2,
        match_score=2,
        mismatch_score=-1,
        gap_start=-2,
        gap_extend=-1
    ):
        """
        Compute Smith-Waterman local alignment between two sequences.

        Uses three matrices:
        - H: Main scoring matrix
        - E: Horizontal gap (insertion in seq2)
        - F: Vertical gap (deletion in seq2)

        Args:
            seq1: First DNA sequence
            seq2: Second DNA sequence
            match_score: Score for matching bases (default: +2)
            mismatch_score: Penalty for mismatches (default: -1)
            gap_start: Penalty for opening a gap (default: -2)
            gap_extend: Penalty for extending a gap (default: -1)

        Returns:
            tuple: (aligned_seq1, aligned_seq2, max_score)
                - aligned_seq1: Aligned version of seq1 (with '-' for gaps)
                - aligned_seq2: Aligned version of seq2 (with '-' for gaps)
                - max_score: Alignment score

        Example:
            >>> sw = SmithWaterman()
            >>> s1, s2, score = sw.compute("ACGT", "AGGT")
            >>> print(f"Score: {score}")
            >>> print(f"Seq1: {s1}")
            >>> print(f"Seq2: {s2}")
        """
        m, n = len(seq1), len(seq2)

        # Initialize scoring matrices
        H = [[0] * (n + 1) for _ in range(m + 1)]
        E = [[0] * (n + 1) for _ in range(m + 1)]
        F = [[0] * (n + 1) for _ in range(m + 1)]

        def score(x, y):
            """Calculate match/mismatch score."""
            return match_score if seq1[x - 1] == seq2[y - 1] else mismatch_score

        # Fill matrices using dynamic programming
        max_score = 0
        max_pos = (0, 0)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Horizontal gap (insertion in seq2)
                E[i][j] = max(
                    H[i][j - 1] + gap_start,
                    E[i][j - 1] + gap_extend
                )

                # Vertical gap (deletion in seq2)
                F[i][j] = max(
                    H[i - 1][j] + gap_start,
                    F[i - 1][j] + gap_extend
                )

                # Main scoring matrix
                H[i][j] = max(
                    0,                              # Local alignment (can restart)
                    H[i - 1][j - 1] + score(i, j),  # Match/mismatch
                    E[i][j],                        # Insertion
                    F[i][j]                         # Deletion
                )

                # Track maximum score position
                if H[i][j] > max_score:
                    max_score = H[i][j]
                    max_pos = (i, j)

        # Backtrack to recover alignment
        aligned_seq1, aligned_seq2 = [], []
        i, j = max_pos

        while H[i][j] > 0:
            if H[i][j] == H[i - 1][j - 1] + score(i, j):
                # Match/mismatch
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif H[i][j] == E[i][j]:
                # Insertion (gap in seq1)
                aligned_seq1.append("-")
                aligned_seq2.append(seq2[j - 1])
                j -= 1
            elif H[i][j] == F[i][j]:
                # Deletion (gap in seq2)
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append("-")
                i -= 1

        # Reverse alignments (backtracking goes backwards)
        return (
            "".join(reversed(aligned_seq1)),
            "".join(reversed(aligned_seq2)),
            max_score
        )

    @staticmethod
    def compute_cigar(aligned_seq1, aligned_seq2):
        """
        Convert aligned sequences to CIGAR string.

        CIGAR (Compact Idiosyncratic Gapped Alignment Report) is the
        standard format for representing alignments in SAM/BAM files.

        Operations:
            M: Match or mismatch
            I: Insertion (gap in reference)
            D: Deletion (gap in query)

        Args:
            aligned_seq1: Aligned sequence 1 (with '-' for gaps)
            aligned_seq2: Aligned sequence 2 (with '-' for gaps)

        Returns:
            str: CIGAR string

        Example:
            >>> sw = SmithWaterman()
            >>> cigar = sw.compute_cigar("ACG--T", "ACGTTT")
            >>> print(cigar)  # "3M2I1M"
        """
        cigar = []
        count = 0
        last_op = ""

        for a, b in zip(aligned_seq1, aligned_seq2):
            # Determine operation
            if a == "-" or b == "-":
                op = "I" if a == "-" else "D"
            else:
                op = "M"  # Match or mismatch (both represented as M)

            # Run-length encode operations
            if op == last_op:
                count += 1
            else:
                if last_op:
                    cigar.append(f"{count}{last_op}")
                count = 1
                last_op = op

        # Add final operation
        if last_op:
            cigar.append(f"{count}{last_op}")

        return "".join(cigar)


# Convenience function for quick alignment
def align_sequences(seq1, seq2):
    """
    Quick alignment of two sequences using default parameters.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        dict: Alignment results with keys:
            - 'aligned_seq1': Aligned first sequence
            - 'aligned_seq2': Aligned second sequence
            - 'score': Alignment score
            - 'cigar': CIGAR string

    Example:
        >>> from models.alignment import align_sequences
        >>> result = align_sequences("ACGTACGT", "ACGTTACGT")
        >>> print(result['cigar'])
    """
    sw = SmithWaterman()
    aligned_seq1, aligned_seq2, score = sw.compute(seq1, seq2)
    cigar = sw.compute_cigar(aligned_seq1, aligned_seq2)

    return {
        'aligned_seq1': aligned_seq1,
        'aligned_seq2': aligned_seq2,
        'score': score,
        'cigar': cigar
    }
