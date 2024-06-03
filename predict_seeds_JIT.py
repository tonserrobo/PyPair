import argparse
from Bio import SeqIO
import pickle
import pandas as pd
from encode_decode import Encode_decode
import pysam  # need to import with bioconda
from numba import jit
import numpy as np


# Smith-Waterman Algorithm
class SmithWaterman:
    @staticmethod
    @jit(nopython=True)
    def compute(
        seq1, seq2, match_score=2, mismatch_score=-1, gap_start=-2, gap_extend=-1
    ):
        m, n = len(seq1), len(seq2)
        H = np.zeros((m + 1, n + 1))
        E = np.zeros((m + 1, n + 1))
        F = np.zeros((m + 1, n + 1))

        max_score = 0
        max_pos = (0, 0)
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                E[i][j] = max(H[i][j - 1] + gap_start, E[i][j - 1] + gap_extend)
                F[i][j] = max(H[i - 1][j] + gap_start, F[i - 1][j] + gap_extend)
                H[i][j] = max(0, H[i - 1][j - 1] + score(i, j), E[i][j], F[i][j])

                if H[i][j] > max_score:
                    max_score = H[i][j]
                    max_pos = (i, j)

        aligned_seq1, aligned_seq2 = [], []
        i, j = max_pos
        while H[i][j] > 0:
            if H[i][j] == H[i - 1][j - 1] + score(i, j):
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif H[i][j] == E[i][j]:
                aligned_seq1.append("-")
                aligned_seq2.append(seq2[j - 1])
                j -= 1
            elif H[i][j] == F[i][j]:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append("-")
                i -= 1

        return (
            "".join(reversed(aligned_seq1)),
            "".join(reversed(aligned_seq2)),
            max_score,
        )

    @staticmethod
    def score(x, y):
        return match_score if seq1[x - 1] == seq2[y - 1] else mismatch_score

    @staticmethod
    def compute_cigar(aligned_seq1, aligned_seq2):
        cigar = []
        count = 0
        last_op = ""

        for a, b in zip(aligned_seq1, aligned_seq2):
            op = ""
            if a == "-" or b == "-":
                op = "I" if a == "-" else "D"
            else:
                op = "M"

            if op == last_op:
                count += 1
            else:
                if last_op:
                    cigar.append(f"{count}{last_op}")
                count = 1
                last_op = op

        if last_op:
            cigar.append(f"{count}{last_op}")

        return "".join(cigar)


# Load Model
def load_model(model_filename):
    try:
        with open(model_filename, "rb") as file:
            return pickle.load(file)
    except FileNotFoundError:
        print(f"Model file {model_filename} not found.")
        return None


# Get Reference at Predicted Location
def get_reference_at_location(ref_df, location, window_size=28):
    try:
        reference = ref_df.loc[ref_df["position"] == location, "sequence"].values[0]
        prev_str = (
            ref_df.loc[ref_df["position"] == location - window_size, "sequence"].values[
                0
            ]
            if location - window_size >= 0
            else ""
        )
        next_str = (
            ref_df.loc[ref_df["position"] == location + window_size, "sequence"].values[
                0
            ]
            if location + window_size < len(ref_df)
            else ""
        )
        return prev_str, reference, next_str
    except IndexError:
        return "", "", ""


def create_bam_file(alignments, output_filename):
    """
    Create a BAM file from sequence alignments.

    :param alignments: A list of alignment data.
    :param output_filename: Name of the output BAM file.
    """

    # Create a SAM file header
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "ref"}]}  # Example header

    with pysam.AlignmentFile(output_filename, "wb", header=header) as bamfile:
        for aln in alignments:
            # Create an AlignmentSegment object
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
            a.query_qualities = pysam.array.array(
                "B", [30] * len(aln["sequence"])
            )  # Assuming high quality
            bamfile.write(a)


def preprocess_input(sequences, encoder_decoder):
    processed_reads = []
    try:
        for seq in sequences:
            encoded_seq = encoder_decoder.encode(seq)
            processed_reads.append(encoded_seq)
    except Exception as e:
        print(f"Error processing input: {e}")
    return processed_reads


def parse_fastq(input_file):
    """Parse .fastq file to extract genomic string data"""
    print("... Parsing input FASTQ file")
    sequences = []
    try:
        for seq_record in SeqIO.parse(input_file, "fastq"):
            sequence_string = str(seq_record.seq).upper()
            sequences.append(sequence_string)
            print(f"Parsing sequence {seq_record.id} with length {len(seq_record)}")
    except Exception as e:
        print(f"Error parsing input sequence file: {e}")
    return sequences


def main(input_file):
    try:
        encoder_decoder = Encode_decode()
        model = load_model("centroid_model.joblib")
        if not model:
            return

        sequences = parse_fastq(input_file)
        input_reads = preprocess_input(sequences, encoder_decoder)
        if not input_reads:
            return

        reference = pd.read_csv("reference_samples/S1_reference.csv")
        predictions = model.predict(input_reads)

        results = []
        for i, pred in enumerate(predictions):
            seq = sequences[i]  # Use sequence from FASTQ file
            prev_str, ref_str, next_str = get_reference_at_location(reference, pred)
            aligned_seq1, aligned_seq2, alignment_score = SmithWaterman.compute(
                seq, ref_str
            )
            cigar = SmithWaterman.compute_cigar(aligned_seq1, aligned_seq2)

            results.append(
                {
                    "Read": seq,
                    "Predicted Location": pred,
                    "Reference": ref_str,
                    "Alignment Score": alignment_score,
                    "Aligned Sequence 1": aligned_seq1,
                    "Aligned Sequence 2": aligned_seq2,
                    "CIGAR": cigar,
                }
            )

        results_df = pd.DataFrame(results)
        results_df.to_csv("Validation/alignment_results.csv", index=False)
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process a FASTQ file for seed prediction."
    )
    parser.add_argument("input_file", type=str, help="Path to the FASTQ file")
    args = parser.parse_args()

    main(args.input_file)
