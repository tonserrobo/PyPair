"""
PyPair: Input/Output Processing Module
=======================================

Handles FASTA and FASTQ file parsing for reference genomes and reads.

Usage:
    from models.io import IO_processing

    io_proc = IO_processing()
    seq_id, seq_str, seq_len = io_proc.parse_reference("reference.fasta", 28)
"""

from Bio import SeqIO
import os


class IO_processing:
    """
    Main IO and data preprocessing class for parsing genomic data.

    Handles:
    - FASTA reference genome parsing
    - FASTQ read file parsing
    - Model weight persistence
    """

    def parse_reference(self, input_file, batch_size):
        """
        Parse FASTA (.fa) file to extract genomic sequence data.

        Args:
            input_file: Path to FASTA file
            batch_size: For future batching functionality (currently unused)

        Returns:
            tuple: (sequence_id, sequence_string, sequence_length)

        Example:
            >>> io_proc = IO_processing()
            >>> seq_id, seq, length = io_proc.parse_reference("chr21.fasta", 28)
            >>> print(f"Loaded {seq_id}: {length:,}bp")
        """
        print("... Parsing input reference genome")

        sequence_id = None
        sequence_string = None
        sequence_length = None

        try:
            for seq_record in SeqIO.parse(input_file, "fasta"):
                print(f"\nParsing reference sequence {input_file}")
                # Note: Assumes single seq_record per file
                # If multiple records exist, last one is retained
                sequence_id = seq_record.id
                sequence_string = str(seq_record.seq).upper()
                sequence_length = len(seq_record)

                print(f"Sequence ID: {sequence_id}")
                print(f"Sequence Length: {sequence_length:,}bp")

        except Exception as e:
            print(f"Error parsing input sequence file: {e}")
            print("Check that input is in FASTA format")

        return sequence_id, sequence_string, sequence_length

    # Legacy method name for backwards compatibility
    def pharse_reference(self, input_file, batch_size):
        """Legacy method name. Use parse_reference() instead."""
        return self.parse_reference(input_file, batch_size)

    def save_model_weights(self, model, input_file):
        """
        Save model training weights as .pkl file.

        Args:
            model: Trained scikit-learn model
            input_file: Base filename (without extension)

        Returns:
            int: 0 on success

        Output:
            pyseed/model_weights/<input_file>.pkl
        """
        import pickle

        model_name = input_file + ".pkl"
        file_path = "pyseed/model_weights/"
        pkl_filename = file_path + model_name

        os.makedirs(file_path, exist_ok=True)

        with open(pkl_filename, 'wb') as file:
            pickle.dump(model, file)

        print(f"[OK] Model weights saved: {pkl_filename}")
        return 0

    def weights_lookup(self, model_training_weights, training_flag):
        """
        Check for pre-trained model weights and trigger retraining if not found.

        Args:
            model_training_weights: Base filename of weights to look up
            training_flag: Current training flag status

        Returns:
            bool: True if training needed, False if weights found
        """
        print("... Looking up available model weights")

        weights_path = 'model_weights/' + model_training_weights + '.pkl'

        if os.path.isfile(weights_path):
            training_flag = False
            print(f"[OK] Found pre-trained weights: {weights_path}")
        else:
            training_flag = True
            print("Training weights not found, activating model training")

        return training_flag


def parse_fastq(input_file):
    """
    Parse FASTQ (.fastq) file to extract read sequences.

    Args:
        input_file: Path to FASTQ file

    Returns:
        list: List of DNA sequence strings (uppercase)

    Example:
        >>> from models.io import parse_fastq
        >>> reads = parse_fastq("reads.fastq")
        >>> print(f"Loaded {len(reads)} reads")
    """
    print("... Parsing input FASTQ file")
    sequences = []

    try:
        for seq_record in SeqIO.parse(input_file, "fastq"):
            sequence_string = str(seq_record.seq).upper()
            sequences.append(sequence_string)
            print(f"Parsing sequence {seq_record.id} with length {len(seq_record)}")
    except Exception as e:
        print(f"Error parsing input sequence file: {e}")

    print(f"[OK] Loaded {len(sequences)} sequences")
    return sequences


def parse_fasta(input_file):
    """
    Parse FASTA (.fasta/.fa) file to extract sequences.

    Args:
        input_file: Path to FASTA file

    Returns:
        list: List of tuples (seq_id, sequence_string)

    Example:
        >>> from models.io import parse_fasta
        >>> seqs = parse_fasta("reference.fasta")
        >>> for seq_id, seq in seqs:
        ...     print(f"{seq_id}: {len(seq)}bp")
    """
    print("... Parsing input FASTA file")
    sequences = []

    try:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            seq_id = seq_record.id
            sequence_string = str(seq_record.seq).upper()
            sequences.append((seq_id, sequence_string))
            print(f"Parsing sequence {seq_id} with length {len(seq_record)}bp")
    except Exception as e:
        print(f"Error parsing input sequence file: {e}")

    print(f"[OK] Loaded {len(sequences)} sequences")
    return sequences
