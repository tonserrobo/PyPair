import os
import polars as pl
from Bio.Seq import Seq
from collections import OrderedDict
from io_processing import IO_processing
from gen_index import GenIndex


def reverse_complement(sequence_df):
    """ Compute reverse complement of k-seed strings using Polars """
    sequence_df = sequence_df.with_columns(
        pl.col("k-seed").apply(lambda x: str(Seq(x).reverse_complement()), return_dtype=pl.Utf8).alias("rc_seeds")
    )
    return sequence_df


def encode(sequence):
    """Encode a string using 2-bit encoding"""
    dictionary = {'$': '', ',': '', 'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    if isinstance(sequence, str):
        return sequence.translate(str.maketrans(dictionary))
    return sequence


def encode_dataframe(df):
    """ Encode DataFrame columns using Polars """
    encoded_df = df.select(
        [pl.col(col).apply(lambda x: encode(x) if isinstance(x, str) else x).alias(col) for col in df.columns]
    )
    return encoded_df


def run_length_encoding(input_string):
    """ Run length encoding of a string """
    count_dict = OrderedDict.fromkeys(input_string, 0)
    for ch in input_string:
        count_dict[ch] += 1

    output = ''.join(f"{key}{val}" for key, val in count_dict.items())
    characters = ''.join(count_dict.keys())
    values = ''.join(map(str, count_dict.values()))
    return output, characters, values


def run_length_encoding_df(df, column_name):
    """ Apply run-length encoding to a DataFrame column """
    outputs, characters, values = zip(
        *df[column_name].apply(lambda x: run_length_encoding(x) if isinstance(x, str) else ("", "", ""))
    )
    return pl.DataFrame({
        f"{column_name}_output": outputs,
        f"{column_name}_character": characters,
        f"{column_name}_value": values,
    })


def extend_dataset(df, config):
    """ Extend training set with repeated instances of class examples """
    extended_rows = []
    unique_classes = df["org_location"].unique().to_list()
    for cls in unique_classes:
        cls_df = df.filter(pl.col("org_location") == cls)
        repeat_times = max(1, config["training_iterations"] - cls_df.shape[0])
        extended_rows.append(cls_df.vstack([cls_df] * (repeat_times - 1)))

    return pl.concat(extended_rows)


def process_fasta_file(file_path, seed_length):
    """ Parse fasta file """
    io_processor = IO_processing()
    return io_processor.pharse_reference(file_path, seed_length)


def index_sequence(sequence_string, seed_length):
    """ Index sequence using GenIndex() """
    index_obj = GenIndex()
    return index_obj.generate_index(sequence_string, seed_length)


class InputEncoder:
    """ Class to handle input processes using the helper functions defined above """

    def __init__(self, config_dict):
        self.seed_length = config_dict['seed_length']
        self.read_length = config_dict['read_length']
        self.training_iterations = config_dict['training_iterations']
        self.folder_path = 'data/reference_samples'

    def process_files(self):
        fasta_files = [f for f in os.listdir(self.folder_path) if f.endswith('.fasta')]

        for i, fasta_file in enumerate(fasta_files, start=1):
            file_path = os.path.join(self.folder_path, fasta_file)

            sequence_id, sequence_string, sequence_length = process_fasta_file(file_path, self.seed_length)
            index = index_sequence(sequence_string, self.seed_length)
            index_df = pl.DataFrame(index)

            # Compute reverse complements
            reverse_complement_df = reverse_complement(index_df)

            # Run length encoding for k-seeds and reverse complements
            kseed_rle = run_length_encoding_df(reverse_complement_df, "k-seed")
            rcseed_rle = run_length_encoding_df(reverse_complement_df, "rc_seeds")
            reverse_complement_df = reverse_complement_df.hstack([kseed_rle, rcseed_rle])

            # Encode dataset
            encoded_index = encode_dataframe(reverse_complement_df)

            # Filter columns based on training data analysis PCA
            training_data_cols = [
                "k-seed", "rc_seeds", "runlength_output", "runlength_rc_output",
                "runlength_rc_character", "runlength_rc_value", "org_location"
            ]
            filtered_data = encoded_index.select(training_data_cols)

            # Extend the dataset
            config = {"training_iterations": self.training_iterations}
            extended_encoded_index = extend_dataset(filtered_data, config)

            # Save the extended encoded index
            output_file_name = os.path.join(self.folder_path, f'S{i}_extended_encoded_reference.csv')
            extended_encoded_index.write_csv(output_file_name)
