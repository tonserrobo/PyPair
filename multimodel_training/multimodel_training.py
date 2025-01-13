import os
import sys
import pandas as pd ## this needs to be a polars lib
from Bio.Seq import Seq
from collections import OrderedDict

# class imports
from io_processing import IO_processing
from gen_index import GenIndex
from input_encoder import InputEncoder
from training_data_preparation import TrainingDataGenerator
from models.benchmark import SeedAlignerBenchmark


if __name__ == '__main__':
    """ Main method for multimodal training"""

    # data preparation
    config = {
        "seed_length": 28,
        "read_length": 100,
        "training_iterations": 8,
    }

    try:
        data_preprocessor = InputEncoder(config['seed_length'], config['read_length'], config['training_iterations'])
    except Exception as e:
        print("Error preparing data for training data generation: ", e)
    try:
        training_data_generator = TrainingDataGenerator()  # Gen training data based upon preprocessed input.
    except Exception as e:
        print("Error generating training data generator: ", e)

    """ Initial test - Binary hash search """
    print("beginning initial benchmark")
    benchmark = SeedAlignerBenchmark()  # error here.

