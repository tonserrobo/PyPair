import pandas as pd
import os


def rotations(s, index):
    """Rotate index string by 1"""
    print("Making rotations")
    rotations = []
    for i in range(len(s)):
        end_number = len(s)
        status = round(((i / end_number) * 100), 2)
        row = s[i:] + s[:i]
        rotations.append(row)
    index['pre_sorted'] = rotations
    del rotations
    index['org_location'] = index['pre_sorted'].index
    return index

def suffix_array(index):
    """ Method for building the suffix array from input reference """
    print("Building suffix array column")
    index['suffix_array'] = index['pre_sorted']
    index = index.sort_values(by=['suffix_array'])
    index = index.drop(['pre_sorted'], axis=1)
    index = index.reset_index(drop=True)
    return index

def first_last(index):
    """ Method for building the first and last(BWT) column of the index file """
    print("Building BWT")
    first = []
    last = []

    for row in index['suffix_array']:
        f = row[:1]
        first.append(f)
    index['first'] = first
    del first

    for row in index['suffix_array']:
        l = row[-1:]
        last.append(l)
    index['BWT'] = last
    del last
    return index

def char_count(index):
    """ Method for populating FM index char counts from BWT """
    print("Populating FM index character counts")
    index['$'] = ''
    index['A'] = ''
    index['C'] = ''
    index['G'] = ''
    index['T'] = ''

    symbols = ['$', 'A', 'C', 'G', 'T']
    for symbol in symbols:
        count = 0
        count_list = []
        for i in index['BWT']:
            if i == symbol:
                count += 1
            count_list.append(count)
        index[symbol] = count_list
        del count_list
    return index

def k_length_seed(k, index):
    """ Method to return k-length seed string without IP-BWT concatenation """
    print("Generating k-length seeds")
    k_length_seed_string = []
    try:
        for i in index['suffix_array']:
            k_seed = i[:k]
            k_length_seed_string.append(k_seed)
    except Exception as e:
        print(f"An error occurred generating the k-length seed from suffix array: {e}")

    index['k-seed'] = k_length_seed_string
    del k_length_seed_string
    return index

def k_length_seed_extension(k, index):
    k_length_seed_extension = []
    try:
        for i in index['suffix_array']:
            end_point = k + k
            k_extend = i[k:end_point]
            k_length_seed_extension.append(k_extend)
    except Exception as e:
        print(f"An error occurred creating an extended full read length seed: {e}")

    index['k-seed-extend'] = k_length_seed_extension
    del k_length_seed_extension
    return index

class GenIndex():
    """
    Genome index class. Should model weights not be available for your genome reference, GenIndex and associated methods are called first. Contains the following methods:
        1. rotations
        2. suffix_array
        3. first_last
        4. char_count
        5. ip_matrix
    """

    def generate_index(self, reference_genome, seed_length):
        """
        Class GenIndex() method call to read reference file and generate reference index.
        """
        seed = "$" + reference_genome  # Modify so will work with passed string sequence

        index = pd.DataFrame()

        # Pass seed through rotations, suffix generation, BWT generation and FM index generation and k-seed generation
        index = rotations(seed, index)
        index = suffix_array(index)
        index = first_last(index)
        index = char_count(index)
        index = k_length_seed(seed_length, index)

        # Full read mapping extension
        index = k_length_seed_extension(seed_length, index)

        output_dir = "data/reference_genomes"
        try:
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, "indexed_reference.csv")
            index.to_csv(output_file, index=False)
        except OSError as e:
            print(f"Cannot create file or directory: {e}")
            output_dir = "data_sets"
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, "indexed_reference.csv")
            index.to_csv(output_file, index=False)

        return index
