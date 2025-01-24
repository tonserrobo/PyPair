import os
import sys
import pandas as pd
from Bio.Seq import Seq
from collections import OrderedDict
from sklearn.model_selection import train_test_split
from io_processing import IO_processing
from gen_index import GenIndex

folder_path = 'data/reference_samples/'

# seed config
config = {
    "seed_length": 28,
    "read_length": 100,
    "training_iterations": 5,
}

def reverse_complement(sequence_set):
    """ Compute reverse compliment of k-seed string """
    sequence_strings = sequence_set['k-seed']

    rc_strings = []
    clean_rc_strings = []
    for string in sequence_strings:
        s_seq = Seq(string)
        rc_string = s_seq.reverse_complement()
        rc_strings.append(rc_string)

    for row in rc_strings:
        new_string = ''.join(row)
        clean_rc_strings.append(new_string)

    sequence_set['rc_seeds'] = clean_rc_strings
    return sequence_set

def encode(_index):
    """Encode the data set using 2-bit encoding"""
    encoded = []
    dictionary = {'$': '', ',':'', 'A': '00', 'C': '01', 'G': '10', 'T': '11'}

    if isinstance(_index, str):
        transTable = _index.maketrans(dictionary)
        txt = _index.translate(transTable)
        encoded = txt
    else:
        for row in _index:
            encoded_row = row.translate(row.maketrans(dictionary))
            encoded.append(encoded_row)
    return encoded

def extend_dataset(_df, _config):
    """ Extend training set with repeated instances of class example """
    classes = _df['org_location'].unique()
    extended_df = pd.DataFrame()

    for cls in classes:
        cls_df = _df[_df['org_location'] == cls]
        repeat_times = max(1, _config['training_iterations'] - cls_df.shape[0])

        extended_cls_df = pd.concat([cls_df] * repeat_times, ignore_index=True)
        extended_df = pd.concat([extended_df, extended_cls_df], ignore_index=True)

    return extended_df

def process_fasta_file(_file_path, _config):
    """ Parse fasta file """
    io_processor = IO_processing()
    _sequence_id, _sequence_string, _sequence_length = io_processor.pharse_reference(_file_path, _config['seed_length'])
    return _sequence_id, _sequence_string, _sequence_length

def index_sequence(_sequence_string, seed_length):
    """ Index sequence using GenIndex() """
    index_obj = GenIndex()
    _index = index_obj.generate_index(_sequence_string, seed_length)
    return _index

def save_index_as_csv(_index, file_name):
    """ Save index function """
    _index_df = pd.DataFrame(_index)
    if os.path.exists(file_name):
        os.remove(file_name)
    _index_df.to_csv(file_name, index=False)

def encode_dataframe(_df):
    """ Encode dataframe """
    encoded_df = _df.copy()
    for column in encoded_df.columns:
        if column != 'org_location':
            encoded_df[column] = encoded_df[column].apply(lambda _x: encode(_x) if isinstance(_x, str) else _x)
    return encoded_df

def runLengthEncoding(_input):
    """ Run length encoding of k-seed string """
    _dict = OrderedDict.fromkeys(_input, 0)

    for ch in _input:
        _dict[ch] += 1
    _output = ''
    characters = ''
    _values = ''
    for key, _value in _dict.items():
        _output = _output + key + str(_value)
        characters = characters + key
        _values = _values + str(_value)
    return _output, characters, _values

fasta_files = [f for f in os.listdir(folder_path) if f.endswith('.fasta')]

for i, fasta_file in enumerate(fasta_files, start=1):
    file_path = os.path.join(folder_path, fasta_file)

    sequence_id, sequence_string, sequence_length = process_fasta_file(file_path, config)
    index = index_sequence(sequence_string, config['seed_length'])
    index_df = pd.DataFrame(index)

    reverse_complement_df = reverse_complement(index_df)

    # Run length encoding
    output = []
    character = []
    value = []

    for seed_string in reverse_complement_df['k-seed']:
        out, char, val = runLengthEncoding(seed_string)
        output.append(out)
        character.append(char)
        value.append(val)

    reverse_complement_df['runlength_output'] = output
    reverse_complement_df['runlength_character'] = character
    reverse_complement_df['runlength_value'] = value

    del output
    del character
    del value

    rc_output = []
    rc_character = []
    rc_value = []

    for rc_seed_string in reverse_complement_df['rc_seeds']:
        rc_out, rc_char, rc_val = runLengthEncoding(rc_seed_string)
        rc_output.append(rc_out)
        rc_character.append(rc_char)
        rc_value.append(rc_val)

    reverse_complement_df['runlength_rc_output'] = rc_output
    reverse_complement_df['runlength_rc_character'] = rc_character
    reverse_complement_df['runlength_rc_value'] = rc_value

    del rc_output
    del rc_character
    del rc_value

    # Encode dataset
    encoded_index = encode_dataframe(reverse_complement_df)

    # Filter columns based on training data analysis PCA
    df_training_data = pd.DataFrame({
        "first": encoded_index['first'],
        "A": encoded_index['A'],
        "C": encoded_index['C'],
        "G": encoded_index['G'],
        "T": encoded_index['T'],
        "rc_seeds": encoded_index['rc_seeds'],
        "runlength_output": encoded_index['runlength_output'],
        "runlength_value": encoded_index['runlength_value'],
        "runlength_rc_output": encoded_index['runlength_rc_output'],
        "runlength_rc_character": encoded_index['runlength_rc_character'],
        "runlength_rc_value": encoded_index['runlength_rc_value'],
        "org_location": encoded_index['org_location'],
    })

    # Extend the dataset
    extended_encoded_index = extend_dataset(df_training_data, config)

    # Save the extended encoded index
    extended_output_file_name = os.path.join(folder_path, f'S{i}_extended_encoded_reference.csv')
    save_index_as_csv(extended_encoded_index, extended_output_file_name)

reference_files = [f for f in os.listdir(folder_path) if f.endswith('extended_encoded_reference.csv')]

for ref_file in reference_files:
    # Data import
    file_path = os.path.join(folder_path, ref_file)
    df = pd.read_csv(file_path)
    df = df.dropna()
    df = df[df['org_location'] != 0.0]  # Handle the error with models treating 0.0 as null
    #df = df.drop(columns=['suffix_array', 'k-seed-extend'])

    # Preparing training and test data
    x = df.loc[:, df.columns != 'org_location']
    y = df.loc[:, 'org_location']
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, shuffle=True, random_state=42)
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, shuffle=True, random_state=42)

    # Saving the datasets in the same folder
    base_filename = ref_file.split('_reference.csv')[0]
    X_train.to_csv(os.path.join(folder_path, f"{base_filename}_x_train.csv"), index=False)
    X_test.to_csv(os.path.join(folder_path, f"{base_filename}_X_test.csv"), index=False)
    y_train.to_csv(os.path.join(folder_path, f"{base_filename}_y_train.csv"), index=False)
    y_test.to_csv(os.path.join(folder_path, f"{base_filename}_y_test.csv"), index=False)
    X_val.to_csv(os.path.join(folder_path, f"{base_filename}_X_val.csv"), index=False)
    y_val.to_csv(os.path.join(folder_path, f"{base_filename}_y_val.csv"), index=False)

    print(f"Processed {ref_file}:")
    print(f"  X training dataset length: {len(X_train)}")
    print(f"  y training dataset length: {len(y_train)}")
