import pandas as pd 
import os
from encode_decode import Encode_decode

@profile
def generate_labeled_training_data(training_data, repetitions):
        '''
        Generate training data based on repeating the number of occurances of k-seed and its assocated org_location. 
        '''
        seed = []
        label = []
        
        for i in range(0,repetitions):
            for i in training_data['k-seed']:
                seed.append(i)
            for i in training_data['org_location']:
                label.append(i)
                
        data = pd.DataFrame()
        data['k-seed'] = seed
        data['position_label'] = label
        del seed
        del label
        return data

@profile
def generate_training_data_full_read_map(training_data, repetitions):
    '''
    Generate training data based on repeating the number of occurances of k-seed and its assocated org_location. 
    '''
    seed = []
    
    for i in range(0,repetitions):
        for i in training_data['k-seed-extend']:
            seed.append(i)
            
    data = pd.DataFrame()
    data['k-seed-extend'] = seed
    del seed
    return data

@profile
def generate_training_data_2(training_data, repetitions):
    """
    Simplified version of generate training data, in which the entire set is repeated in a single instance
    - this method will become the defult for batch processing 
    """
    data = []
    temp_range = len(training_data)
    batch_size = 5000 #this should prob be passed via main as part of the main options
    total_data_size = len(training_data)
    number_of_batches = round(total_data_size/batch_size)

    for i in range(0, repetitions):
        data_temp = training_data[0:temp_range]
        data.append(data_temp)

    comp_data = pd.DataFrame(data)
    return comp_data

@profile
def reverse_k_seed_string(training_data):
    """
    Method to reverse k-seed string and include within training set as additional feature
    """
    data = [] 
    for i in training_data['k-seed']:
        rev_seed = i[::-1]
        data.append(rev_seed)

    training_data['rev_k-seed'] = data
    del data
    return training_data

@profile
def reverse_k_seed_string_full_map(training_data):
    """
    Method to reverse k-seed string and include within training set as additional feature
    """
    data = [] 
    for i in training_data['k-seed-extend']:
        rev_seed = i[::-1]
        data.append(rev_seed)

    training_data['rev_k-seed-extend'] = data
    del data
    return training_data

"""
Training data generation class and method returns a shuffled dataframe of k-seeds with position label, repeated x times.
"""
@profile
def generate_training_data(reference_genome, training_iterations):
    """
    Function call to GenTraining() class for training data generation. The action of generating training data will change the training
    flag to True, causing the model to retrain for the given training set. 
    """
    # Read/Generate reference index file
    try:
        index = pd.read_csv("data_sets/indexed_reference.csv")
    except FileNotFoundError:
        print("indexed_reference.csv file not found") #not calling index gen from here
        #index = generate_index(reference_genome) 

    encoder_decoder = Encode_decode()
    
    # Generate and encode foward training data
    training_set = generate_labeled_training_data(index, training_iterations)
    training_k_seeds = list(training_set['k-seed'])
    training_set['k-seed'] = encoder_decoder.encode(training_k_seeds) # 2bit encoding

    # Generate and encode reverse seed string
    training_set = reverse_k_seed_string(training_set)
    training_k_seeds_rev = list(training_set['rev_k-seed'])
    training_set['rev_k-seed'] = encoder_decoder.encode(training_k_seeds_rev)

    #for full read mapping - fwd
    training_set_full_map_fwd = generate_training_data_full_read_map(index, training_iterations)
    training_k_seeds_extend_fwd = list(training_set_full_map_fwd['k-seed-extend'])
    training_set['k-seed-extend'] = encoder_decoder.encode(training_k_seeds_extend_fwd) # 2bit encoding

    #for full read mapping - rev
    training_set_full_map_rev = reverse_k_seed_string_full_map(training_set)
    training_k_seeds_extend_rev = list(training_set_full_map_rev['rev_k-seed-extend'])
    training_set['rev_k-seed-extend'] = encoder_decoder.encode(training_k_seeds_extend_rev)
    try:
        training_set.to_csv("data_sets/training_data.csv", index = False)
    except OSError:
        print("Cannot find data_sets directory, creating one...")
        os.makedirs('data_sets')
        training_set.to_csv("data_sets/training_data.csv", index = False)

    training_flag = True
    return training_flag, training_set

@profile
def main():
    training_iterations = 4   
    reference_genome_file = r"C:\Users\tonyr\OneDrive - Ulster University\Phd - FPGA acceleration in genomics\2020-21 MainResearch\Seed_generation\rapid-SMEM\Rapid-SMEM\gnb_seed_generation_2\code_profiling\5000bp_ref.txt"
    training_flag, training_data = generate_training_data(reference_genome_file, training_iterations)

if __name__=='__main__':
    main()
