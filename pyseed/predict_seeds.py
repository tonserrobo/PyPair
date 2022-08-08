import pickle
import pandas as pd 
import numpy as np

'''
Load synth short reads from neat gen reads on a given genome. 
Load genome weights as pkl 
predict all synth short reads locations 
val function to look up predicted locations on ref string and iterate through each <= seed length 
    could do this in the left + right direction to create list of mems/SMEMs 
    if mismatch > 2 drop seed. therefore threshold == 0.92 min to qualify valid seed for > 20bp seed 

need to populate this into a data frame for further analysis. 
i.e. 
location, seed_no, match_accuracy, min_bp_mapped, max_bp_mapped, 

psudo for MEM_val:
def val(min_seed_length, mismatch_threshold, seed, map_location, seed_length):
    total_mismatch = 0
    total_bp_mapped = 0

    # all parts of the dataframe to analyse
    seed_location_start = []
    seed_location_end = []
    seed_length = []
    seed_mismatches = []
    seed_bp_mapped = []

    seed_start = map_location
    seed_end = map_location + seed_length
    map_location_sequence = reference[seed_start:seed_end]


    if len(map_location_sequence) == len(seed):
        while (total_mismatch < mismatch_threshold): # prob shouldnt be a while loop 
            for i in range(len(seed)):
                ref_char = map_location_sequence[i]
                seed_char = seed[i]
                if ref_char == seed_char:
                    # match 
                    total_bp_mapped += 1
                else:
                    total_mismatch += 1 
            





'''
def load_test(file):
    model_path = r"C:\Users\tonyr\OneDrive - Ulster University\Phd - FPGA acceleration in genomics\2020-21 MainResearch\Seed_generation\rapid-SMEM\Rapid-SMEM\saved_models\Centoid.pkl"
    #model_filename = file + ".pkl"
    read = '1' #k=3
    read_s = np.array(list(read), dtype=int)
    resaped_read = read_s.reshape(1, -1)
    correct_res = '15'
    pkl_filename = model_path #+ model_filename
    return pkl_filename, resaped_read

def predict(pickle_model, read):
    Ypredict = pickle_model.predict(read)
    return Ypredict

def main():
    files = ['dt', 'knn']
    filename, read = load_test('Centoid')
    # Load from file
    with open(filename, 'rb') as file:
        model = pickle.load(file)

    res = predict(model, read)
    print(res)

if __name__ == '__main__':
    main()