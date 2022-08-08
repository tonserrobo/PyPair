from msilib import sequence
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import glob
import sys
import os

@profile
def batch_split(batch_size, sequence_string, sequence_length):
    '''Generate reference batches for large ref string processesing'''
    number_of_batches = sequence_length/batch_size
    print("Total number of temp batches to generate are: " + str(number_of_batches))
    batches = list(map(''.join, zip(*[iter(sequence_string)]*batch_size))) #generate batches from sequence - does this inc tails?
    
    path = r"gnb_seed_generation_2\data_sets\_temp\reference_" # temp dir file path
    for i in range(len(batches)):
        filename = path  + "batch_" + str(i) + ".txt"
        file = open(filename, "w")
        file.write(batches[i] + "\n")
        file.close()
    pass

'''main IO and data preprocessing class for pharsing fasta data.'''
@profile
def pharse_reference(input_file, batch_size):
    '''Pharse .fa file to extract genomic string data'''
    try:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            print("pharsing reference sequence " + input_file + "\n")
            # note this assumes there is a single seq_record in each file - currently, last record carries
            sequence_id = seq_record.id
            sequence_string = str(seq_record.seq)
            sequence_string = sequence_string.upper()
            sequence_length = len(seq_record)
            print("sequence ID: " + str(sequence_id))
            print("sequence Length: " + str(sequence_length))

            # generate batches based on batch size
            batch_split(batch_size, sequence_string, sequence_length)
    except UnboundLocalError:
        print("error pharsing input requence sequence, check input is fasta format")
    return sequence_id, sequence_string, sequence_length

@profile
def save_model_weights(model):
    '''Function to save model training weights as .pkl'''
    import pickle
    try:
        name = str(input("Please provide a name for these training weights: "))
    except:
        print("please define a name for the training file")
    #try:
    model_name = "\human.pkl"
    file_path = r"C:\Users\tonyr\OneDrive - Ulster University\Phd - FPGA acceleration in genomics\2020-21 MainResearch\Seed_generation\rapid-SMEM\Rapid-SMEM\gnb_seed_generation_2\code_profiling\sample2.fasta"
    pkl_filename = file_path + model_name
    with open(pkl_filename, 'wb') as file:
        pickle.dump(model, file)
    # except:
    #     print("error saving model weights....")
    return 0

@profile
def weights_lookup(model_training_weights, training_flag):
    '''Function to return list of pretrained weights and trigger retraining in unknown'''
    
    temp_dir_files = glob.glob("gnb_seed_generation_2\model_weights\weigths_*.pkl")
    print(len(temp_dir_files)) #test

    for file in temp_dir_files:
        print(file)
        print(model_training_weights)
        if file == model_training_weights:
            training_flag == False
        else:
            training_flag == True
    return training_flag

def main():
    batch_size = 1
    reference_genome_file = r"C:\Users\tonyr\OneDrive - Ulster University\Phd - FPGA acceleration in genomics\2020-21 MainResearch\Seed_generation\rapid-SMEM\Rapid-SMEM\gnb_seed_generation_2\code_profiling\5000bp_ref.txt"
    pharse_reference(reference_genome_file, batch_size)

if __name__=='__main__':
    main()


