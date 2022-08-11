from msilib import sequence
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import glob
import sys
import os

def batch_split(batch_size, sequence_string, sequence_length):
    # '''Generate reference batches for large ref string processesing'''
    # number_of_batches = sequence_length/batch_size
    # print("Total number of temp batches to generate are: " + str(number_of_batches))
    # batches = list(map(''.join, zip(*[iter(sequence_string)]*batch_size))) #generate batches from sequence - does this inc tails?
    
    # path = r"gnb_seed_generation_2\data_sets\_temp\reference_" # temp dir file path
    # for i in range(len(batches)):
    #     filename = path  + "batch_" + str(i) + ".txt"
    #     file = open(filename, "w")
    #     file.write(batches[i] + "\n")
    #     file.close()
    pass

class IO_processing:
    '''main IO and data preprocessing class for pharsing fasta data.'''
    def pharse_reference(self, input_file, batch_size):
        '''Pharse .fa file to extract genomic string data'''
        print("... Pharsing input reference genome")
        try:
            for seq_record in SeqIO.parse(input_file, "fasta"):
                print("\npharsing reference sequence " + input_file)
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

    def save_model_weights(self, model, input_file):
        '''Function to save model training weights as .pkl'''
        import pickle
        model_name = input_file + ".pkl"
        file_path = "pyseed/model_weights/"
        pkl_filename = file_path + model_name
        with open(pkl_filename, 'wb') as file:
            pickle.dump(model, file)
        return 0

    def weights_lookup(self, model_training_weights, training_flag):
        '''Function to return list of pretrained weights and trigger retraining in unknown'''
        print("... Looking up available model weights")
        if os.path.isfile('model_weights/' + model_training_weights + '.pkl'):
            training_flag = False
        else:
            training_flag = True
            print("Training weights not found, activating model training")
        return training_flag
                


# def pharse_reads():
#     count = 0
#     total_len = 0
#     with open("ls_orchid.fasta") as in_handle:
#         for title, seq in SimpleFastaParser(in_handle):
#             count += 1
#             total_len += len(seq)

# print("%i records with total sequence length %i" % (count, total_len))


