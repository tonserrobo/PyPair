import pandas as pd


# works fine
def sequence_variation(sequence):
    """ add single character variation into the sequence string """
    sequence_variants = []
    options = ["A", "T", "C", "G"]
    for count, char in enumerate(sequence):
        temp_list = list(sequence.strip())
        for j in range(len(options)):
            character = options[j]
            temp_list[count] = character
            new_string = "".join(temp_list)
            sequence_variants.append(new_string)
    return sequence_variants


sequences = pd.read_csv('data_sets/indexed_reference.csv')
#alignment_locations = sequences['org_location']
#seeds = sequences['k-seed']
#total_length = len(seeds[0] * 4)
#
# seed_var = []
# alignment_location_var = []

for seed in range(sequences):
    seed_var = []
    print(seed)
    _alignment_location = seed['org_location']
    _seq = seed['k-seed']

    _seeds = sequence_variation(_seq)
    #seed_var.extend(_seeds)

print(seed_var)

# print("Total var per sequence: " + str(total_length))
# print("Total var list len expected: " + str(total_length * len(seeds)))
# print("Total var list produced: " + str(len(seed_var)))

# apply label to correct seed variation.
