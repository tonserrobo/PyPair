import sys
import logging

import warnings

warnings.filterwarnings('ignore')  # ignore f1score warning on zero dev

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_fscore_support as f1_score
from sklearn.naive_bayes import GaussianNB

# Importing the required libraries
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate
from numpy import arange
from sklearn.model_selection import cross_val_score

# import class libs
from gen_index import GenIndex
from gen_training_data import GenTraining
from encode_decode import Encode_decode
from io_processing import IO_processing

from numpy import mean
from numpy import std

config = {
    'seed_length': 20,
    'read_length': 100,
    'training_iterations': 15,
}

logging.basicConfig(level=logging.INFO)


# def test_reads_generation(string, read_length):
#     '''Function to generate test reads from the reference sequence - for debug only'''
#     reads = []
#     for i in range(len(string) - read_length):
#         read_end = i + read_length
#         temp_read = string[i : read_end]
#         reads.append(temp_read)
#     test_reads = pd.DataFrame(reads)
#     test_reads.to_csv('data_sets/test_reads.csv', index = False)

def index_and_train(sequence_string, seed_length, training_iterations, read_length, model):
    """Function to issue genome indexing, training data preparation and model training"""

    logging.info('...Indexing reference genome and training model: Started')

    # generate index
    index_obj = GenIndex()
    index_obj.generate_index(sequence_string, seed_length)  # maybe this should come first.

    # generate training data
    training_data_obj = GenTraining()
    training_data = training_data_obj.generate_training_data(training_iterations, seed_length,
                                                             read_length)  # will gen index if does not already exsist
    df = training_data
    x = df.loc[:, df.columns != 'position_label']
    y = df.loc[:, 'position_label']
    X_train, X_test, y_train, y_test = train_test_split(x, y, shuffle=True, test_size=0.2)
    X_Val, X_test2, y_val, y_test2 = train_test_split(X_test, y_test, shuffle=True,
                                                      test_size=0.2)  # takes a further 20% as blind

    ## for training loop debugging and profiling
    import wandb
    wandb.init(
        project="gnb_seed_generation",
        notes="reference_index as training data",
        tags=["training data modification", "paper2"]
    )
    # model training
    model.fit(X_train, y_train)
    wandb.log({'Training Score': model.score(X_train, y_train), 'Test Score': model.score(X_test, y_test),
               'Classification report': classification_report(y_test, model.predict(X_test))})

    '''
    Training Evaluation
    '''
    # Printing Accuracy on Training and Test sets
    print(f"Training Set Score : {model.score(X_train, y_train) * 100} %")
    print(f"Test Set Score : {model.score(X_test, y_test) * 100} %")
    print(f"Val Set Score : {model.score(X_Val, y_val) * 100} %")

    # cross_validate also allows to specify metrics which you want to see
    for i, score in enumerate(cross_validate(model, X_test, y_test, cv=3)["test_score"]):
        print(f"Accuracy for the fold no. {i} on the test set: {score}")
        # wandb.log({'K-Fold': i, 'k-Fold CV': score })

    # evaluate model
    predicted = model.predict(X_test)
    scores = cross_val_score(model, X_test, y_test, scoring='accuracy', cv=3, n_jobs=-1)
    print('Accuracy: %.3f (%.3f)' % (mean(scores), std(scores)))
    wandb.log({'K-Fold mean': mean(scores),
               'k-Fold std': std(scores)
               })
    precision, recall, fscore, support = f1_score(y_test, predicted)
    wandb.log({'Precision': mean(precision),
               'Recall': mean(recall),
               'F1 Score': mean(fscore),
               'Class Support min': min(support),
               'Class Support Max': max(support),
               'Class Support mean:': mean(support)
               })
    wandb.finish()
    logging.info('...Indexing reference genome and training model: Compleate')
    return model


def main():
    """ for commonly used genome reference such as human, fish, mouse, rat - need something here to append the name of the ref used  """
    # config 
    seed_length = config['seed_length']
    read_length = config['read_length']
    training_iterations = config[
        'training_iterations']  # 30 try 10, then 20 then 30 to find which is giving the best tradeoff. This can be added to analysis
    training_flag = False

    # imput arguments
    if len(sys.argv) > 1:
        reference_genome_file = sys.argv[0]
        model_training_weights = sys.argv[0]
    else:
        reference_genome_file = 'reference_samples/s1.fasta'
        model_training_weights = "mouse"

    data = IO_processing()
    batch_size = len(reference_genome_file)  # using full length as batch
    sequence_id, sequence_string, sequence_length = data.pharse_reference(reference_genome_file, batch_size)
    training_flag = data.weights_lookup(model_training_weights, training_flag)

    # Gen test reads 
    # test_reads_generation(sequence_string, read_length)

    # model assignment
    model = GaussianNB()

    if training_flag:
        logging.info('Training flag active, please wait...')
        trained_model = index_and_train(sequence_string, seed_length, training_iterations, read_length, model)
        data.save_model_weights(trained_model, model_training_weights)
        logging.info('Model trained')
    else:
        print("this should just call predict")


if __name__ == '__main__':
    main()
