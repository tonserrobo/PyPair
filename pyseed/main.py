import sys

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

def main():
    """ for commonly used genome reference such as human, fish, mouse, rat - need something here to append the name of the ref used  """
    # config 
    seed_length = 28 #defult seed length as used by BWA-MEM
    training_iterations = 30 #looks optimum from sweeps
    training_flag = False # used to trigger model retrain

    # imput arguments
    if len(sys.argv) > 1:
        reference_genome_file = sys.argv[0]
        model_training_weights = sys.argv[0]
    else:
        reference_genome_file = 'reference_samples/s1.fasta'
        model_training_weights = "mouse"
    
    print(len(sys.argv))
    print(str(sys.argv))

    # # pharse reference sequence
    data = IO_processing()
    batch_size = len(reference_genome_file) #using full length as batch
    sequence_id, sequence_string, sequence_length = data.pharse_reference(reference_genome_file, batch_size) #pharse reference sequence
    
    # # model weights lookup needs to happen here, or first just 
    training_flag = data.weights_lookup(model_training_weights, training_flag)

    # model assignment
    model = GaussianNB() #model def

    import warnings
    warnings.filterwarnings('ignore') #ignore f1score warning on zero dev

    if training_flag:
        # generate index
        index_obj = GenIndex() 
        index_obj.generate_index(sequence_string, seed_length) #maybe this should come first. 

        # generate training data
        training_data_obj = GenTraining()
        
        ## for training loop debugging and profiling
        import wandb
        wandb.init(
        project="gnb_seed_generation",
        notes="reference_index as training data",
        tags=["training data modification", "paper2"]
        )

        training_data = training_data_obj.generate_training_data(reference_genome_file, training_iterations) #will gen index if does not already exsist
 
        df = training_data
        x = df.loc[:,df.columns != 'position_label']
        y = df.loc[:, 'position_label']

        # Training/Test/Val split (80/20)
        X_train, X_test, y_train, y_test = train_test_split(x,y, shuffle=True, test_size=0.2)
        X_Val, X_test2, y_val, y_test2 = train_test_split(X_test, y_test, shuffle=True, test_size=0.2) #takes a further 20% as blind

        # model training
        model.fit(X_train, y_train)
        wandb.log({'Training Score': model.score(X_train, y_train), 'Test Score': model.score(X_test, y_test), 'Classification report': classification_report(y_test, model.predict(X_test)) })

        '''
        Training Evaluation
        '''
        # Printing Accuracy on Training and Test sets
        print(f"Training Set Score : {model.score(X_train, y_train) * 100} %")
        print(f"Test Set Score : {model.score(X_test, y_test) * 100} %")

        # cross_validate also allows to specify metrics which you want to see
        for i, score in enumerate(cross_validate(model, X_test,y_test, cv=3)["test_score"]):
            print(f"Accuracy for the fold no. {i} on the test set: {score}")
            #wandb.log({'K-Fold': i, 'k-Fold CV': score })

        # evaluate model
        predicted = model.predict(X_test)
        scores = cross_val_score(model, X_test, y_test, scoring='accuracy', cv=3, n_jobs=-1)
        print('Accuracy: %.3f (%.3f)' % (mean(scores), std(scores)))
        wandb.log({'K-Fold mean': mean(scores), 'k-Fold std': std(scores) })
        precision, recall, fscore, support = f1_score(y_test, predicted)
        wandb.log({'Precision': mean(precision), 'Recall': mean(recall), 'F1 Score':mean(fscore), 'Class Support min':min(support), 'Class Support Max': max(support), 'Class Support mean:':mean(support)})
        wandb.finish()

        data.save_model_weights(model) #save model weights usign IOclass

    else:
        print("this should just call predict")

if __name__ == '__main__':
    main()