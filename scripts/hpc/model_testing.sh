# do we need to create an initial job?

## import and install libraries
# activate conda envrionment 'pypair', create if not already exsist
# install packages detailed within requirements.txt to 'pypair' env
# confirm that env is created

# data set preparation
echo "running multimodel_training.py"
multimodel_training.py
echo "multimodel_training.py compleate, training and test data generated"

# model tests
echo "beginning model test"

# Nearest neighbour models
echo "testing nearest neighbour model archetype"

# Nearest centroid model
echo "Nearest centroid model training job"
# create a new job
# run nearest centroid model training
nearest_centroid.py
echo "Nearest centroid model training complete"

# random forest model
echo "Random forest model training"
# create new job
# run random_forest.py within the new job -> actually will need to ensure page
random_forest.py
echo "Random forest model training complete"

nearest_centroid.py
kNeighbours_classifier.py
random_forest.py
extratrees_classifier.py
decisiontrees_classifier.py
gaussianNaiveBayes_classifier.py
LGBM_classifier.py
XGB_classifier.py

archetype_aggregation.py
