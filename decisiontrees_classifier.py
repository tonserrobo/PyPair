import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import precision_recall_fscore_support, accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit, GridSearchCV
from imblearn.over_sampling import SMOTE
from memory_profiler import memory_usage
import joblib
import time
import gc
import warnings
import shap

# read training & Test set
folder_path = 'data/reference_samples/'
S1_X_train = pd.read_csv(folder_path + 'S1_extended_encoded_X_train.csv')
S1_X_test = pd.read_csv(folder_path + 'S1_extended_encoded_X_test.csv')
S1_y_train = pd.read_csv(folder_path + 'S1_extended_encoded_y_train.csv')
S1_y_test = pd.read_csv(folder_path + 'S1_extended_encoded_y_test.csv')

S2_X_train = pd.read_csv(folder_path + 'S2_extended_encoded_X_train.csv')
S2_X_test = pd.read_csv(folder_path + 'S2_extended_encoded_X_test.csv')
S2_y_train = pd.read_csv(folder_path + 'S2_extended_encoded_y_train.csv')
S2_y_test = pd.read_csv(folder_path + 'S2_extended_encoded_y_test.csv')

S3_X_train = pd.read_csv(folder_path + 'S3_extended_encoded_X_train.csv')
S3_X_test = pd.read_csv(folder_path + 'S3_extended_encoded_X_test.csv')
S3_y_train = pd.read_csv(folder_path + 'S3_extended_encoded_y_train.csv')
S3_y_test = pd.read_csv(folder_path + 'S3_extended_encoded_y_test.csv')

S4_X_train = pd.read_csv(folder_path + 'S4_extended_encoded_X_train.csv')
S4_X_test = pd.read_csv(folder_path + 'S4_extended_encoded_X_test.csv')
S4_y_train = pd.read_csv(folder_path + 'S4_extended_encoded_y_train.csv')
S4_y_test = pd.read_csv(folder_path + 'S4_extended_encoded_y_test.csv')


# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# Function to preprocess the data
def preprocess_data(X_train, X_test):
    # Convert all values to numerical types, replacing non-numerical values with NaN
    X_train = pd.DataFrame(X_train).apply(pd.to_numeric, errors='coerce')
    X_test = pd.DataFrame(X_test).apply(pd.to_numeric, errors='coerce')

    # Replace NaNs, infinities, and too large values
    X_train = X_train.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    X_test = X_test.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # Scale the data
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    return X_train, X_test

# Function to process data in batches
def process_in_batches(model, X_test, batch_size=500):
    predictions = []
    total_memory = 0
    total_time = 0

    for start in range(0, len(X_test), batch_size):
        end = start + batch_size
        batch = X_test[start:end]

        # Measure memory usage for this batch
        try:
            mem_usage = memory_usage((model.predict, (batch,)), interval=0.1)
            batch_mem_usage = max(mem_usage) - min(mem_usage)
            if batch_mem_usage == 0:
                batch_mem_usage = np.nan  # Ensure no zero values
        except Exception as e:
            batch_mem_usage = np.nan
        total_memory += batch_mem_usage

        # Measure execution time for this batch
        start_time = time.time()
        preds = model.predict(batch)
        total_time += time.time() - start_time

        predictions.extend(preds)

        # Explicitly clear memory
        del batch
        gc.collect()

    return np.array(predictions), total_memory, total_time

# Assuming datasets is a list of tuples with training and testing data
datasets = [(S1_X_train, S1_X_test, S1_y_train, S1_y_test),
            (S2_X_train, S2_X_test, S2_y_train, S2_y_test),
            (S3_X_train, S3_X_test, S3_y_train, S3_y_test),
            (S4_X_train, S4_X_test, S4_y_train, S4_y_test)]

aggregate_metrics = {'Precision': [],
                     'Recall': [],
                     'F1 Score': [],
                     'Accuracy': [],
                     'Memory (MB) 50k batch': [],
                     'Execution Time (s) 50k batch': [],
                     'Memory (MB) per read mapping': [],
                     'Execution Time (s) per read mapping': [],
                     'Mean CV F1 Score': []}

# Create a placeholder for the best model
best_model = None

for X_train, X_test, y_train, y_test in datasets:
    # Preprocess the data
    X_train, X_test = preprocess_data(X_train, X_test)

    # Ensure data is in float32 format
    X_train = X_train.astype(np.float32)
    X_test = X_test.astype(np.float32)

    # Convert y_train and y_test to 1D array
    y_train = y_train.squeeze()
    y_test = y_test.squeeze()

    # Handle class imbalance with SMOTE
    min_samples = min([sum(y_train == cls) for cls in np.unique(y_train)])
    if min_samples > 1:
        k_neighbors = min(5, min_samples - 1)
        smote = SMOTE(k_neighbors=k_neighbors)
        X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
    else:
        # If a class has only one sample, duplicate it to allow splitting
        class_counts = pd.Series(y_train).value_counts()
        for cls in class_counts[class_counts == 1].index:
            X_train = np.vstack([X_train, X_train[y_train == cls]])
            y_train = np.hstack([y_train, y_train[y_train == cls]])
        X_train_resampled, y_train_resampled = X_train, y_train

    # Ensure the resampled data is in numpy array format
    X_train_resampled = np.array(X_train_resampled)
    y_train_resampled = np.array(y_train_resampled)

    # Cross-validation to detect overfitting
    sss = StratifiedShuffleSplit(n_splits=5, test_size=0.2)

    # Adjust the test_size to avoid errors
    if len(np.unique(y_train_resampled)) > int(0.2 * len(y_train_resampled)):
        sss = StratifiedShuffleSplit(n_splits=5, test_size=len(np.unique(y_train_resampled)))

    # Hyperparameter tuning
    param_grid = {
        'max_depth': [None, 10, 20],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4],
        'criterion': ['gini', 'entropy']
    }

    model = DecisionTreeClassifier(random_state=0)

    try:
        grid_search = GridSearchCV(model, param_grid, scoring='f1_weighted', cv=sss)
        grid_search.fit(X_train_resampled, y_train_resampled)
        best_model = grid_search.best_estimator_
    except Exception as e:
        print(f"Grid search failed: {e}")
        best_model = model

    cv_scores = []
    for train_index, test_index in sss.split(X_train_resampled, y_train_resampled):
        X_cv_train, X_cv_test = X_train_resampled[train_index], X_train_resampled[test_index]
        y_cv_train, y_cv_test = y_train_resampled[train_index], y_train_resampled[test_index]

        best_model.fit(X_cv_train, y_cv_train)
        y_cv_pred = best_model.predict(X_cv_test)
        _, _, f1, _ = precision_recall_fscore_support(y_cv_test, y_cv_pred, average='weighted')
        cv_scores.append(f1)

    print(f"Cross-validation F1 scores: {cv_scores}")
    print(f"Mean CV F1 score: {np.mean(cv_scores)}")

    best_model.fit(X_train_resampled, y_train_resampled)

    # Process in batches to handle memory usage
    predicted, mem_usage, execution_time = process_in_batches(best_model, X_test, batch_size=500)

    # Metrics Calculation
    precision, recall, fscore, support = precision_recall_fscore_support(y_test, predicted, average='weighted')
    accuracy = accuracy_score(y_test, predicted)

    execution_time_per_read = execution_time / len(X_test)  # per read mapping execution time
    memory_usage_per_read = mem_usage / len(X_test)  # memory required per read mapping

    # Aggregating metrics
    aggregate_metrics['Precision'].append(precision)
    aggregate_metrics['Recall'].append(recall)
    aggregate_metrics['F1 Score'].append(fscore)
    aggregate_metrics['Accuracy'].append(accuracy)
    aggregate_metrics['Mean CV F1 Score'].append(np.mean(cv_scores))
    aggregate_metrics['Memory (MB) 50k batch'].append(mem_usage)
    aggregate_metrics['Execution Time (s) 50k batch'].append(execution_time)
    aggregate_metrics['Memory (MB) per read mapping'].append(memory_usage_per_read)
    aggregate_metrics['Execution Time (s) per read mapping'].append(execution_time_per_read)

# Save the last trained Decision Tree model
filename = 'models/model_weights/decision_tree_model.joblib'
joblib.dump(best_model, filename)

# SHAP Explainer
explainer = shap.TreeExplainer(best_model)
shap_values = explainer.shap_values(X_train_resampled[:100])  # Use a subset of the training data

# SHAP summary plot with feature names
feature_names = [f"Feature {i}" for i in range(X_train_resampled.shape[1])]
shap.summary_plot(shap_values,
                  X_train_resampled[:100],
                  feature_names=feature_names,
                  show=False)
plt.savefig("results/decision_tree_shap_values.png", dpi=300, bbox_inches='tight')
plt.close()

DecisionTreeClassifier_shap_values = pd.DataFrame(shap_values)
DecisionTreeClassifier_shap_values['Classifier'] = 'DecisionTreeClassifier'

DecisionTreeClassifier_results = pd.DataFrame(aggregate_metrics)
DecisionTreeClassifier_results['Model'] = 'DecisionTreeClassifier'
DecisionTreeClassifier_results['Sample'] = ['S1', 'S2', 'S3', 'S4']
DecisionTreeClassifier_results = DecisionTreeClassifier_results[['Model',
                                                                 'Sample',
                                                                 'Mean CV F1 Score',
                                                                 'Precision',
                                                                 'Recall',
                                                                 'F1 Score',
                                                                 'Memory (MB) per read mapping',
                                                                 'Execution Time (s) per read mapping',
                                                                 'Memory (MB) 50k batch',
                                                                 'Execution Time (s) 50k batch']]

DecisionTreeClassifier_results.to_csv("results/decisionTreeClassifier_results.csv", index=False)
DecisionTreeClassifier_shap_values.to_csv("results/DecisionTreeClassifier_shap_values.csv", index=False)

del best_model, X_test, X_train, y_train, y_test
gc.collect()