import pandas as pd
import numpy as np
import time
import joblib
import warnings
import shap
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import precision_recall_fscore_support
from memory_profiler import memory_usage
from sklearn.model_selection import StratifiedShuffleSplit, GridSearchCV
from imblearn.over_sampling import SMOTE
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import NotFittedError
import matplotlib.pyplot as plt

# read training and test data
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

# Assuming datasets is a list of tuples with training and testing data
datasets = [(S1_X_train, S1_X_test, S1_y_train, S1_y_test),
            (S2_X_train, S2_X_test, S2_y_train, S2_y_test),
            (S3_X_train, S3_X_test, S3_y_train, S3_y_test),
            (S4_X_train, S4_X_test, S4_y_train, S4_y_test)]

aggregate_metrics = {'Precision': [],
                     'Recall': [],
                     'F1 Score': [],
                     'Memory (MB) 50k batch': [],
                     'Execution Time (s) 50k batch': [],
                     'Memory (MB) per read mapping': [],
                     'Execution Time (s) per read mapping': [],
                     'Mean CV F1 Score': []}

def profile_memory(_model, _X_test):
    _mem_usage = memory_usage((_model.predict, (_X_test,)), interval=0.1)
    return max(_mem_usage) - min(_mem_usage)

# Create a placeholder for the best model
best_model = None

for X_train, X_test, y_train, y_test in datasets:
    # Convert y_train and y_test to 1D array
    y_train = y_train.squeeze()
    y_test = y_test.squeeze()

    # Ensure all data is numeric
    X_train = pd.DataFrame(X_train).apply(pd.to_numeric, errors='coerce')
    X_test = pd.DataFrame(X_test).apply(pd.to_numeric, errors='coerce')
    y_train = pd.Series(y_train).apply(pd.to_numeric, errors='coerce')
    y_test = pd.Series(y_test).apply(pd.to_numeric, errors='coerce')

    # Handle NaN values
    X_train.ffill(inplace=True)
    X_train.dropna(inplace=True)
    X_test.ffill(inplace=True)
    X_test.dropna(inplace=True)

    y_train = pd.Series(y_train).ffill().dropna().values
    y_test = pd.Series(y_test).ffill().dropna().values

    # Handle infinite values
    X_train.replace([np.inf, -np.inf], np.nan, inplace=True)
    X_test.replace([np.inf, -np.inf], np.nan, inplace=True)
    y_train = np.where(np.isinf(y_train), np.nan, y_train)
    y_test = np.where(np.isinf(y_test), np.nan, y_test)

    X_train.ffill(inplace=True)
    X_train.dropna(inplace=True)
    X_test.ffill(inplace=True)
    X_test.dropna(inplace=True)

    y_train = pd.Series(y_train).ffill().dropna().values
    y_test = pd.Series(y_test).ffill().dropna().values

    # Ensure data consistency
    feature_names = X_train.columns.tolist()  # Save column names for SHAP
    X_train = X_train.values
    X_test = X_test.values
    y_train = y_train
    y_test = y_test

    # Feature Scaling
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Handle class imbalance with SMOTE
    min_samples = min([sum(y_train == cls) for cls in np.unique(y_train)])
    if min_samples > 1:
        k_neighbors = min(5, min_samples - 1)
        smote = SMOTE(k_neighbors=k_neighbors)
        X_train_resampled, y_train_resampled = smote.fit_resample(X_train_scaled, y_train)
    else:
        # If a class has only one sample, duplicate it to allow splitting
        class_counts = pd.Series(y_train).value_counts()
        for cls in class_counts[class_counts == 1].index:
            X_train_scaled = np.vstack([X_train_scaled, X_train_scaled[y_train == cls]])
            y_train = np.hstack([y_train, y_train[y_train == cls]])
        X_train_resampled, y_train_resampled = X_train_scaled, y_train

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
        'n_neighbors': [3, 5, 7],
        'metric': ['minkowski', 'euclidean', 'manhattan'],
        'weights': ['uniform', 'distance']
    }

    model = KNeighborsClassifier()

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

    # Measure memory usage during prediction
    try:
        mem_usage = profile_memory(best_model, X_test_scaled)
        if mem_usage == 0:
            mem_usage = np.nan  # Ensure no zero values
    except NotFittedError:
        mem_usage = np.nan

    # Measure execution time during prediction
    start_time = time.time()
    predicted = best_model.predict(X_test_scaled)
    execution_time = time.time() - start_time

    # Metrics Calculation
    precision, recall, fscore, support = precision_recall_fscore_support(y_test, predicted, average='weighted')

    execution_time_per_read = execution_time / len(X_test_scaled)  # per read mapping execution time
    memory_usage_per_read = mem_usage / len(X_test_scaled) if not np.isnan(mem_usage) else np.nan  # memory required per read mapping

    # Aggregating metrics
    aggregate_metrics['Precision'].append(precision)
    aggregate_metrics['Recall'].append(recall)
    aggregate_metrics['F1 Score'].append(fscore)
    aggregate_metrics['Mean CV F1 Score'].append(np.mean(cv_scores))
    aggregate_metrics['Memory (MB) 50k batch'].append(mem_usage)
    aggregate_metrics['Execution Time (s) 50k batch'].append(execution_time)
    aggregate_metrics['Memory (MB) per read mapping'].append(memory_usage_per_read)
    aggregate_metrics['Execution Time (s) per read mapping'].append(execution_time_per_read)

# Save the last trained KNN model
filename = 'models/model_weights/knn_model.joblib'
joblib.dump(best_model, filename)

# SHAP Explainer
explainer = shap.KernelExplainer(best_model.predict, X_train_resampled[:100])  # Use a subset of the training data
shap_values = explainer.shap_values(X_train_resampled[:100])  # Use a subset of the training data

# SHAP summary plot with feature names
shap.summary_plot(shap_values,
                  X_train_resampled[:100],
                  feature_names=feature_names,
                  show=False)
plt.savefig("results/knn_shap_values.png", dpi=300, bbox_inches='tight')
plt.close()

KNeighborsClassifier_shap_values = pd.DataFrame(shap_values)
KNeighborsClassifier_shap_values['Classifier'] = 'KNeighborsClassifier'

KNeighborsClassifier_results = pd.DataFrame(aggregate_metrics)
KNeighborsClassifier_results['Model'] = 'KNeighborsClassifier'
KNeighborsClassifier_results['Sample'] = ['S1', 'S2', 'S3', 'S4']
KNeighborsClassifier_results = KNeighborsClassifier_results[['Model',
                                                             'Sample',
                                                             'Mean CV F1 Score',
                                                             'Precision',
                                                             'Recall',
                                                             'F1 Score',
                                                             'Memory (MB) per read mapping',
                                                             'Execution Time (s) per read mapping',
                                                             'Memory (MB) 50k batch',
                                                             'Execution Time (s) 50k batch']]

KNeighborsClassifier_results.to_csv("results/KNeighborsClassifier_results.csv", index=False)
KNeighborsClassifier_shap_values.to_csv("results/KNeighborsClassifier_shap_values.csv", index=False)

del best_model, X_test, X_train, y_train, y_test