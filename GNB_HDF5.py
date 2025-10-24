import matplotlib.pyplot as plt
import numpy as np
import vaex  # <-- Vaex for large/out-of-core data
import time
import gc
import warnings
import shap
from sklearn.metrics import precision_recall_fscore_support, accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit, GridSearchCV
from imblearn.over_sampling import SMOTE
from memory_profiler import memory_usage
import joblib
import pandas as pd  # still used for some transformations and final data frames

# We need these imports to replicate scikit-learn's interface
from sklearn.base import BaseEstimator, ClassifierMixin

# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

###############################################################################
# 1. Implement a Custom Gaussian Naive Bayes Class (no change from pandas script)
###############################################################################
class MyGaussianNB(BaseEstimator, ClassifierMixin):
    """
    A custom implementation of Gaussian Naive Bayes, emulating scikit-learn's
    interface so it can be used with GridSearchCV, etc.
    """
    def __init__(self, var_smoothing=1e-9):
        self.var_smoothing = var_smoothing

    def fit(self, X, y):
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y)

        self.classes_, y_indices = np.unique(y, return_inverse=True)
        n_classes = len(self.classes_)
        n_features = X.shape[1]

        self.theta_ = np.zeros((n_classes, n_features), dtype=np.float64)
        self.sigma_ = np.zeros((n_classes, n_features), dtype=np.float64)
        self.class_count_ = np.zeros(n_classes, dtype=np.float64)

        for idx, c in enumerate(self.classes_):
            X_c = X[y_indices == idx, :]
            self.class_count_[idx] = X_c.shape[0]
            self.theta_[idx, :] = np.mean(X_c, axis=0)
            self.sigma_[idx, :] = np.var(X_c, axis=0)

        self.class_prior_ = self.class_count_ / self.class_count_.sum()
        epsilon = self.var_smoothing * np.var(X, axis=0).max()
        self.sigma_ += epsilon

        return self

    def predict(self, X):
        X = np.asarray(X, dtype=np.float64)
        jll = self._joint_log_likelihood(X)
        return self.classes_[np.argmax(jll, axis=1)]

    def _joint_log_likelihood(self, X):
        n_samples, n_features = X.shape
        n_classes = len(self.classes_)
        log_likelihood = np.zeros((n_samples, n_classes), dtype=np.float64)

        log_prior = np.log(self.class_prior_)
        var = self.sigma_
        mean = self.theta_

        for c in range(n_classes):
            log_prob_const = -0.5 * np.sum(np.log(2.0 * np.pi * var[c]))
            log_prob_x = -0.5 * np.sum(((X - mean[c]) ** 2) / (var[c]), axis=1)
            log_likelihood[:, c] = log_prior[c] + log_prob_const + log_prob_x

        return log_likelihood

    def get_params(self, deep=True):
        return {"var_smoothing": self.var_smoothing}

    def set_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)
        return self


###############################################################################
# 2. Data Loading with Vaex
###############################################################################
folder_path = 'data/reference_samples/'

# Example: from_csv parameters can vary depending on delimiter, chunk size, etc.
# For large files, set convert=True to create a .hdf5 file that loads faster next time.

S1_X_train_vaex = vaex.from_csv(folder_path + 'S1_extended_encoded_X_train.csv', convert=False)
S1_X_test_vaex = vaex.from_csv(folder_path + 'S1_extended_encoded_X_test.csv', convert=False)
S1_y_train_vaex = vaex.from_csv(folder_path + 'S1_extended_encoded_y_train.csv', convert=False)
S1_y_test_vaex = vaex.from_csv(folder_path + 'S1_extended_encoded_y_test.csv', convert=False)

S2_X_train_vaex = vaex.from_csv(folder_path + 'S2_extended_encoded_X_train.csv', convert=False)
S2_X_test_vaex = vaex.from_csv(folder_path + 'S2_extended_encoded_X_test.csv', convert=False)
S2_y_train_vaex = vaex.from_csv(folder_path + 'S2_extended_encoded_y_train.csv', convert=False)
S2_y_test_vaex = vaex.from_csv(folder_path + 'S2_extended_encoded_y_test.csv', convert=False)

S3_X_train_vaex = vaex.from_csv(folder_path + 'S3_extended_encoded_X_train.csv', convert=False)
S3_X_test_vaex = vaex.from_csv(folder_path + 'S3_extended_encoded_X_test.csv', convert=False)
S3_y_train_vaex = vaex.from_csv(folder_path + 'S3_extended_encoded_y_train.csv', convert=False)
S3_y_test_vaex = vaex.from_csv(folder_path + 'S3_extended_encoded_y_test.csv', convert=False)

S4_X_train_vaex = vaex.from_csv(folder_path + 'S4_extended_encoded_X_train.csv', convert=False)
S4_X_test_vaex = vaex.from_csv(folder_path + 'S4_extended_encoded_X_test.csv', convert=False)
S4_y_train_vaex = vaex.from_csv(folder_path + 'S4_extended_encoded_y_train.csv', convert=False)
S4_y_test_vaex = vaex.from_csv(folder_path + 'S4_extended_encoded_y_test.csv', convert=False)

datasets_vaex = [
    (S1_X_train_vaex, S1_X_test_vaex, S1_y_train_vaex, S1_y_test_vaex),
    (S2_X_train_vaex, S2_X_test_vaex, S2_y_train_vaex, S2_y_test_vaex),
    (S3_X_train_vaex, S3_X_test_vaex, S3_y_train_vaex, S3_y_test_vaex),
    (S4_X_train_vaex, S4_X_test_vaex, S4_y_train_vaex, S4_y_test_vaex)
]


###############################################################################
# 3. Helper Functions
###############################################################################
def preprocess_vaex_data(X_train_vaex, X_test_vaex):
    """
    Convert Vaex DataFrames to numeric, replace inf/nan (out-of-core if needed),
    then apply scaling. Because scikit-learn scaling requires in-memory arrays,
    we'll convert to NumPy at the last step.
    """
    # 1) Basic cleanup in Vaex (replace NaN/Inf).
    #    Vaex uses expressions; we can do in-place replacements with fillna().

    # For example: for each column in X_train_vaex, fill missing with 0
    for col in X_train_vaex.column_names:
        # skip virtual or string columns if any
        if X_train_vaex[col].dtype.kind not in ['f','i','u']:
            continue
        X_train_vaex = X_train_vaex.fillna(col, value=0.0)
        X_train_vaex = X_train_vaex.replace([np.inf, -np.inf], value=0.0)

    for col in X_test_vaex.column_names:
        if X_test_vaex[col].dtype.kind not in ['f','i','u']:
            continue
        X_test_vaex = X_test_vaex.fillna(col, value=0.0)
        X_test_vaex = X_test_vaex.replace([np.inf, -np.inf], value=0.0)

    # 2) Convert to NumPy arrays
    X_train_np = X_train_vaex.to_pandas_df().apply(pd.to_numeric, errors='coerce').fillna(0.0).values
    X_test_np = X_test_vaex.to_pandas_df().apply(pd.to_numeric, errors='coerce').fillna(0.0).values

    # 3) Scale
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train_np)
    X_test_scaled = scaler.transform(X_test_np)

    return X_train_scaled, X_test_scaled

def process_in_batches(model, X_test, batch_size=1000):
    """
    Process test data in batches to measure memory usage and execution time.
    (Same as the pandas-based function)
    """
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
                batch_mem_usage = np.nan
        except Exception:
            batch_mem_usage = np.nan
        total_memory += batch_mem_usage

        # Measure execution time for this batch
        start_time = time.time()
        preds = model.predict(batch)
        total_time += time.time() - start_time

        predictions.extend(preds)
        del batch
        gc.collect()

    return np.array(predictions), total_memory, total_time


###############################################################################
# 4. Main Execution and Model Training (with Vaex + in-memory steps)
###############################################################################
aggregate_metrics = {
    'Precision': [],
    'Recall': [],
    'F1 Score': [],
    'Accuracy': [],
    'Memory (MB) 50k batch': [],
    'Execution Time (s) 50k batch': [],
    'Memory (MB) per read mapping': [],
    'Execution Time (s) per read mapping': [],
    'Mean CV F1 Score': []
}

best_model = None  # Will store the last 'best_model' from the loop

for (X_train_vaex, X_test_vaex, y_train_vaex, y_test_vaex) in datasets_vaex:
    # -----------------------
    # Preprocess data
    # -----------------------
    # Convert y to NumPy
    y_train = y_train_vaex.to_pandas_df().iloc[:, 0].values.ravel()
    y_test = y_test_vaex.to_pandas_df().iloc[:, 0].values.ravel()

    # Convert X and apply scaling
    X_train, X_test = preprocess_vaex_data(X_train_vaex, X_test_vaex)

    # Handle class imbalance with SMOTE (in-memory)
    min_samples = min([sum(y_train == cls) for cls in np.unique(y_train)])
    if min_samples > 1:
        k_neighbors = min(5, min_samples - 1)
        smote = SMOTE(k_neighbors=k_neighbors)
        X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
    else:
        # If a class has only one sample, duplicate it
        class_counts = pd.Series(y_train).value_counts()
        for cls in class_counts[class_counts == 1].index:
            X_train = np.vstack([X_train, X_train[y_train == cls]])
            y_train = np.hstack([y_train, y_train[y_train == cls]])
        X_train_resampled, y_train_resampled = X_train, y_train

    X_train_resampled = X_train_resampled.astype(np.float64)
    y_train_resampled = y_train_resampled.astype(np.float64)

    # -----------------------
    # Cross-validation
    # -----------------------
    sss = StratifiedShuffleSplit(n_splits=5, test_size=0.2)
    if len(np.unique(y_train_resampled)) > int(0.2 * len(y_train_resampled)):
        sss = StratifiedShuffleSplit(n_splits=5,
                                     test_size=len(np.unique(y_train_resampled)))

    param_grid = {
        'var_smoothing': np.logspace(0, -9, num=100)
    }

    model = MyGaussianNB()

    try:
        grid_search = GridSearchCV(model, param_grid, scoring='f1_weighted', cv=sss)
        grid_search.fit(X_train_resampled, y_train_resampled)
        best_model = grid_search.best_estimator_
    except Exception as e:
        print(f"Grid search failed: {e}")
        best_model = model

    # Manual CV to compute F1 scores
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

    # Refit on entire resampled data
    best_model.fit(X_train_resampled, y_train_resampled)

    # -----------------------
    # Batch prediction
    # -----------------------
    predicted, mem_usage, execution_time = process_in_batches(best_model, X_test, batch_size=1000)

    # -----------------------
    # Compute final metrics
    # -----------------------
    precision, recall, fscore, support = precision_recall_fscore_support(y_test, predicted, average='weighted')
    accuracy = accuracy_score(y_test, predicted)

    execution_time_per_read = execution_time / len(X_test)
    memory_usage_per_read = mem_usage / len(X_test)

    aggregate_metrics['Precision'].append(precision)
    aggregate_metrics['Recall'].append(recall)
    aggregate_metrics['F1 Score'].append(fscore)
    aggregate_metrics['Accuracy'].append(accuracy)
    aggregate_metrics['Mean CV F1 Score'].append(np.mean(cv_scores))
    aggregate_metrics['Memory (MB) 50k batch'].append(mem_usage)
    aggregate_metrics['Execution Time (s) 50k batch'].append(execution_time)
    aggregate_metrics['Memory (MB) per read mapping'].append(memory_usage_per_read)
    aggregate_metrics['Execution Time (s) per read mapping'].append(execution_time_per_read)

    # Cleanup in loop
    del X_train, X_test, y_train, y_test, X_train_vaex, X_test_vaex, y_train_vaex, y_test_vaex
    gc.collect()

###############################################################################
# 5. Saving Model and Creating SHAP Explanations
###############################################################################
filename = 'models/model_weights/gaussian_nb_native_vaex.joblib'
joblib.dump(best_model, filename)

# For SHAP, we do an approximate KernelExplainer
# Using a small subset for speed
explainer = shap.KernelExplainer(best_model.predict, X_train_resampled[:100])
shap_values = explainer.shap_values(X_train_resampled[:100])

feature_names = [f"Feature {i}" for i in range(X_train_resampled.shape[1])]
shap.summary_plot(shap_values,
                  X_train_resampled[:100],
                  feature_names=feature_names,
                  show=False)
plt.savefig("results/gaussian_nb_vaex_shap_values.png", dpi=300, bbox_inches='tight')
plt.close()

NaiveBayes_shap_values = pd.DataFrame(shap_values)
NaiveBayes_shap_values['Classifier'] = 'NaiveBayesClassifier_Vaex'

NaiveBayes_results = pd.DataFrame(aggregate_metrics)
NaiveBayes_results['Model'] = 'NaiveBayesClassifier_Vaex'
NaiveBayes_results['Sample'] = ['S1', 'S2', 'S3', 'S4']
NaiveBayes_results = NaiveBayes_results[[
    'Model',
    'Sample',
    'Mean CV F1 Score',
    'Precision',
    'Recall',
    'F1 Score',
    'Memory (MB) per read mapping',
    'Execution Time (s) per read mapping',
    'Memory (MB) 50k batch',
    'Execution Time (s) 50k batch'
]]

NaiveBayes_results.to_csv("results/NaiveBayes_vaex_results.csv", index=False)
NaiveBayes_shap_values.to_csv("results/NaiveBayes_vaex_shap_values.csv", index=False)

# Final cleanup
del best_model, X_train_resampled, y_train_resampled
gc.collect()
