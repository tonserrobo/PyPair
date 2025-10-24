import matplotlib.pyplot as plt
import numpy as np
import polars as pl  # <-- Polars
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
import pandas as pd  # still used for final output (csv) or partial manipulations

# We need these imports to replicate scikit-learn's interface
from sklearn.base import BaseEstimator, ClassifierMixin

# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

###############################################################################
# 1. Implement the same Custom Gaussian Naive Bayes Class
###############################################################################
class MyGaussianNB(BaseEstimator, ClassifierMixin):
    """
    A custom implementation of Gaussian Naive Bayes, emulating scikit-learn's
    interface so it can be used with GridSearchCV, etc.
    """
    def __init__(self, var_smoothing=1e-9):
        """
        Parameters
        ----------
        var_smoothing : float
            Portion of the largest variance of all features that is added to
            variances for calculation stability.
        """
        self.var_smoothing = var_smoothing

    def fit(self, X, y):
        """
        Fit the Gaussian Naive Bayes model according to the given training data.
        """
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y)

        # Identify unique classes and their indices
        self.classes_, y_indices = np.unique(y, return_inverse=True)

        # Calculate the number of classes and features
        n_classes = len(self.classes_)
        n_features = X.shape[1]

        # Preallocate arrays for means and variances
        self.theta_ = np.zeros((n_classes, n_features), dtype=np.float64)
        self.sigma_ = np.zeros((n_classes, n_features), dtype=np.float64)
        self.class_count_ = np.zeros(n_classes, dtype=np.float64)

        # Calculate mean, variance, and class count for each class
        for idx, c in enumerate(self.classes_):
            X_c = X[y_indices == idx, :]
            self.class_count_[idx] = X_c.shape[0]
            self.theta_[idx, :] = np.mean(X_c, axis=0)
            self.sigma_[idx, :] = np.var(X_c, axis=0)

        # Calculate class prior probabilities
        self.class_prior_ = self.class_count_ / self.class_count_.sum()

        # Add var_smoothing to variances (similar to sklearn's approach)
        # Typically, var_smoothing is multiplied by the overall data variance's max
        epsilon = self.var_smoothing * np.var(X, axis=0).max()
        self.sigma_ += epsilon

        return self

    def predict(self, X):
        """
        Perform classification on an array of test vectors X.
        """
        X = np.asarray(X, dtype=np.float64)
        # Compute the posterior log probability for each class
        jll = self._joint_log_likelihood(X)
        # Choose the class with the highest log posterior
        return self.classes_[np.argmax(jll, axis=1)]

    def _joint_log_likelihood(self, X):
        """
        Calculate the joint log likelihood of each sample in X
        belonging to each of the classes.
        """
        n_samples, n_features = X.shape
        n_classes = len(self.classes_)

        # We will accumulate the log probabilities in log_likelihood
        log_likelihood = np.zeros((n_samples, n_classes), dtype=np.float64)

        # Precompute log of class priors
        log_prior = np.log(self.class_prior_)

        # For numerical stability in log of Gaussian formula
        var = self.sigma_
        mean = self.theta_

        # Compute log probabilities
        # log P(x_i | c) = - 1/2 log(2 pi sigma_c) - (x_i - mu_c)^2 / (2 sigma_c)
        for c in range(n_classes):
            # constant term per feature dimension
            log_prob_const = -0.5 * np.sum(np.log(2.0 * np.pi * var[c]))
            # variable term for each sample
            log_prob_x = -0.5 * np.sum(((X - mean[c]) ** 2) / (var[c]), axis=1)
            # total log likelihood for class c is the sum of the constant,
            # variable part, plus the log of the prior
            log_likelihood[:, c] = log_prior[c] + log_prob_const + log_prob_x

        return log_likelihood

    def get_params(self, deep=True):
        """
        Return parameters for this estimator (for GridSearchCV).
        """
        return {"var_smoothing": self.var_smoothing}

    def set_params(self, **params):
        """
        Set parameters of this estimator (for GridSearchCV).
        """
        for key, value in params.items():
            setattr(self, key, value)
        return self


###############################################################################
# 2. Data Loading with Polars
###############################################################################
folder_path = 'data/reference_samples/'

# Instead of pd.read_csv, we do pl.read_csv
# By default, Polars will infer column types.
# If your dataset is huge, consider polars.scan_csv(...) + lazy transforms.

S1_X_train = pl.read_csv(folder_path + 'S1_extended_encoded_X_train.csv', infer_schema_length=10000, ignore_errors=True)
S1_X_test  = pl.read_csv(folder_path + 'S1_extended_encoded_X_test.csv', infer_schema_length=10000, ignore_errors=True)
S1_y_train = pl.read_csv(folder_path + 'S1_extended_encoded_y_train.csv', infer_schema_length=10000, ignore_errors=True)
S1_y_test  = pl.read_csv(folder_path + 'S1_extended_encoded_y_test.csv', infer_schema_length=10000, ignore_errors=True)

S2_X_train = pl.read_csv(folder_path + 'S2_extended_encoded_X_train.csv', infer_schema_length=10000, ignore_errors=True)
S2_X_test  = pl.read_csv(folder_path + 'S2_extended_encoded_X_test.csv', infer_schema_length=10000, ignore_errors=True)
S2_y_train = pl.read_csv(folder_path + 'S2_extended_encoded_y_train.csv', infer_schema_length=10000, ignore_errors=True)
S2_y_test  = pl.read_csv(folder_path + 'S2_extended_encoded_y_test.csv', infer_schema_length=10000, ignore_errors=True)

S3_X_train = pl.read_csv(folder_path + 'S3_extended_encoded_X_train.csv', infer_schema_length=10000, ignore_errors=True)
S3_X_test  = pl.read_csv(folder_path + 'S3_extended_encoded_X_test.csv', infer_schema_length=10000, ignore_errors=True)
S3_y_train = pl.read_csv(folder_path + 'S3_extended_encoded_y_train.csv', infer_schema_length=10000, ignore_errors=True)
S3_y_test  = pl.read_csv(folder_path + 'S3_extended_encoded_y_test.csv', infer_schema_length=10000, ignore_errors=True)

S4_X_train = pl.read_csv(folder_path + 'S4_extended_encoded_X_train.csv', infer_schema_length=10000, ignore_errors=True)
S4_X_test  = pl.read_csv(folder_path + 'S4_extended_encoded_X_test.csv', infer_schema_length=10000, ignore_errors=True)
S4_y_train = pl.read_csv(folder_path + 'S4_extended_encoded_y_train.csv', infer_schema_length=10000, ignore_errors=True)
S4_y_test  = pl.read_csv(folder_path + 'S4_extended_encoded_y_test.csv', infer_schema_length=10000, ignore_errors=True)


###############################################################################
# 3. Helper Functions: Polars-based pre-processing
###############################################################################
def preprocess_data_polars(X_train_pl: pl.DataFrame, X_test_pl: pl.DataFrame):
    """
    Convert columns to numeric (if needed), replace inf/nan, and then
    convert to NumPy for scaling with scikit-learn.

    Returns (X_train_np, X_test_np) as scaled NumPy arrays.
    """

    # 1) Replace NaNs, inf, -inf with 0.0 using Polars expressions
    #    If your dataset is large, you could do this lazily or chunked.
    #    For demonstration, we do it eagerly.
    #    Polars doesn't do in-place; we reassign after each expression.


    # For each numeric column, fill nulls and replace infinite with 0
    numeric_cols = [col for col, dtype in X_train_pl.schema.items()
                    if dtype in (pl.Float32, pl.Float64, pl.Int32, pl.Int64, pl.UInt32, pl.UInt64)]



    # Fill null with 0
    X_train_pl = X_train_pl.fill_null(0.0)
    X_test_pl  = X_test_pl.fill_null(0.0)

    # Polars doesn’t have a direct "replace inf with 0" built-in yet.
    # We can do a trick: if abs(value) is very large, set to 0.
    # Or convert to a pandas df to do the replacement quickly (like below).
    X_train_pd = X_train_pl.to_pandas()
    X_test_pd  = X_test_pl.to_pandas()

    X_train_pd = X_train_pd.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    X_test_pd  = X_test_pd.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # 2) Scale using StandardScaler
    scaler = StandardScaler()
    X_train_np = scaler.fit_transform(X_train_pd)
    X_test_np  = scaler.transform(X_test_pd)

    return X_train_np, X_test_np


def process_in_batches(model, X_test: np.ndarray, batch_size=1000):
    """
    Identical function from the pandas script to measure memory usage,
    execution time, and generate predictions in chunks.
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

        # Measure execution time
        start_time = time.time()
        preds = model.predict(batch)
        total_time += time.time() - start_time

        predictions.extend(preds)
        del batch
        gc.collect()

    return np.array(predictions), total_memory, total_time


###############################################################################
# 4. Main Execution & Model Training
###############################################################################
datasets = [
    (S1_X_train, S1_X_test, S1_y_train, S1_y_test),
    (S2_X_train, S2_X_test, S2_y_train, S2_y_test),
    (S3_X_train, S3_X_test, S3_y_train, S3_y_test),
    (S4_X_train, S4_X_test, S4_y_train, S4_y_test)
]

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

best_model = None

for X_train_pl, X_test_pl, y_train_pl, y_test_pl in datasets:
    # 1) Convert Polars DataFrame -> NumPy for y
    #    We assume the label column is the single column in the Polars DF.
    y_train_np = y_train_pl.to_numpy().ravel()
    y_test_np  = y_test_pl.to_numpy().ravel()

    # 2) Preprocess Polars data for X
    X_train_np, X_test_np = preprocess_data_polars(X_train_pl, X_test_pl)
    X_train_np = X_train_np.astype(np.float64)
    X_test_np  = X_test_np.astype(np.float64)

    # 3) Handle class imbalance with SMOTE (in-memory)
    min_samples = min([sum(y_train_np == cls) for cls in np.unique(y_train_np)])
    if min_samples > 1:
        k_neighbors = min(5, min_samples - 1)
        smote = SMOTE(k_neighbors=k_neighbors)
        X_train_resampled, y_train_resampled = smote.fit_resample(X_train_np, y_train_np)
    else:
        # If a class has only one sample, duplicate it
        class_counts = pd.Series(y_train_np).value_counts()
        for cls_ in class_counts[class_counts == 1].index:
            X_train_np = np.vstack([X_train_np, X_train_np[y_train_np == cls_]])
            y_train_np = np.hstack([y_train_np, y_train_np[y_train_np == cls_]])
        X_train_resampled, y_train_resampled = X_train_np, y_train_np

    # 4) Cross-validation
    sss = StratifiedShuffleSplit(n_splits=5, test_size=0.2)
    if len(np.unique(y_train_resampled)) > int(0.2 * len(y_train_resampled)):
        sss = StratifiedShuffleSplit(n_splits=5,
                                     test_size=len(np.unique(y_train_resampled)))

    # Hyperparameter search
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

    # 5) Manual cross-validation to compute F1
    cv_scores = []
    for train_idx, test_idx in sss.split(X_train_resampled, y_train_resampled):
        X_cv_train, X_cv_test = X_train_resampled[train_idx], X_train_resampled[test_idx]
        y_cv_train, y_cv_test = y_train_resampled[train_idx], y_train_resampled[test_idx]

        best_model.fit(X_cv_train, y_cv_train)
        y_cv_pred = best_model.predict(X_cv_test)
        _, _, f1, _ = precision_recall_fscore_support(y_cv_test, y_cv_pred, average='weighted')
        cv_scores.append(f1)

    print(f"Cross-validation F1 scores: {cv_scores}")
    print(f"Mean CV F1 score: {np.mean(cv_scores)}")

    # 6) Refit on entire data
    best_model.fit(X_train_resampled, y_train_resampled)

    # 7) Batch prediction
    predicted, mem_usage, execution_time = process_in_batches(best_model, X_test_np, batch_size=1000)

    precision, recall, fscore, support = precision_recall_fscore_support(y_test_np, predicted, average='weighted')
    accuracy = accuracy_score(y_test_np, predicted)

    execution_time_per_read = execution_time / len(X_test_np)
    memory_usage_per_read = mem_usage / len(X_test_np)

    aggregate_metrics['Precision'].append(precision)
    aggregate_metrics['Recall'].append(recall)
    aggregate_metrics['F1 Score'].append(fscore)
    aggregate_metrics['Accuracy'].append(accuracy)
    aggregate_metrics['Mean CV F1 Score'].append(np.mean(cv_scores))
    aggregate_metrics['Memory (MB) 50k batch'].append(mem_usage)
    aggregate_metrics['Execution Time (s) 50k batch'].append(execution_time)
    aggregate_metrics['Memory (MB) per read mapping'].append(memory_usage_per_read)
    aggregate_metrics['Execution Time (s) per read mapping'].append(execution_time_per_read)

    # Cleanup references
    del (X_train_pl, X_test_pl, y_train_pl, y_test_pl,
         X_train_np, X_test_np, X_train_resampled, y_train_resampled)
    gc.collect()

###############################################################################
# 5. Save Model & Summarize
###############################################################################
filename = 'models/model_weights/gaussian_nb_polars_model.joblib'
joblib.dump(best_model, filename)

# # Example: SHAP analysis (sampling first 100 rows of training)
# # If the dataset is huge, you might want a smaller subset
# # The exact lines can remain the same as in your original script:
# explainer = shap.KernelExplainer(best_model.predict, X_train_resampled[:100])
# shap_values = explainer.shap_values(X_train_resampled[:100])
#
# feature_names = [f"Feature {i}" for i in range(X_train_resampled.shape[1])]
# shap.summary_plot(
#     shap_values,
#     X_train_resampled[:100],
#     feature_names=feature_names,
#     show=False
# )
# plt.savefig("results/gaussian_nb_polars_shap_values.png", dpi=300, bbox_inches='tight')
# plt.close()

# # Collect results
# NaiveBayes_shap_values = pd.DataFrame(shap_values)
# NaiveBayes_shap_values['Classifier'] = 'NaiveBayesClassifier_Polars'

NaiveBayes_results = pd.DataFrame(aggregate_metrics)
NaiveBayes_results['Model'] = 'NaiveBayesClassifier_Polars'
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

# Save results to CSV
NaiveBayes_results.to_csv("results/NaiveBayes_polars_results.csv", index=False)
# NaiveBayes_shap_values.to_csv("results/NaiveBayes_polars_shap_values.csv", index=False)

del best_model
gc.collect()
