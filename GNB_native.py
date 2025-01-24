import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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

# We need these imports to replicate scikit-learn's interface
from sklearn.base import BaseEstimator, ClassifierMixin

# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

###############################################################################
# 1. Implement a Custom Gaussian Naive Bayes Class (without sklearn.naive_bayes)
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

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors.
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        self : object
            Returns the instance itself.
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

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        Returns
        -------
        C : ndarray of shape (n_samples,)
            Predicted target values for X.
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
        Return parameters for this estimator.
        Required for compatibility with GridSearchCV.
        """
        return {"var_smoothing": self.var_smoothing}

    def set_params(self, **params):
        """
        Set the parameters of this estimator.
        Required for compatibility with GridSearchCV.
        """
        for key, value in params.items():
            setattr(self, key, value)
        return self

###############################################################################
# 2. Data Loading
###############################################################################
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

###############################################################################
# 3. Helper Functions
###############################################################################
def preprocess_data(X_train, X_test):
    """
    Convert all values to numerical, replace inf/nan, then scale.
    """
    # Convert to numeric
    X_train = pd.DataFrame(X_train).apply(pd.to_numeric, errors='coerce')
    X_test = pd.DataFrame(X_test).apply(pd.to_numeric, errors='coerce')

    # Replace NaNs, infinities, etc.
    X_train = X_train.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    X_test = X_test.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # Scale
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    return X_train, X_test

def process_in_batches(model, X_test, batch_size=1000):
    """
    Process test data in batches to measure memory usage and execution time.
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
                batch_mem_usage = np.nan  # Ensure no zero values
        except Exception:
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

###############################################################################
# 4. Main Execution and Model Training
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

best_model = None  # Will store the last 'best_model' from the loop

for X_train, X_test, y_train, y_test in datasets:
    # Preprocess data
    X_train, X_test = preprocess_data(X_train, X_test)
    X_train = X_train.astype(np.float64)
    X_test = X_test.astype(np.float64)

    # Convert y to 1D
    y_train = y_train.squeeze()
    y_test = y_test.squeeze()

    # Handle class imbalance with SMOTE
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

    X_train_resampled = np.array(X_train_resampled)
    y_train_resampled = np.array(y_train_resampled)

    # Cross-validation
    sss = StratifiedShuffleSplit(n_splits=5, test_size=0.2)

    # Adjust test_size if necessary (this part is somewhat ad-hoc)
    if len(np.unique(y_train_resampled)) > int(0.2 * len(y_train_resampled)):
        sss = StratifiedShuffleSplit(n_splits=5, test_size=len(np.unique(y_train_resampled)))

    # Hyperparameter search
    param_grid = {
        'var_smoothing': np.logspace(0, -9, num=100)
    }

    # Instantiate our custom MyGaussianNB
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

    # Batch prediction
    predicted, mem_usage, execution_time = process_in_batches(best_model, X_test, batch_size=1000)

    # Metrics
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

###############################################################################
# 5. Saving Model and Creating SHAP Explanations
###############################################################################
filename = 'models/model_weights/gaussian_nb_native_model.joblib'
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
plt.savefig("results/gaussian_nb_native_shap_values.png", dpi=300, bbox_inches='tight')
plt.close()

NaiveBayes_shap_values = pd.DataFrame(shap_values)
NaiveBayes_shap_values['Classifier'] = 'NaiveBayesClassifier'

NaiveBayes_results = pd.DataFrame(aggregate_metrics)
NaiveBayes_results['Model'] = 'NaiveBayesClassifier'
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

NaiveBayes_results.to_csv("results/NaiveBayes_native_results.csv", index=False)
NaiveBayes_shap_values.to_csv("results/NaiveBayes_native_shap_values.csv", index=False)

# Cleanup
del best_model, X_test, X_train, y_train, y_test
gc.collect()
