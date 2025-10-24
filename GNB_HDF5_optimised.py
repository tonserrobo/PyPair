import vaex
import numpy as np
import pandas as pd
import gc
import time
import warnings
from sklearn.metrics import confusion_matrix

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

###############################################################################
# 1. Define a Vaex-based Gaussian Naive Bayes
###############################################################################
class VaexGaussianNB:
    """
    Minimal Gaussian Naive Bayes implementation for Vaex DataFrames.
    - Fits by grouping on the target and computing count, mean, var for each feature.
    - Predicts chunk-by-chunk and accumulates confusion matrix.
    """

    def __init__(self, feature_columns, target_column, var_smoothing=1e-9, chunk_size=100_000):
        """
        Parameters
        ----------
        feature_columns : list of str
            Names of the numeric feature columns in the Vaex DataFrame.
        target_column : str
            Name of the target column in the Vaex DataFrame.
        var_smoothing : float, optional
            Smoothing factor to avoid zero variances.
        chunk_size : int, optional
            Number of rows to process per chunk in prediction.
        """
        self.feature_columns = feature_columns
        self.target_column = target_column
        self.var_smoothing = var_smoothing
        self.chunk_size = chunk_size

        # Learned parameters
        self.classes_ = None
        self.class_priors_ = None
        self.theta_ = None   # means
        self.sigma_ = None   # variances
        self._class_to_idx = {}
        self.fitted_ = False

    def fit(self, df: vaex.DataFrame):
        """
        Fit the Naive Bayes model using Vaex groupby to get count, mean, var.

        Parameters
        ----------
        df : vaex.DataFrame
            Must contain all features + the target column.
        """
        # 1) Group by the target to compute count, mean, var for each feature
        agg_list = [vaex.agg.count('*')]
        for feat in self.feature_columns:
            agg_list.append(vaex.agg.mean(feat))
            agg_list.append(vaex.agg.var(feat))

        grouped = df.groupby(by=self.target_column, agg=agg_list)

        # 2) Extract classes & counts
        self.classes_ = grouped[self.target_column].values
        counts = grouped['count_*'].values
        total_count = counts.sum()
        self.class_priors_ = counts / total_count

        # 3) Build arrays for means & variances
        n_classes = len(self.classes_)
        n_features = len(self.feature_columns)
        means = np.zeros((n_classes, n_features), dtype=np.float64)
        variances = np.zeros((n_classes, n_features), dtype=np.float64)

        # Because we appended mean, var in order, parse them consistently
        for feat_index, feat in enumerate(self.feature_columns):
            mean_col = f'mean_{feat}'
            var_col = f'var_{feat}'
            means[:, feat_index] = grouped[mean_col].values
            variances[:, feat_index] = grouped[var_col].values

        # 4) var_smoothing
        global_max_var = variances.max()
        epsilon = self.var_smoothing * global_max_var
        variances += epsilon

        # 5) Store learned parameters
        self.theta_ = means
        self.sigma_ = variances
        self._class_to_idx = {cls: i for i, cls in enumerate(self.classes_)}
        self.fitted_ = True

    def predict_chunk(self, X_chunk: np.ndarray) -> np.ndarray:
        """
        Predict labels for a chunk of features in memory (NumPy array).

        Parameters
        ----------
        X_chunk : np.ndarray of shape (chunk_size, n_features)

        Returns
        -------
        preds : np.ndarray of shape (chunk_size,)
            Predicted class labels.
        """
        if not self.fitted_:
            raise RuntimeError("Model has not been fitted.")

        n_samples, n_features = X_chunk.shape
        n_classes = len(self.classes_)

        # Precompute log priors
        log_prior = np.log(self.class_priors_)

        # For convenience
        mean = self.theta_
        var = self.sigma_

        # We'll build [n_samples x n_classes] log-likelihood matrix
        log_likelihood = np.zeros((n_samples, n_classes), dtype=np.float64)

        # Precompute the Gaussian constant term: sum(log(2*pi*var_c)) per class
        gauss_const = 0.5 * np.sum(np.log(2.0 * np.pi * var), axis=1)  # shape (n_classes,)

        for c in range(n_classes):
            # For each class c, compute:
            # - 1/2 sum(log(2 * pi * var_c)) [constant part]
            # - sum((x_i - mu_c)^2 / (2*var_c)) [sample part]
            diff = X_chunk - mean[c]
            exponent_term = 0.5 * np.sum((diff**2) / var[c], axis=1)
            log_likelihood[:, c] = log_prior[c] - gauss_const[c] - exponent_term

        best_class_idx = np.argmax(log_likelihood, axis=1)
        return self.classes_[best_class_idx]

    def evaluate_chunked(self, df: vaex.DataFrame) -> dict:
        """
        Predict on the entire Vaex DataFrame in chunks, accumulate a confusion matrix,
        and compute overall accuracy, precision, recall, F1 (micro-averaged).

        Parameters
        ----------
        df : vaex.DataFrame
            Must contain the features + the target column.

        Returns
        -------
        metrics : dict
            Dictionary with 'accuracy', 'precision', 'recall', 'f1' (micro).
        """
        if not self.fitted_:
            raise RuntimeError("Model has not been fitted.")

        # Build a confusion matrix
        n_classes = len(self.classes_)
        cm = np.zeros((n_classes, n_classes), dtype=np.int64)
        label_to_idx = {cls: i for i, cls in enumerate(self.classes_)}

        # We'll read the data in chunk_size blocks
        columns_needed = self.feature_columns + [self.target_column]

        for i1, i2, chunk in df.to_items(columns=columns_needed, chunk_size=self.chunk_size):
            # chunk is a dict: {col_name: np_array, ...}
            X_chunk = np.column_stack([chunk[col] for col in self.feature_columns])
            y_chunk = chunk[self.target_column]

            preds = self.predict_chunk(X_chunk)

            # Update confusion matrix
            for truth, pred in zip(y_chunk, preds):
                if truth in label_to_idx and pred in label_to_idx:
                    cm[label_to_idx[truth], label_to_idx[pred]] += 1

            # Clean up
            del X_chunk, y_chunk, chunk, preds
            gc.collect()

        return self._metrics_from_confusion_matrix(cm)

    def _metrics_from_confusion_matrix(self, cm: np.ndarray) -> dict:
        """
        Compute micro-averaged precision, recall, f1, and accuracy from confusion matrix.
        """
        total_samples = cm.sum()
        correct = np.trace(cm)
        accuracy = correct / total_samples if total_samples else 0.0

        # Micro-averaged P, R, F1 in multi-class => all become 'accuracy'
        precision = accuracy
        recall = accuracy
        if precision + recall > 0:
            f1 = 2 * precision * recall / (precision + recall)
        else:
            f1 = 0.0

        return {
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1': f1
        }


###############################################################################
# 2. Load the SAME data as the Pandas script, but with Vaex
###############################################################################
folder_path = 'data/reference_samples/'

S1_X_train_vaex = vaex.from_csv(folder_path + 'S1_extended_encoded_X_train.csv', convert=False)
S1_X_test_vaex  = vaex.from_csv(folder_path + 'S1_extended_encoded_X_test.csv', convert=False)
S1_y_train_vaex = vaex.from_csv(folder_path + 'S1_extended_encoded_y_train.csv', convert=False)
S1_y_test_vaex  = vaex.from_csv(folder_path + 'S1_extended_encoded_y_test.csv', convert=False)

S2_X_train_vaex = vaex.from_csv(folder_path + 'S2_extended_encoded_X_train.csv', convert=False)
S2_X_test_vaex  = vaex.from_csv(folder_path + 'S2_extended_encoded_X_test.csv', convert=False)
S2_y_train_vaex = vaex.from_csv(folder_path + 'S2_extended_encoded_y_train.csv', convert=False)
S2_y_test_vaex  = vaex.from_csv(folder_path + 'S2_extended_encoded_y_test.csv', convert=False)

S3_X_train_vaex = vaex.from_csv(folder_path + 'S3_extended_encoded_X_train.csv', convert=False)
S3_X_test_vaex  = vaex.from_csv(folder_path + 'S3_extended_encoded_X_test.csv', convert=False)
S3_y_train_vaex = vaex.from_csv(folder_path + 'S3_extended_encoded_y_train.csv', convert=False)
S3_y_test_vaex  = vaex.from_csv(folder_path + 'S3_extended_encoded_y_test.csv', convert=False)

S4_X_train_vaex = vaex.from_csv(folder_path + 'S4_extended_encoded_X_train.csv', convert=False)
S4_X_test_vaex  = vaex.from_csv(folder_path + 'S4_extended_encoded_X_test.csv', convert=False)
S4_y_train_vaex = vaex.from_csv(folder_path + 'S4_extended_encoded_y_train.csv', convert=False)
S4_y_test_vaex  = vaex.from_csv(folder_path + 'S4_extended_encoded_y_test.csv', convert=False)


###############################################################################
# 3. Helper: Combine X + y into a single Vaex DataFrame
###############################################################################
def combine_features_and_target(X_vaex, y_vaex, target_name='target'):
    """
    Join X and y into a single Vaex DataFrame. We assume:
    - X and y have the same row count and corresponding order.
    - y has 1 column (the target).
    """
    # Convert y to a DataFrame with a known column name
    y_column_name = y_vaex.get_column_names()[0]

    # For convenience, rename that y column to 'target' (or a user-specified name)
    y_vaex_renamed = y_vaex.rename(columns={y_column_name: target_name})

    # Concatenate along columns => they must have the same row count
    # vaex DataFrames do not have a direct "add_column from another DF" if row order is guaranteed.
    # The simplest approach is: create an in-memory array from y and assign it to X.
    y_array = y_vaex_renamed[target_name].values  # in-memory
    X_vaex = X_vaex.copy()  # to avoid mutating original
    X_vaex[target_name] = y_array  # now X_vaex includes the target as a new column

    return X_vaex


###############################################################################
# 4. Training & Evaluation on Each Dataset
#    (No SMOTE or CV; we do out-of-core training/eval in Vaex)
###############################################################################
# We'll store metrics in a dictionary of lists, similar to your original approach
aggregate_metrics = {
    'Sample': [],
    'Accuracy': [],
    'Precision': [],
    'Recall': [],
    'F1': []
}

# We loop over S1..S4 in the same style
for sample_name, (X_train_v, X_test_v, y_train_v, y_test_v) in {
    'S1': (S1_X_train_vaex, S1_X_test_vaex, S1_y_train_vaex, S1_y_test_vaex),
    'S2': (S2_X_train_vaex, S2_X_test_vaex, S2_y_train_vaex, S2_y_test_vaex),
    'S3': (S3_X_train_vaex, S3_X_test_vaex, S3_y_train_vaex, S3_y_test_vaex),
    'S4': (S4_X_train_vaex, S4_X_test_vaex, S4_y_train_vaex, S4_y_test_vaex),
}.items():
    print(f"\n=== Processing {sample_name} ===")

    # 1) Combine X + y into single DF
    #    We'll call the target column "label" for clarity
    train_df = combine_features_and_target(X_train_v, y_train_v, target_name='label')
    test_df  = combine_features_and_target(X_test_v,  y_test_v,  target_name='label')

    # 2) Basic pre-processing in Vaex (similar to fillna, etc.)
    #    For large datasets, do out-of-core transformations
    for df_ in [train_df, test_df]:
        # e.g., replace NaNs with 0 in each numeric column
        for col in df_.column_names:
            # If it's numeric, fill NaNs
            if df_[col].dtype.kind in ['f','i','u']:
                df_ = df_.fillna(col, value=0.0)
                df_ = df_.replace([np.inf, -np.inf], value=0.0)

    # 3) Identify the feature columns
    #    (Assume everything except 'label' is a feature)
    feature_columns = [col for col in train_df.column_names if col != 'label']

    # 4) Fit the Vaex-based GNB
    model = VaexGaussianNB(feature_columns=feature_columns, target_column='label',
                           var_smoothing=1e-9, chunk_size=50_000)
    model.fit(train_df)

    # 5) Evaluate chunk-by-chunk on the test set
    metrics = model.evaluate_chunked(test_df)
    print(f"Metrics for {sample_name}: {metrics}")

    # 6) Store results
    aggregate_metrics['Sample'].append(sample_name)
    aggregate_metrics['Accuracy'].append(metrics['accuracy'])
    aggregate_metrics['Precision'].append(metrics['precision'])
    aggregate_metrics['Recall'].append(metrics['recall'])
    aggregate_metrics['F1'].append(metrics['f1'])

    # Cleanup references
    del train_df, test_df
    gc.collect()

###############################################################################
# 5. Output aggregated metrics
###############################################################################
result_df = pd.DataFrame(aggregate_metrics)
print("\n=== Final Aggregated Metrics ===")
print(result_df)

# Save to CSV if desired
result_df.to_csv("results/NaiveBayes_vaex_results.csv", index=False)
