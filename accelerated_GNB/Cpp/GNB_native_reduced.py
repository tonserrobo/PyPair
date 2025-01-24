import numpy as np
import pandas as pd
import warnings
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
import joblib

# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

###############################################################################
# 1. Custom Gaussian Naive Bayes (without sklearn.naive_bayes)
###############################################################################
class MyGaussianNB:
    """
    A custom implementation of Gaussian Naive Bayes.
    """
    def __init__(self, var_smoothing=1e-9):
        self.var_smoothing = var_smoothing

    def fit(self, X, y):
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y)

        self.classes_, y_indices = np.unique(y, return_inverse=True)
        n_classes = len(self.classes_)
        n_features = X.shape[1]

        # Mean and variance per class
        self.theta_ = np.zeros((n_classes, n_features), dtype=np.float64)
        self.sigma_ = np.zeros((n_classes, n_features), dtype=np.float64)
        self.class_count_ = np.zeros(n_classes, dtype=np.float64)

        for idx, c in enumerate(self.classes_):
            X_c = X[y_indices == idx]
            self.class_count_[idx] = X_c.shape[0]
            self.theta_[idx, :] = np.mean(X_c, axis=0)
            self.sigma_[idx, :] = np.var(X_c, axis=0)

        self.class_prior_ = self.class_count_ / self.class_count_.sum()

        # Stabilize variances
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
            # Constant term per feature
            log_prob_const = -0.5 * np.sum(np.log(2.0 * np.pi * var[c]))
            # Variable term per sample
            log_prob_x = -0.5 * np.sum(((X - mean[c]) ** 2) / var[c], axis=1)
            log_likelihood[:, c] = log_prior[c] + log_prob_const + log_prob_x

        return log_likelihood

###############################################################################
# 2. Data Loading (S1 example)
###############################################################################
folder_path = 'data/reference_samples/'
S1_X_train = pd.read_csv(folder_path + 'S1_extended_encoded_X_train.csv')
S1_X_test = pd.read_csv(folder_path + 'S1_extended_encoded_X_test.csv')
S1_y_train = pd.read_csv(folder_path + 'S1_extended_encoded_y_train.csv')
S1_y_test = pd.read_csv(folder_path + 'S1_extended_encoded_y_test.csv')

###############################################################################
# 3. Helper Functions
###############################################################################
def preprocess_data(X_train, X_test):
    """
    Convert values to numeric, replace inf/nan, then scale.
    """
    X_train = pd.DataFrame(X_train).apply(pd.to_numeric, errors='coerce')
    X_test = pd.DataFrame(X_test).apply(pd.to_numeric, errors='coerce')

    X_train = X_train.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    X_test = X_test.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    return X_train_scaled, X_test_scaled

def batch_predict(model, X_test, batch_size=1000):
    """
    Perform prediction on test data in batches.
    """
    predictions = []
    for start in range(0, len(X_test), batch_size):
        end = start + batch_size
        batch = X_test[start:end]
        preds = model.predict(batch)
        predictions.extend(preds)
    return np.array(predictions)

###############################################################################
# 4. Main Execution
###############################################################################
# Preprocess
X_train, X_test = preprocess_data(S1_X_train, S1_X_test)
y_train = S1_y_train.squeeze()
y_test = S1_y_test.squeeze()

# Handle class imbalance with SMOTE (optional)
min_samples = min([sum(y_train == cls) for cls in np.unique(y_train)])
if min_samples > 1:
    k_neighbors = min(5, min_samples - 1)
    smote = SMOTE(k_neighbors=k_neighbors)
    X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
else:
    # If a class has only one sample, replicate it
    import numpy as np
    class_counts = pd.Series(y_train).value_counts()
    for cls in class_counts[class_counts == 1].index:
        X_train = np.vstack([X_train, X_train[y_train == cls]])
        y_train = np.hstack([y_train, y_train[y_train == cls]])
    X_train_resampled, y_train_resampled = X_train, y_train

# Train model
model = MyGaussianNB()
model.fit(X_train_resampled, y_train_resampled)

# Predict on test set in batches
y_pred = batch_predict(model, X_test, batch_size=1000)

# Save the trained model
model_filename = 'models/model_weights/gaussian_nb_S1.joblib'
joblib.dump(model, model_filename)

print(f"Done. Model saved to: {model_filename}")
