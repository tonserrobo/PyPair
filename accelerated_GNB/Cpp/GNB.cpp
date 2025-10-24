#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <map>

/**
 * Utility function to read a CSV file into a 2D std::vector<double>.
 * - Assumes no header
 * - Splits on commas
 * - Ignores empty lines
 */
std::vector<std::vector<double>> readCSV(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<double> row;
        std::stringstream ss(line);
        std::string value;
        while (std::getline(ss, value, ',')) {
            try {
                row.push_back(std::stod(value));
            } catch (...) {
                row.push_back(0.0);
            }
        }
        data.push_back(row);
    }
    file.close();
    return data;
}

/**
 * StandardScaler-like function:
 * - Computes mean and stddev per column in X
 * - Scales each column to (x - mean) / stddev
 */
void standardScale(std::vector<std::vector<double>>& X) {
    if (X.empty()) return;

    size_t n_samples = X.size();
    size_t n_features = X[0].size();

    // Compute mean, stddev per feature
    std::vector<double> means(n_features, 0.0);
    std::vector<double> stds(n_features, 0.0);

    // Sum up each column
    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < n_features; ++j) {
            means[j] += X[i][j];
        }
    }
    for (size_t j = 0; j < n_features; ++j) {
        means[j] /= static_cast<double>(n_samples);
    }

    // Compute variance -> stddev
    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < n_features; ++j) {
            double diff = X[i][j] - means[j];
            stds[j] += diff * diff;
        }
    }
    for (size_t j = 0; j < n_features; ++j) {
        stds[j] = std::sqrt(stds[j] / static_cast<double>(n_samples));
        if (stds[j] == 0.0) {
            stds[j] = 1e-9; // prevent divide by zero
        }
    }

    // Scale in place
    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < n_features; ++j) {
            X[i][j] = (X[i][j] - means[j]) / stds[j];
        }
    }
}

/**
 * Minimal “SMOTE-like” oversampling for classes that have very few samples.
 * - Full SMOTE is more advanced (random interpolation), but here we simply replicate
 *   minority classes if they have less than "min_required" samples.
 *
 *   y: class vector, same length as X
 *   X: NxM matrix (N samples, M features)
 */
void oversampleRareClasses(std::vector<std::vector<double>>& X,
                           std::vector<double>& y,
                           size_t min_required = 2)
{
    // Count how many times each class appears
    std::map<double, int> class_counts;
    for (double label : y) {
        class_counts[label]++;
    }

    // For each class that has < min_required, replicate those rows
    for (auto& kv : class_counts) {
        double cls = kv.first;
        int count = kv.second;
        if (count < static_cast<int>(min_required)) {
            // replicate enough times
            size_t needed = min_required - count;
            for (size_t i = 0; i < y.size(); ++i) {
                if (y[i] == cls) {
                    // replicate row
                    for (size_t rep = 0; rep < needed; ++rep) {
                        X.push_back(X[i]);
                        y.push_back(cls);
                    }
                }
            }
        }
    }
}

/**
 * A custom Gaussian Naive Bayes class (similar to the Python MyGaussianNB).
 */
class MyGaussianNB {
public:
    explicit MyGaussianNB(double var_smoothing=1e-9)
        : var_smoothing_(var_smoothing) {}

    /**
     * Fit the model with training data X, y.
     */
    void fit(const std::vector<std::vector<double>>& X,
             const std::vector<double>& y)
    {
        // Identify unique classes
        classes_.clear();
        {
            std::map<double, bool> unique_map;
            for (auto val : y) {
                unique_map[val] = true;
            }
            for (auto& kv : unique_map) {
                classes_.push_back(kv.first);
            }
            // Sort classes so the order is deterministic
            std::sort(classes_.begin(), classes_.end());
        }
        size_t n_classes = classes_.size();
        size_t n_samples = X.size();
        if (n_samples == 0) return;
        size_t n_features = X[0].size();

        // Create storage
        theta_.assign(n_classes, std::vector<double>(n_features, 0.0));
        sigma_.assign(n_classes, std::vector<double>(n_features, 0.0));
        class_count_.assign(n_classes, 0.0);

        // For quick indexing of classes
        std::map<double, int> class_index;
        for (size_t i = 0; i < n_classes; ++i) {
            class_index[classes_[i]] = static_cast<int>(i);
        }

        // Sum up for mean and variance
        for (size_t i = 0; i < n_samples; ++i) {
            double label = y[i];
            int idx = class_index[label];
            class_count_[idx] += 1.0;
            for (size_t j = 0; j < n_features; ++j) {
                theta_[idx][j] += X[i][j];
            }
        }
        // Compute means
        for (size_t c = 0; c < n_classes; ++c) {
            if (class_count_[c] > 0) {
                for (size_t j = 0; j < n_features; ++j) {
                    theta_[c][j] /= class_count_[c];
                }
            }
        }
        // Compute variance
        for (size_t i = 0; i < n_samples; ++i) {
            double label = y[i];
            int idx = class_index[label];
            for (size_t j = 0; j < n_features; ++j) {
                double diff = X[i][j] - theta_[idx][j];
                sigma_[idx][j] += diff * diff;
            }
        }
        // Finalize variance
        for (size_t c = 0; c < n_classes; ++c) {
            if (class_count_[c] > 1) {
                for (size_t j = 0; j < n_features; ++j) {
                    sigma_[c][j] /= class_count_[c];
                }
            }
        }
        // Class priors
        class_prior_.assign(n_classes, 0.0);
        double total_count = std::accumulate(class_count_.begin(), class_count_.end(), 0.0);
        for (size_t c = 0; c < n_classes; ++c) {
            class_prior_[c] = class_count_[c] / total_count;
        }

        // Add var_smoothing
        double global_var = 0.0;
        // compute overall variance across entire dataset
        {
            // We'll do an approximate approach: sum up sigma and pick max
            double max_var = 0.0;
            for (size_t c = 0; c < n_classes; ++c) {
                for (size_t j = 0; j < n_features; ++j) {
                    if (sigma_[c][j] > max_var) {
                        max_var = sigma_[c][j];
                    }
                }
            }
            global_var = max_var;
        }
        double epsilon = var_smoothing_ * global_var;
        for (size_t c = 0; c < n_classes; ++c) {
            for (size_t j = 0; j < n_features; ++j) {
                sigma_[c][j] += epsilon;
            }
        }
    }

    /**
     * Predict class labels for X.
     */
    std::vector<double> predict(const std::vector<std::vector<double>>& X) const {
        std::vector<double> preds;
        preds.reserve(X.size());
        for (const auto& row : X) {
            preds.push_back(predictSingle(row));
        }
        return preds;
    }

private:
    double var_smoothing_;
    std::vector<double> classes_;                    // Unique class labels
    std::vector<double> class_prior_;                // Prior probabilities per class
    std::vector<std::vector<double>> theta_;         // Means per class (c, feature)
    std::vector<std::vector<double>> sigma_;         // Variances per class (c, feature)
    std::vector<double> class_count_;                // Number of samples per class

    /**
     * Predict single sample
     */
    double predictSingle(const std::vector<double>& x) const {
        // Compute log-likelihood for each class
        size_t n_classes = classes_.size();
        if (n_classes == 0) {
            return 0.0; // fallback
        }

        std::vector<double> log_probs(n_classes, 0.0);

        for (size_t c = 0; c < n_classes; ++c) {
            // log of prior
            double lp = std::log(class_prior_[c]);
            // sum of gaussian log-likelihood
            double sum_log = 0.0;
            for (size_t j = 0; j < theta_[c].size(); ++j) {
                double diff = x[j] - theta_[c][j];
                double var_cj = sigma_[c][j];
                // -0.5 * (log(2*pi*var) + diff^2 / var)
                sum_log += -0.5 * (std::log(2.0 * M_PI * var_cj) + (diff * diff) / var_cj);
            }
            log_probs[c] = lp + sum_log;
        }

        // Argmax
        auto max_iter = std::max_element(log_probs.begin(), log_probs.end());
        size_t idx = static_cast<size_t>(std::distance(log_probs.begin(), max_iter));
        return classes_[idx];
    }
};

/**
 * Helper for batch predictions.
 */
std::vector<double> batchPredict(const MyGaussianNB& model,
                                 const std::vector<std::vector<double>>& X,
                                 size_t batch_size = 1000)
{
    std::vector<double> all_preds;
    all_preds.reserve(X.size());

    for (size_t start = 0; start < X.size(); start += batch_size) {
        size_t end = std::min(start + batch_size, X.size());
        std::vector<std::vector<double>> batch(X.begin() + start, X.begin() + end);
        std::vector<double> preds = model.predict(batch);
        all_preds.insert(all_preds.end(), preds.begin(), preds.end());
    }
    return all_preds;
}

/**
 * Example save-model function:
 * - Writes model parameters to a text file (very naive approach).
 */
void saveModel(const MyGaussianNB& model, const std::string& filepath) {
    // In a real implementation, store all relevant arrays:
    // classes_, class_prior_, theta_, sigma_, etc.
    // We'll keep this placeholder to demonstrate usage.
    std::ofstream out(filepath);
    if (out.is_open()) {
        out << "NaiveBayesModel - placeholder serialization" << std::endl;
        out.close();
    }
    std::cout << "Model saved to " << filepath << std::endl;
}

int main() {
    // 1. Load CSV data for S1
    // Adjust these filenames to match your system if needed
    std::string folder_path = "data/reference_samples/";
    auto S1_X_train = readCSV(folder_path + "S1_extended_encoded_X_train.csv");
    auto S1_X_test  = readCSV(folder_path + "S1_extended_encoded_X_test.csv");
    auto ytrain_2d  = readCSV(folder_path + "S1_extended_encoded_y_train.csv");
    auto ytest_2d   = readCSV(folder_path + "S1_extended_encoded_y_test.csv");

    // Convert y from 2D -> 1D
    std::vector<double> S1_y_train, S1_y_test;
    for (const auto& row : ytrain_2d) {
        if (!row.empty()) S1_y_train.push_back(row[0]);
    }
    for (const auto& row : ytest_2d) {
        if (!row.empty()) S1_y_test.push_back(row[0]);
    }

    // 2. Preprocess (standard scale)
    standardScale(S1_X_train);
    standardScale(S1_X_test);

    // 3. Simple SMOTE-like oversampling for rare classes
    //    (If a class has fewer than min_required samples, replicate them)
    size_t min_required = 2; // or however many you want
    oversampleRareClasses(S1_X_train, S1_y_train, min_required);

    // 4. Train MyGaussianNB
    MyGaussianNB model(1e-9);
    model.fit(S1_X_train, S1_y_train);

    // 5. Predict on test set in batches
    auto y_pred = batchPredict(model, S1_X_test, 1000);

    // 6. Save the model to disk (very simplistic)
    saveModel(model, "gaussian_nb_S1.txt");

    std::cout << "Done. Predictions size: " << y_pred.size() << std::endl;
    return 0;
}
