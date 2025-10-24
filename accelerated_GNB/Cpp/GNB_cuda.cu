#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

// =========== CUDA ERROR CHECK MACRO =============
#define CUDA_CHECK(call)                                            \
    do {                                                            \
        cudaError_t err = call;                                     \
        if (err != cudaSuccess) {                                   \
            std::cerr << "CUDA error: " << cudaGetErrorString(err)  \
                      << " in " << __FILE__                         \
                      << " at line " << __LINE__ << std::endl;      \
            exit(EXIT_FAILURE);                                     \
        }                                                           \
    } while (0)

// ==============================================================================
// 1) Utility: Read CSV (no header) into a 2D std::vector<double>
// ==============================================================================
std::vector<std::vector<double>> readCSV(const std::string& filepath) {
    std::ifstream file(filepath);
    std::vector<std::vector<double>> data;

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filepath << std::endl;
        return data; // empty
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string val;
        std::vector<double> row;
        while (std::getline(ss, val, ',')) {
            try {
                row.push_back(std::stod(val));
            } catch (...) {
                // If parsing fails, treat as 0.0
                row.push_back(0.0);
            }
        }
        data.push_back(row);
    }
    file.close();
    return data;
}

// ==============================================================================
// 2) Utility: Standard Scale (column-wise). Similar to Python's StandardScaler.
// ==============================================================================
void standardScale(std::vector<std::vector<double>>& X) {
    if (X.empty()) return;

    const size_t n_samples = X.size();
    const size_t n_features = X[0].size();

    // Compute mean and stddev per column
    std::vector<double> means(n_features, 0.0);
    std::vector<double> stds(n_features, 0.0);

    // Sum up columns
    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < n_features; ++j) {
            means[j] += X[i][j];
        }
    }
    for (size_t j = 0; j < n_features; ++j) {
        means[j] /= static_cast<double>(n_samples);
    }

    // Compute variances -> stddev
    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < n_features; ++j) {
            double diff = X[i][j] - means[j];
            stds[j] += diff * diff;
        }
    }
    for (size_t j = 0; j < n_features; ++j) {
        stds[j] = std::sqrt(stds[j] / static_cast<double>(n_samples));
        if (stds[j] == 0.0) {
            stds[j] = 1e-9; // prevent division by zero
        }
    }

    // Apply scaling
    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < n_features; ++j) {
            X[i][j] = (X[i][j] - means[j]) / stds[j];
        }
    }
}

// ==============================================================================
// 3) Minimal SMOTE-like oversampling: replicate rows for any class < min_required
//    (For a real SMOTE, you would need interpolation, etc.)
// ==============================================================================
void oversampleRareClasses(std::vector<std::vector<double>>& X,
                           std::vector<double>& y,
                           size_t min_required = 2)
{
    // Count how many times each class appears
    std::unordered_map<double, size_t> class_counts;
    for (auto cls : y) {
        class_counts[cls] += 1;
    }

    // For each class that has fewer than min_required samples,
    // replicate those samples to reach min_required.
    for (auto& kv : class_counts) {
        double cls   = kv.first;
        size_t count = kv.second;
        if (count < min_required) {
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

// ==============================================================================
// 4) GPU Kernel for prediction
//    - Each thread processes exactly one sample from X
// ==============================================================================
__global__
void predictKernel(const double* d_X,          // NxM samples in row-major
                   const double* d_theta,      // (n_classes * n_features) row-major
                   const double* d_sigma,      // (n_classes * n_features) row-major
                   const double* d_classPrior, // (n_classes)
                   const double* d_classes,    // (n_classes)
                   double* d_out,              // (n_samples) -> predicted class label
                   size_t n_samples,
                   size_t n_features,
                   size_t n_classes)
{
    size_t i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n_samples) return;

    // Compute log-prob for each class c
    double max_log_prob = -1e308;
    double best_class   = 0.0;

    for (size_t c = 0; c < n_classes; ++c) {
        // start with log(prior)
        double log_prob = log(d_classPrior[c]);

        // offset for class c in theta_, sigma_
        size_t offset = c * n_features;

        // sum of Gaussian log-likelihood across features
        //   -0.5 * [ log(2*pi*var) + (x - mean)^2 / var ]
        double sum_gauss = 0.0;
        for (size_t j = 0; j < n_features; ++j) {
            double x_ij = d_X[i * n_features + j];
            double diff = x_ij - d_theta[offset + j];
            double var_cj = d_sigma[offset + j];
            sum_gauss += -0.5 * (log(2.0 * M_PI * var_cj) + (diff * diff) / var_cj);
        }
        log_prob += sum_gauss;

        // track best
        if (log_prob > max_log_prob) {
            max_log_prob = log_prob;
            best_class   = d_classes[c];
        }
    }

    d_out[i] = best_class;
}

// ==============================================================================
// 5) Naive Bayes Class
//    - CPU-based fit
//    - GPU-based predict
// ==============================================================================
class MyGaussianNB {
public:
    MyGaussianNB(double var_smoothing=1e-9)
        : var_smoothing_(var_smoothing),
          d_theta_(nullptr), d_sigma_(nullptr),
          d_class_prior_(nullptr), d_classes_(nullptr) {}

    ~MyGaussianNB() {
        freeDeviceMemory();
    }

    // Fit on CPU
    void fit(const std::vector<std::vector<double>>& X,
             const std::vector<double>& y);

    // Predict on GPU
    std::vector<double> predict(const std::vector<std::vector<double>>& X) const;

private:
    double var_smoothing_;

    // Host-side
    std::vector<double> classes_;       // unique class labels
    std::vector<double> class_prior_;   // length = n_classes
    std::vector<double> theta_;         // length = n_classes * n_features
    std::vector<double> sigma_;         // length = n_classes * n_features
    size_t n_classes_  = 0;
    size_t n_features_ = 0;

    // Device-side copies
    double* d_theta_;
    double* d_sigma_;
    double* d_class_prior_;
    double* d_classes_;

    void copyToDevice();
    void freeDeviceMemory();
};

void MyGaussianNB::fit(const std::vector<std::vector<double>>& X,
                       const std::vector<double>& y)
{
    // 1) Identify unique classes
    {
        std::vector<double> tmp = y;
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        classes_ = tmp;
    }
    n_classes_ = classes_.size();
    size_t n_samples = X.size();
    if (n_samples == 0) {
        std::cerr << "No training data.\n";
        return;
    }
    n_features_ = X[0].size();

    class_prior_.assign(n_classes_, 0.0);
    theta_.assign(n_classes_ * n_features_, 0.0);
    sigma_.assign(n_classes_ * n_features_, 0.0);

    // 2) Map class -> index
    std::unordered_map<double, size_t> class_index;
    for (size_t i = 0; i < n_classes_; ++i) {
        class_index[classes_[i]] = i;
    }

    // 3) Sum up for means
    std::vector<size_t> class_counts(n_classes_, 0);
    for (size_t i = 0; i < n_samples; ++i) {
        double label = y[i];
        size_t cidx = class_index[label];
        class_counts[cidx]++;
        for (size_t j = 0; j < n_features_; ++j) {
            theta_[cidx * n_features_ + j] += X[i][j];
        }
    }
    // 4) Compute mean
    for (size_t c = 0; c < n_classes_; ++c) {
        if (class_counts[c] > 0) {
            for (size_t j = 0; j < n_features_; ++j) {
                theta_[c * n_features_ + j] /= static_cast<double>(class_counts[c]);
            }
        }
    }

    // 5) Compute variance
    for (size_t i = 0; i < n_samples; ++i) {
        double label = y[i];
        size_t cidx = class_index[label];
        for (size_t j = 0; j < n_features_; ++j) {
            double diff = X[i][j] - theta_[cidx * n_features_ + j];
            sigma_[cidx * n_features_ + j] += diff * diff;
        }
    }
    for (size_t c = 0; c < n_classes_; ++c) {
        if (class_counts[c] > 1) {
            for (size_t j = 0; j < n_features_; ++j) {
                sigma_[c * n_features_ + j] /= static_cast<double>(class_counts[c]);
            }
        }
    }

    // 6) Class priors
    double total = static_cast<double>(n_samples);
    for (size_t c = 0; c < n_classes_; ++c) {
        class_prior_[c] = static_cast<double>(class_counts[c]) / total;
    }

    // 7) var_smoothing
    //    Multiply global max variance by var_smoothing_ as an epsilon
    double max_var = 0.0;
    for (auto v : sigma_) {
        if (v > max_var) max_var = v;
    }
    double eps = var_smoothing_ * max_var;
    for (auto &val : sigma_) {
        val += eps;
    }

    // 8) Copy host -> device
    copyToDevice();
}

std::vector<double> MyGaussianNB::predict(const std::vector<std::vector<double>>& X) const {
    // Check for empty input or model
    if (X.empty() || n_classes_ == 0 || n_features_ == 0) {
        return {};
    }

    const size_t n_samples = X.size();

    // Flatten X into 1D array in row-major
    std::vector<double> h_X(n_samples * n_features_);
    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < n_features_; ++j) {
            h_X[i*n_features_ + j] = X[i][j];
        }
    }

    // Allocate device memory for X
    double* d_X = nullptr;
    CUDA_CHECK(cudaMalloc(&d_X, n_samples * n_features_ * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_X, h_X.data(),
                          n_samples * n_features_ * sizeof(double),
                          cudaMemcpyHostToDevice));

    // Allocate output buffer
    double* d_out = nullptr;
    CUDA_CHECK(cudaMalloc(&d_out, n_samples * sizeof(double)));

    // Launch kernel
    int blockSize = 256;
    int gridSize  = (n_samples + blockSize - 1) / blockSize;

    predictKernel<<<gridSize, blockSize>>>(
        d_X, d_theta_, d_sigma_, d_class_prior_, d_classes_,
        d_out, n_samples, n_features_, n_classes_);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    std::vector<double> h_out(n_samples);
    CUDA_CHECK(cudaMemcpy(h_out.data(), d_out,
                          n_samples * sizeof(double),
                          cudaMemcpyDeviceToHost));

    // Clean up
    CUDA_CHECK(cudaFree(d_X));
    CUDA_CHECK(cudaFree(d_out));

    return h_out;
}

void MyGaussianNB::copyToDevice() {
    // Free old memory if re-fitting
    freeDeviceMemory();

    if (n_classes_ == 0 || n_features_ == 0) return;

    // Allocate
    size_t szCf = n_classes_ * n_features_; // for theta_, sigma_
    size_t szCl = n_classes_;              // for priors, classes

    CUDA_CHECK(cudaMalloc(&d_theta_,       szCf * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_sigma_,       szCf * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_class_prior_, szCl * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_classes_,     szCl * sizeof(double)));

    // Copy
    CUDA_CHECK(cudaMemcpy(d_theta_,       theta_.data(),       szCf * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_sigma_,       sigma_.data(),       szCf * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_class_prior_, class_prior_.data(), szCl * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_classes_,     classes_.data(),     szCl * sizeof(double), cudaMemcpyHostToDevice));
}

void MyGaussianNB::freeDeviceMemory() {
    if (d_theta_)       { CUDA_CHECK(cudaFree(d_theta_));       d_theta_ = nullptr; }
    if (d_sigma_)       { CUDA_CHECK(cudaFree(d_sigma_));       d_sigma_ = nullptr; }
    if (d_class_prior_) { CUDA_CHECK(cudaFree(d_class_prior_)); d_class_prior_ = nullptr; }
    if (d_classes_)     { CUDA_CHECK(cudaFree(d_classes_));     d_classes_ = nullptr; }
}

// ==============================================================================
// 6) MAIN: Full pipeline for S1 data
// ==============================================================================
int main() {
    // --- 1) Read CSVs for S1 (adjust folder_path as needed) ---
    std::string folder_path = "data/reference_samples/";
    auto X_train = readCSV(folder_path + "S1_extended_encoded_X_train.csv");
    auto X_test  = readCSV(folder_path + "S1_extended_encoded_X_test.csv");
    auto ytrain2d = readCSV(folder_path + "S1_extended_encoded_y_train.csv");
    auto ytest2d  = readCSV(folder_path + "S1_extended_encoded_y_test.csv");

    // Convert y from 2D -> 1D
    std::vector<double> y_train;
    y_train.reserve(ytrain2d.size());
    for (auto& row : ytrain2d) {
        if (!row.empty()) {
            y_train.push_back(row[0]);
        }
    }
    std::vector<double> y_test;
    y_test.reserve(ytest2d.size());
    for (auto& row : ytest2d) {
        if (!row.empty()) {
            y_test.push_back(row[0]);
        }
    }

    // Quick checks
    std::cout << "Train samples: " << X_train.size() << ", features: "
              << (X_train.empty() ? 0 : X_train[0].size()) << std::endl;
    std::cout << "Test  samples: " << X_test.size()  << ", features: "
              << (X_test.empty() ? 0 : X_test[0].size()) << std::endl;

    // --- 2) Standard Scale train/test ---
    standardScale(X_train);
    standardScale(X_test);

    // --- 3) Minimal SMOTE-like oversampling (replicate minority classes) ---
    oversampleRareClasses(X_train, y_train, 2);

    // --- 4) Fit model on CPU ---
    MyGaussianNB model;
    model.fit(X_train, y_train);
    std::cout << "Model fit complete." << std::endl;

    // --- 5) Predict on GPU in batches (or all at once) ---
    //     We'll just do it all at once here:
    auto y_pred = model.predict(X_test);

    // y_pred now contains the predicted classes for each row in X_test
    std::cout << "Prediction complete. Predicted " << y_pred.size() << " samples.\n";

    // Optional: Inspect first 5 predictions
    for (size_t i = 0; i < std::min<size_t>(5, y_pred.size()); ++i) {
        std::cout << "Test sample " << i << ": predicted class " << y_pred[i]
                  << ", actual class " << y_test[i] << std::endl;
    }

    // DONE
    return 0;
}
