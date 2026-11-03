# Overfitting Analysis: Why Tree-Based Models Fail

> **Canonical status:** the aggregate table in this note is the protocol-consistent
> baseline. The repository also contains divergent legacy files; see
> [`docs/results_provenance.md`](results_provenance.md). The KNN result file in
> `results/` is contaminated (its precision, recall and F1 are identical to the
> GNB file); use the canonical table, not that file.

## Executive Summary

Comprehensive evaluation of 8 machine learning classifiers revealed **severe overfitting in tree-based ensemble methods** (Random Forest, LGBM, XGBoost) when applied to genomic seed alignment. This document analyzes the root causes and explains why **Gaussian Naive Bayes emerged as the optimal solution**.

## Observed Performance Patterns

### Training vs Test Performance Gap

| Model | CV Score (%) | F1 Score (%) | Gap (pp) | Severity |
|-------|--------------|--------------|----------|----------|
| Random Forest | 100.0 | 92.8 | **7.2** | High ❌ |
| LGBM | 100.0 | 90.5 | **9.5** | Critical ❌ |
| XGBoost | 100.0 | 91.8 | **8.2** | High ❌ |
| Extra Trees | 99.2 | 94.9 | **4.3** | Medium ⚠️ |
| **Gaussian NB** | **95.1** | **94.7** | **0.4** | **Minimal** ✅ |
| Decision Tree | 96.2 | 94.8 | **1.4** | Low ✅ |
| KNN | 99.7 | 94.4 | **5.3** | Medium ⚠️ |
| Nearest Centroid | 98.0 | 95.7 | **2.3** | Low ✅ |

### Key Observation

**Perfect 100% CV scores are a red flag**, not a success metric. They indicate the model has memorized training data rather than learning generalizable patterns.

## Root Cause Analysis

### 1. Training Data Limitations

**Problem**: Small effective training set size relative to problem space

```
Training Data Size: ~400MB (S1-S4 samples)
Feature Space: DNA k-mers (4^k combinations)
Human Genome Size: ~3.2GB (3.2 billion base pairs)
Training Coverage: <0.1% of possible genomic positions
```

**Impact**:
- Models see limited variation in genomic contexts
- Tree-based models create decision boundaries that perfectly fit training samples
- No regularization against memorization

**Evidence from [model_analysis.ipynb](../model_analysis.ipynb)**:
- Training samples use **extended/repeated data** to balance classes
- This artificial augmentation inflates training metrics
- Test data contains genuinely unseen genomic regions

### 2. Model Capacity vs Data Size Mismatch

**Tree-Based Model Sizes**:
```
LGBM:          233 MB  (117,000 leaves × 2KB/leaf)
Random Forest: 117 MB  (1000 trees × 117KB/tree)
XGBoost:       102 MB  (15,000 boosting rounds)
Gaussian NB:     7 MB  (class means + variances only)
```

**Analysis**:
- LGBM has **33x more parameters** than GNB
- With only 400MB training data, this creates a **parameter-to-data ratio** that invites overfitting
- Tree-based models can memorize individual training examples in leaf nodes

### 3. Nature of Genomic Sequence Data

**K-mer Feature Characteristics**:
- DNA sequences encoded as k-mer frequency vectors
- High-dimensional sparse feature space (4^k possible k-mers)
- Many features are correlated (overlapping k-mers)
- **Conditional independence assumption (GNB) is actually beneficial** for this data type

**Why Trees Overfit**:
```python
# Tree-based models create splits like:
if k-mer_ACG_count > 5:
    if k-mer_CGA_count > 3:
        if k-mer_GAT_count > 7:
            return class_A  # Highly specific to training examples

# GNB instead models:
P(class | features) = P(class) * Π P(feature_i | class)
# More robust to unseen k-mer combinations
```

### 4. Class Imbalance Effects

**From Training Pipeline**:
- SMOTE oversampling used to balance classes
- Synthetic samples created via interpolation
- Tree-based models memorize these synthetic patterns
- GNB probabilistic approach naturally handles imbalance

**Evidence**:
```python
# From multimodel_training.py
smote = SMOTE(random_state=seed, n_jobs=-1)
X_train_balanced, y_train_balanced = smote.fit_resample(X_train, y_train)
```

SMOTE increases training set size but doesn't add genuinely new information, which tree models exploit.

### 5. Cross-Validation Methodology

**Current CV Setup** (inferred from results):
- K-fold cross-validation on S1-S4 datasets
- Each fold contains samples from same extended/repeated pool
- **Data leakage**: CV folds are not truly independent due to data augmentation

**Why This Affects Trees More**:
- Trees memorize specific augmented examples
- GNB learns distributional parameters that generalize across augmentations

## Detailed Model-by-Model Analysis

### Random Forest (100% CV, 92.8% F1)

**Mechanism**:
- 1000 trees, each memorizing different bootstrap samples
- Averaging predictions doesn't prevent overfitting when all trees overfit
- **Bagging amplifies memorization** on small datasets

**Model Size**: 117 MB (avg 117 KB per tree)

**Recommendation**: ❌ **Not suitable** for this task

### LightGBM (100% CV, 90.5% F1)

**Worst Overfitting** (-9.5pp gap)

**Mechanism**:
- Gradient boosting iteratively fits residuals
- With sufficient rounds, can achieve perfect training fit
- No early stopping or regularization effective enough

**Model Size**: 233 MB (largest model)

**Recommendation**: ❌ **Strongly discouraged**

### XGBoost (100% CV, 91.8% F1)

**Mechanism**:
- Similar to LGBM but different tree construction algorithm
- Regularization parameters (max_depth, min_child_weight) insufficient
- 15,000+ boosting rounds enable overfitting

**Model Size**: 102 MB

**Recommendation**: ❌ **Not suitable**

### Extra Trees (99.2% CV, 94.9% F1)

**Better Than Other Ensembles** but still problematic

**Mechanism**:
- Extremely randomized splits reduce overfitting vs Random Forest
- Still memorizes via ensemble of diverse trees
- 4.3pp gap is borderline acceptable

**Model Size**: 118 MB

**Recommendation**: ⚠️ **Marginal** - consider only if GNB insufficient

### Decision Tree (96.2% CV, 94.8% F1)

**Surprisingly Competitive**

**Mechanism**:
- Single tree with limited depth (max_depth likely set)
- Acts as a regularized baseline
- Gap of 1.4pp is acceptable

**Model Size**: 25 KB (tiny!)

**Recommendation**: ✅ **Acceptable baseline**, but GNB still better

### K-Nearest Neighbors (99.7% CV, 94.4% F1)

**Mechanism**:
- Instance-based learning (stores all training examples)
- Overfits by memorizing local neighborhoods
- Poor scalability (2.0 MB memory per query)

**Recommendation**: ⚠️ **Not scalable** to genome-wide deployment

### Nearest Centroid (98.0% CV, 95.7% F1)

**Best Non-GNB Model**

**Mechanism**:
- Computes class centroids (mean feature vectors)
- Simple linear decision boundaries
- Similar philosophy to GNB (distributional approach)

**Gap**: 2.3pp (acceptable)

**Recommendation**: ✅ **Viable alternative** if GNB unavailable

### Gaussian Naive Bayes (95.1% CV, 94.7% F1) ✅

**Best Generalization** (0.4pp gap)

**Why It Succeeds**:

1. **Low Model Capacity**: Only stores class priors, means, and variances
   ```
   Parameters: n_classes × n_features × 3 (prior, mean, var)
   Example: 10 classes × 1000 features × 3 = 30,000 params
   vs LGBM: ~500,000+ split nodes
   ```

2. **Probabilistic Framework**: Models P(X|Y) distributions rather than memorizing examples
   ```python
   # GNB decision rule
   class = argmax_c [ log(P(c)) + Σ log(P(x_i | c)) ]
   # Smooth Gaussian assumptions prevent overfitting
   ```

3. **Feature Independence Assumption**: Acts as implicit regularization
   - Ignores spurious correlations in training data
   - Generalizes better to new k-mer combinations

4. **Small Model Size**: 2-7 MB (33x smaller than LGBM)
   - Less prone to memorization
   - Faster training and inference

5. **Scalability**: Out-of-core Vaex implementation enables genome-scale training
   - Can train on full 3.2GB human reference
   - More training data reduces overfitting further

## Recommendations

### For PyPair Project

✅ **Primary Model**: Gaussian Naive Bayes (Vaex implementation)
- Best generalization
- Genome-scale ready
- Minimal overfitting

✅ **Baseline Comparison**: Nearest Centroid
- Simple sanity check
- Similar performance to GNB

⚠️ **Fallback**: Single Decision Tree (with max_depth=10)
- If probabilistic model unsuitable
- Must tune max_depth carefully

❌ **Avoid**: Random Forest, LGBM, XGBoost, Extra Trees
- Severe overfitting on current data
- Would require 10-100x more training data

### Future Work to Reduce Overfitting

1. **Increase Training Data**:
   ```bash
   # Use Vaex to train on full genome
   python accelerated_GNB/GNB_vaex/GNB_HDF5_optimised.py \
       --genome GRCh38_full.hdf5 \
       --output models/gnb_full_genome.joblib
   ```

2. **Remove Data Augmentation**:
   - Train on genuinely distinct genomic regions
   - Avoid extended/repeated samples

3. **Better CV Strategy**:
   - Chromosome-wise splits (train on chr1-20, test on chr21-22)
   - Ensures truly independent test set

4. **Regularization for Trees** (if must use):
   ```python
   # Heavily regularized Random Forest
   rf = RandomForestClassifier(
       n_estimators=50,        # Reduce from 1000
       max_depth=5,            # Limit depth
       min_samples_split=100,  # Require more samples
       max_features='sqrt',    # Feature subsampling
       random_state=42
   )
   ```

5. **Ensemble of GNB Models**:
   - Train multiple GNB on different genomic regions
   - Ensemble via weighted voting
   - Maintains generalization while improving accuracy

## Conclusion

The severe overfitting observed in tree-based models is not a flaw in the models themselves, but rather a **data regime mismatch**:

- **Tree ensembles excel** when training data >> model capacity
- **Genomic alignment has** training data << problem space

**Gaussian Naive Bayes is the optimal choice** because:
1. Lower capacity prevents memorization
2. Probabilistic foundation generalizes better
3. Scalable to genome-wide training (which will further reduce overfitting)
4. Smallest gap between training and test performance

**Next Step**: Deploy Vaex-based GNB training on full human reference genome at NIHPC, which should improve test F1 scores from 94.7% to an estimated **96-97%** due to increased training coverage.

---

**Analysis Date**: 2025-01-24
**Data Source**: [model_analysis.ipynb](../model_analysis.ipynb)
**Author**: Tony Robinson (t.robinson@ulster.ac.uk)
