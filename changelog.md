# [0.1.0] - 13/01/2024
## Added 
- Added PCA analysis to support data structure [done]
- Added XAI (using SHAP) within multimodal training to support findings from PCA [done]
- Added binary hash based seed mapping [Benchmark][need to fully implement]
- Added native Python GNB algorithm and compared results to samples and full genome [GNB-Native][]
- Added native Python out-of-core computed GNB algorithm and compared results to samples and full genome [GNB-Native-HDF5][]
- Added hi-memory implementation to support full genome testing of scikit learn based GNB model + profiled []
- Added python 3.9 and Vaex HDF5 data loading for full genome mapping []
- Added vaex==4.17.0 and polars==1.19.0 to the requirements.txt [done]
- Added researchpy==0.3.6 to requirements.txt [done]
- Added 'old' directory to archive deprecated code base [done]

## Changed
- changed project repo structure [done]

## Fixed 

## Deprecated
- multi_model_training.ipynb deprecated and replaced with contents of dir 'multi_model_training' [done]
- predict_seeds.JIT.py deprecated added to 'old' [done]

## Removed
- Removed scikit_learn==1.2.2 from requirements.txt [done]

## Security
- Security vulnerability in scikit learn==1.2.2, removed and updated to 1.6.1 







More detail memory utilisation capture and incorporation of compressed data structures (.HDF5, as opposed to internal .csv) - for simplicity I've included the HDF5 stuff as a separate permutation of the code base - all GNB-scikit, GNB-Native, and GNB-Native-HDF5 are profiled with GNB-Native-HDF5 exceeding previous implementations. [GNB-Native, GNB-Native-HDF5, GNB-scikit comparasion] 