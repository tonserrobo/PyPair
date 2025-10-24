#!/bin/bash

###############################################################################
# run_all.sh
# Orchestrates multiple separate jobs on Slurm for data prep and model training.
###############################################################################

# Submit data prep job (multimodel_training.py)
jobid_data=$(sbatch data_prep.sh | awk '{print $4}')
if [[ -z "$jobid_data" ]]; then
    echo "Error: Data preparation job submission failed."
    exit 1
fi
echo "Submitted data prep job with ID: $jobid_data"

# Submit nearest_centroid job
jobid_nc=$(sbatch --dependency=afterok:${jobid_data} nearest_centroid.sh | awk '{print $4}')
if [[ -z "$jobid_nc" ]]; then
    echo "Error: nearest_centroid job submission failed."
    exit 1
fi
echo "Submitted nearest_centroid job with ID: $jobid_nc"

# Submit kneighbours_classifier job
jobid_kn=$(sbatch --dependency=afterok:${jobid_data} kNeighbours_classifier.sh | awk '{print $4}')
if [[ -z "$jobid_kn" ]]; then
    echo "Error: kNeighbours_classifier job submission failed."
    exit 1
fi
echo "Submitted kNeighbours_classifier job with ID: $jobid_kn"

# Submit random_forest job
jobid_rf=$(sbatch --dependency=afterok:${jobid_data} random_forest.sh | awk '{print $4}')
if [[ -z "$jobid_rf" ]]; then
    echo "Error: random_forest job submission failed."
    exit 1
fi
echo "Submitted random_forest job with ID: $jobid_rf"

# Submit extratrees_classifier job
jobid_et=$(sbatch --dependency=afterok:${jobid_data} extratrees_classifier.sh | awk '{print $4}')
if [[ -z "$jobid_et" ]]; then
    echo "Error: extratrees_classifier job submission failed."
    exit 1
fi
echo "Submitted extratrees_classifier job with ID: $jobid_et"

# Submit DecisionTrees job
jobid_dt=$(sbatch --dependency=afterok:${jobid_data} decisiontrees_classifier.sh | awk '{print $4}')
if [[ -z "$jobid_dt" ]]; then
    echo "Error: decisiontrees_classifier job submission failed."
    exit 1
fi
echo "Submitted decisiontrees_classifier job with ID: $jobid_dt"

# Submit GNaiveBayes job
jobid_nb=$(sbatch --dependency=afterok:${jobid_data} gaussianNaiveBayes_classifier.sh | awk '{print $4}')
if [[ -z "$jobid_nb" ]]; then
    echo "Error: GNaiveBayes_classifier job submission failed."
    exit 1
fi
echo "Submitted GNaiveBayes_classifier job with ID: $jobid_nb"

# Submit LGBM Classifier job
jobid_lg=$(sbatch --dependency=afterok:${jobid_data} LGBM_classifier.sh | awk '{print $4}')
if [[ -z "$jobid_lg" ]]; then
    echo "Error: LGBM Classifier job submission failed."
    exit 1
fi
echo "Submitted LGBM Classifier job with ID: $jobid_lg"

# Submit XGB Classifier job
jobid_xg=$(sbatch --dependency=afterok:${jobid_data} XGB_classifier.sh | awk '{print $4}')
if [[ -z "$jobid_xg" ]]; then
    echo "Error: XGB classifier job submission failed."
    exit 1
fi
echo "Submitted XGB classifier job with ID: $jobid_xg"

echo "All jobs submitted successfully."

# Cleanup Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda deactivate
echo "Conda environment deactivated."
