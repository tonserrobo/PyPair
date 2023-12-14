# PyPair: Python-Based Short Read Alignment Software

## Introduction
PyPair is a Python-based software tool for aligning short DNA sequences (reads) to a reference genome. It employs a learned index for exact-match seed generation, enhancing the efficiency and accuracy of the alignment process. The software integrates several bioinformatics tools and techniques, such as Biopython for parsing sequence data and PySAM for alignment handling, combined with machine learning models for improved seed generation.

## Getting Started
### Installation
To install PyPair, clone the repository and install the required dependencies:

```
git clone https://github.com/your-repository/PyPair.git` 
cd PyPair
pip install -r requirements.txt
```
### Usage
PyPair operates through a series of steps to align reads to a reference genome:

1. *FASTA Reference and FASTQ Reads Parsing:* Utilizes Biopython to read FASTA formatted reference genomes and FASTQ formatted read files.
2. *Seed Generation and Mapping:* Generates mapping locations from FASTQ seeds and builds them into an internal dictionary.
3. *Sequence Alignment:* Employs the Smith-Waterman algorithm via PySAM for precise sequence alignment.
   
To start the seed generation and mapping process, run the following:

```
python predict_seeds.py <path-to-FASTQ-file>
```
### Script Overview
The main script, predict_seeds.py, orchestrates the alignment process. It includes:

1. Parsing of FASTQ files for read sequences.
2. A Smith-Waterman algorithm implementation for sequence alignment.
3. CIGAR string computation from sequence alignments.
4. Model loading and sequence processing utilities.
5. Generation of BAM files from sequence alignments.
   
### Dependencies
PyPair requires Python 3.x and several libraries, including Biopython, PySAM, Pandas, and Pickle. Ensure these are installed as per 1requirements.txt`.

### Project Structure
The PyPair project is structured as follows:
```
PyPair/
├── model_weights/         # Trained model weights
├── reference_samples/     # Reference genome samples
├── Validation/            # Validation datasets
├── wandb/                 # Weights & Biases logging (optional)
├── encode_decode.py       # Encoding/decoding utilities
├── gen_index.py           # Index generation script
├── gen_training_data.py   # Training data generation script
├── io_processing.py       # I/O processing utilities
├── multi_model_training.ipynb  # Jupyter notebook for model training
├── predict_seeds.py       # Main prediction script
├── readme                 # Project README
└── requirements.txt       # Python package requirements
```
### Support
For issues, questions, or contributions, please open an issue or pull request in the GitHub repository.
