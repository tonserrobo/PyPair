# binary model benchmark code
import hashlib
from memory_profiler import memory_usage
import time
import numpy as np
import pandas as pd

# Benchmark Seed Lookup Aligner with Profiling
class SeedAlignerBenchmark:
    def __init__(self, seed_sequences):
        """
        Initialize the seed aligner with a list of seed sequences.
        """
        # Convert sequences into a hash table for fast lookup
        self.seed_dict = {self._hash_sequence(seq): seq for seq in seed_sequences}

    def _hash_sequence(self, sequence):
        """
        Hash a sequence for consistent lookups.
        """
        return hashlib.sha1(sequence.encode()).hexdigest()

    def align(self, query_sequence):
        """
        Attempt to align a query sequence with known seeds by binary lookup.
        Returns True if an exact match is found, else False.
        """
        query_hash = self._hash_sequence(query_sequence)
        return query_hash in self.seed_dict

    def approximate_align(self, query_sequence, threshold=0.8):
        """
        Attempt to align a query sequence approximately by calculating the percentage match.
        Returns True if match score meets or exceeds threshold, else False.
        """
        for seed in self.seed_dict.values():
            match_score = self._calculate_match_score(query_sequence, seed)
            if match_score >= threshold:
                return True  # Found a match that meets the threshold
        return False  # No match found above the threshold

    def _calculate_match_score(self, seq1, seq2):
        """
        Calculate a basic match score between two sequences as the proportion of matching characters.
        """
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / max(len(seq1), len(seq2))

    def evaluate(self, query_sequences, mode="exact", threshold=0.8, batch_size=1000):
        """
        Evaluate memory usage and execution time for a batch of query sequences.
        """
        predictions = []
        total_memory = 0
        total_time = 0

        for i in range(0, len(query_sequences), batch_size):
            batch = query_sequences[i:i + batch_size]
            # Measure memory and time for the current batch
            try:
                if mode == "exact":
                    mem_usage = memory_usage((self._batch_align, (batch,)), interval=0.1)
                elif mode == "approximate":
                    mem_usage = memory_usage((self._batch_approximate_align, (batch, threshold)), interval=0.1)
                batch_mem_usage = max(mem_usage) - min(mem_usage)
            except Exception as e:
                batch_mem_usage = np.nan

            total_memory += batch_mem_usage
            start_time = time.time()

            # Perform alignment based on the mode
            if mode == "exact":
                preds = self._batch_align(batch)
            elif mode == "approximate":
                preds = self._batch_approximate_align(batch, threshold)

            batch_time = time.time() - start_time
            total_time += batch_time

            predictions.extend(preds)

        # Metrics Capture
        avg_memory = total_memory / len(query_sequences)
        avg_time_per_query = total_time / len(query_sequences)

        metrics = {
            'Memory (MB)': total_memory,
            'Execution Time (s)': total_time,
            'Memory (MB) per query': avg_memory,
            'Execution Time (s) per query': avg_time_per_query,
        }

        return metrics, predictions

    def _batch_align(self, batch):
        """Batch processing for exact matches."""
        return [self.align(query) for query in batch]

    def _batch_approximate_align(self, batch, threshold):
        """Batch processing for approximate matches."""
        return [self.approximate_align(query, threshold) for query in batch]


test_data = pd.read_csv('research/multimodel_training/data/reference_genomes/indexed_reference.csv')
seed_sequences = test_data['suffix_array']
query_sequences = test_data['k-seed']


aligner_benchmark = SeedAlignerBenchmark(seed_sequences)

# Run exact alignment benchmark
exact_metrics, exact_predictions = aligner_benchmark.evaluate(query_sequences, mode="exact")
print("Exact Match Metrics:", exact_metrics)

# Run approximate alignment benchmark with threshold
approximate_metrics, approximate_predictions = aligner_benchmark.evaluate(query_sequences, mode="approximate", threshold=0.8)
print("Approximate Match Metrics:", approximate_metrics)

# would be cool to see if we could capture the vairance of the results etc. So maybe record results per run or something?
# what other metrics would be useful to capture?? 
