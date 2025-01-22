import hashlib
from memory_profiler import memory_usage
import time
import numpy as np
import pandas as pd

class SeedAlignerBenchmark:
    def __init__(self, seed_sequences):
        """
        Initialize the seed aligner with a list of seed sequences.
        """
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
        Returns True if match score meets or exceeds a threshold, else False.
        """
        for seed in self.seed_dict.values():
            match_score = self._calculate_match_score(query_sequence, seed)
            if match_score >= threshold:
                return True
        return False

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
        if not query_sequences:
            raise ValueError("Query sequences list is empty.")

        batch_memory_usage = []
        batch_time_usage = []
        predictions = []

        for i in range(0, len(query_sequences), batch_size):
            batch = query_sequences[i:i + batch_size]
            try:
                if mode == "exact":
                    mem_usage = memory_usage((self._batch_align, (batch,)), interval=0.1)
                elif mode == "approximate":
                    mem_usage = memory_usage((self._batch_approximate_align, (batch, threshold)), interval=0.1)
                batch_mem_usage = max(mem_usage) if mem_usage else 0
            except Exception as e:
                batch_mem_usage = np.nan

            batch_memory_usage.append(batch_mem_usage)
            start_time = time.time()

            if mode == "exact":
                preds = self._batch_align(batch)
            elif mode == "approximate":
                preds = self._batch_approximate_align(batch, threshold)

            batch_time = time.time() - start_time
            batch_time_usage.append(batch_time)
            predictions.extend(preds)

        total_memory = sum(batch_memory_usage)
        total_time = sum(batch_time_usage)
        total_time = total_time if total_time > 0 else 1e-6  # Prevent division by zero

        metrics = {
            'Memory (MB)': total_memory,
            'Execution Time (s)': total_time,
            'Memory (MB) per query': total_memory / len(query_sequences),
            'Execution Time (s) per query': total_time / len(query_sequences),
            'Memory Variance': np.var(batch_memory_usage),
            'Execution Time Variance': np.var(batch_time_usage),
            'Memory Std Dev': np.std(batch_memory_usage),
            'Execution Time Std Dev': np.std(batch_time_usage),
            'Throughput (queries/sec)': len(query_sequences) / total_time
        }
        return pd.DataFrame([metrics]), predictions

    def _batch_align(self, batch):
        return [self.align(query) for query in batch]

    def _batch_approximate_align(self, batch, threshold):
        return [self.approximate_align(query, threshold) for query in batch]

base_path = "C:/Users/Tony/OneDrive - Ulster University (1)/Phd - FPGA acceleration in genomics/Publications + Writing/PhD_reports/PhD Thesis/cp3. Seed generation via ML/research"
test_data = pd.read_csv(base_path + '/multimodel_training/data/reference_genomes/indexed_reference.csv')
seed_sequences = test_data['suffix_array'].astype(str).tolist()
query_sequences = test_data['k-seed'].astype(str).tolist()

aligner_benchmark = SeedAlignerBenchmark(seed_sequences)

exact_metrics_df, exact_predictions = aligner_benchmark.evaluate(query_sequences, mode="exact")
print("Exact Match Metrics:")
print(exact_metrics_df)
exact_metrics_df.to_csv("exact_match_hash_benchmark.csv", index=False)

approximate_metrics_df, approximate_predictions = aligner_benchmark.evaluate(query_sequences, mode="approximate", threshold=0.8)
print("Approximate Match Metrics:")
print(approximate_metrics_df)
approximate_metrics_df.to_csv("approximate_match_hash_benchmark.csv", index=False)
