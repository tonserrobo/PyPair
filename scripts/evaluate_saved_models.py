#!/usr/bin/env python3
"""Evaluate the saved per-sample models on the cleaned test sets.

This exists because the older `KNeighborsClassifier_results.csv` was a stale
copy of the GNB metrics. The saved `*_s1..s4.joblib` models are genuine and
were trained by `scripts/comprehensive_model_analysis.py` on the cleaned data.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import joblib
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score

DATASETS = ["S1", "S2", "S3", "S4"]
MODELS = ["knn", "gaussian_nb"]


def evaluate(dataset: str, model_name: str, data_dir: Path, weights_dir: Path) -> dict | None:
    model_path = weights_dir / f"{model_name}_{dataset.lower()}.joblib"
    if not model_path.exists():
        print(f"  skip {model_name} {dataset}: {model_path} not found")
        return None

    X_test = pd.read_csv(data_dir / f"{dataset}_extended_encoded_X_test.csv")
    y_test = pd.read_csv(data_dir / f"{dataset}_extended_encoded_y_test.csv")["org_location"]

    model = joblib.load(model_path)
    start = time.perf_counter()
    pred = model.predict(X_test.values)
    elapsed = time.perf_counter() - start

    return {
        "model": model_name,
        "dataset": dataset,
        "n_test": len(X_test),
        "n_features": X_test.shape[1],
        "test_f1": f1_score(y_test, pred, average="weighted"),
        "test_precision": precision_score(y_test, pred, average="weighted"),
        "test_recall": recall_score(y_test, pred, average="weighted"),
        "test_accuracy": accuracy_score(y_test, pred),
        "predict_time_ms_per_sample": elapsed / len(X_test) * 1000,
        "model_size_mb": model_path.stat().st_size / (1024 ** 2),
        "protocol": "cleaned_test_set_saved_model",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate saved KNN/GNB models on cleaned test sets")
    parser.add_argument("--data-dir", type=Path, default=Path("data/reference_samples/cleaned"))
    parser.add_argument("--weights-dir", type=Path, default=Path("models/model_weights"))
    parser.add_argument("--out", type=Path, default=Path("results/cleaned_model_evaluation.csv"))
    args = parser.parse_args(argv)

    rows = []
    for dataset in DATASETS:
        for model_name in MODELS:
            print(f"evaluating {model_name} {dataset}")
            row = evaluate(dataset, model_name, args.data_dir, args.weights_dir)
            if row:
                rows.append(row)

    if not rows:
        raise SystemExit("no models evaluated")

    frame = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)
    print(frame[["model", "dataset", "test_f1", "predict_time_ms_per_sample"]].to_string(index=False))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
