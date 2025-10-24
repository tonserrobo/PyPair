#!/usr/bin/env python3
"""
Clean S1 dataset by converting string-encoded features to numeric types.

This script fixes data type issues in the S1 extended encoded reference data:
- Converts binary string features to integers
- Hashes overly long string features
- Converts numeric strings to proper int/float types
- Ensures all features are ML-ready

Author: Tony Robinson (t.robinson@ulster.ac.uk)
Date: 2025-01-24
"""

import pandas as pd
import numpy as np
from pathlib import Path
import hashlib


def binary_string_to_int(binary_str):
    """Convert binary string like '110110...' to integer."""
    if pd.isna(binary_str) or binary_str == '':
        return 0
    try:
        # Limit to first 63 bits to avoid overflow (int64 max)
        binary_str = str(binary_str)[:63]
        return int(binary_str, 2)
    except (ValueError, TypeError):
        return 0


def hash_string_to_int(string_val):
    """Hash a string to a 64-bit integer."""
    if pd.isna(string_val) or string_val == '':
        return 0
    try:
        # Use first 16 hex digits of MD5 hash (64 bits)
        hash_obj = hashlib.md5(str(string_val).encode())
        hash_hex = hash_obj.hexdigest()[:16]
        return int(hash_hex, 16) % (2**63)  # Keep within int64 range
    except (ValueError, TypeError):
        return 0


def numeric_string_to_int(num_str):
    """Convert numeric string to integer, handling overflow."""
    if pd.isna(num_str) or num_str == '':
        return 0
    try:
        # Truncate very long numbers to prevent overflow
        num_str = str(num_str)[:18]  # Max ~18 digits for int64
        return int(num_str)
    except (ValueError, TypeError, OverflowError):
        return 0


def clean_dataframe(df, dataset_name="dataset"):
    """Clean a single DataFrame by converting string features to numeric."""
    print(f"\nCleaning {dataset_name}...")
    print(f"  Original shape: {df.shape}")
    print(f"  Original dtypes:")
    for col in df.columns:
        print(f"    {col}: {df[col].dtype}")

    df_clean = df.copy()

    # Convert binary string features to integers
    if 'rc_seeds' in df_clean.columns:
        print("  Converting 'rc_seeds' (binary string -> int64)...")
        df_clean['rc_seeds'] = df_clean['rc_seeds'].apply(binary_string_to_int)

    if 'runlength_rc_character' in df_clean.columns and df_clean['runlength_rc_character'].dtype == 'object':
        print("  Converting 'runlength_rc_character' (binary string -> int64)...")
        df_clean['runlength_rc_character'] = df_clean['runlength_rc_character'].apply(binary_string_to_int)

    # Hash long string features that would overflow int64
    if 'runlength_output' in df_clean.columns and df_clean['runlength_output'].dtype == 'object':
        print("  Converting 'runlength_output' (hash -> int64)...")
        df_clean['runlength_output'] = df_clean['runlength_output'].apply(hash_string_to_int)

    if 'runlength_rc_output' in df_clean.columns and df_clean['runlength_rc_output'].dtype == 'object':
        print("  Converting 'runlength_rc_output' (hash -> int64)...")
        df_clean['runlength_rc_output'] = df_clean['runlength_rc_output'].apply(hash_string_to_int)

    # Convert numeric string features to integers
    if 'runlength_value' in df_clean.columns and df_clean['runlength_value'].dtype == 'object':
        print("  Converting 'runlength_value' (numeric string -> int64)...")
        df_clean['runlength_value'] = df_clean['runlength_value'].apply(numeric_string_to_int)

    if 'runlength_rc_value' in df_clean.columns and df_clean['runlength_rc_value'].dtype == 'object':
        print("  Converting 'runlength_rc_value' (numeric string -> int64)...")
        df_clean['runlength_rc_value'] = df_clean['runlength_rc_value'].apply(numeric_string_to_int)

    # Ensure numeric columns are proper types
    if 'first' in df_clean.columns:
        df_clean['first'] = df_clean['first'].astype('int32')
    if 'A' in df_clean.columns:
        df_clean['A'] = df_clean['A'].astype('int32')
    if 'C' in df_clean.columns:
        df_clean['C'] = df_clean['C'].astype('int32')
    if 'G' in df_clean.columns:
        df_clean['G'] = df_clean['G'].astype('int32')
    if 'T' in df_clean.columns:
        df_clean['T'] = df_clean['T'].astype('int32')

    # Handle target variable if present
    if 'org_location' in df_clean.columns:
        print("  Converting 'org_location' (float -> int32)...")
        df_clean['org_location'] = df_clean['org_location'].astype('int32')

    print(f"  Cleaned shape: {df_clean.shape}")
    print(f"  Cleaned dtypes:")
    for col in df_clean.columns:
        print(f"    {col}: {df_clean[col].dtype}")

    # Check for any remaining non-numeric types
    object_cols = df_clean.select_dtypes(include=['object']).columns
    if len(object_cols) > 0:
        print(f"  WARNING: {len(object_cols)} columns still have object dtype: {list(object_cols)}")

    # Check for inf/nan values
    inf_counts = np.isinf(df_clean.select_dtypes(include=[np.number])).sum()
    if inf_counts.any():
        print(f"  WARNING: Found infinity values:")
        for col, count in inf_counts[inf_counts > 0].items():
            print(f"    {col}: {count} inf values")

    return df_clean


def main():
    print("="*60)
    print("S1-S4 Dataset Cleaning")
    print("="*60)

    data_path = Path("data/reference_samples")
    output_path = Path("data/reference_samples/cleaned")
    output_path.mkdir(exist_ok=True)

    # Clean all datasets for S1-S4
    for sample_num in ['S1', 'S2', 'S3', 'S4']:
        print(f"\n{'#'*60}")
        print(f"# Cleaning {sample_num} Dataset")
        print(f"{'#'*60}")

        datasets = {
            'x_train': f'{sample_num}_extended_encoded_x_train.csv',
            'X_test': f'{sample_num}_extended_encoded_X_test.csv',
            'y_train': f'{sample_num}_extended_encoded_y_train.csv',
            'y_test': f'{sample_num}_extended_encoded_y_test.csv',
            'X_val': f'{sample_num}_extended_encoded_X_val.csv',
            'y_val': f'{sample_num}_extended_encoded_y_val.csv'
        }

        for dataset_name, filename in datasets.items():
            input_file = data_path / filename
            output_file = output_path / filename

            if not input_file.exists():
                print(f"\nSkipping {filename} (not found)")
                continue

            print(f"\n{'='*60}")
            print(f"Processing: {sample_num} {dataset_name}")
            print(f"{'='*60}")

            # Load
            print(f"Loading {input_file}...")
            df = pd.read_csv(input_file)

            # Clean
            df_clean = clean_dataframe(df, f"{sample_num}_{dataset_name}")

            # Save
            print(f"Saving to {output_file}...")
            df_clean.to_csv(output_file, index=False)
            print(f"[OK] Saved cleaned {sample_num} {dataset_name}")

    print(f"\n{'='*60}")
    print("CLEANING COMPLETE!")
    print(f"{'='*60}")
    print(f"Cleaned datasets saved to: {output_path}")
    print(f"\nNext step: Run comprehensive training on all cleaned data:")
    print(f"  python scripts/comprehensive_model_analysis.py")


if __name__ == "__main__":
    main()
