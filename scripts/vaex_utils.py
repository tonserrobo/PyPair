"""Small helpers for Vaex/local dataframes used by the chromosome scripts."""

from __future__ import annotations

import numpy as np


def to_numpy(value):
    """Convert a Vaex DataFrame, expression, pandas object or array to numpy."""
    if hasattr(value, "to_pandas_df"):
        return value.to_pandas_df().to_numpy()
    if hasattr(value, "to_pandas_series"):
        return value.to_pandas_series().to_numpy()
    if hasattr(value, "to_numpy"):
        return value.to_numpy()
    return np.asarray(value)
