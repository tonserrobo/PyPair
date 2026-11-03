import numpy as np
import pandas as pd

from scripts.vaex_utils import to_numpy


def test_to_numpy_from_dataframe_like():
    class FrameLike:
        def to_pandas_df(self):
            return pd.DataFrame({"a": [1, 2], "b": [3, 4]})

    result = to_numpy(FrameLike())

    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 2)


def test_to_numpy_from_series_like():
    class SeriesLike:
        def to_pandas_series(self):
            return pd.Series([1, 2, 3])

    result = to_numpy(SeriesLike())

    assert isinstance(result, np.ndarray)
    assert result.tolist() == [1, 2, 3]


def test_to_numpy_from_plain_list():
    result = to_numpy([1, 2, 3])

    assert isinstance(result, np.ndarray)
    assert result.tolist() == [1, 2, 3]
