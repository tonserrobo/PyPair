import math

from scripts.noninferiority_testing import noninferiority_ttest_from_summary


def test_non_inferior_when_ci_below_threshold():
    result = noninferiority_ttest_from_summary(
        mean=1.0, std=0.1, n=5, bwa_time=1.0, margin=0.3
    )

    assert result["threshold"] == 1.3
    assert result["p_value"] < 0.01
    assert result["ci_upper_bound"] < 1.3
    assert result["non_inferior"] is True


def test_not_non_inferior_when_mean_is_above_threshold():
    result = noninferiority_ttest_from_summary(
        mean=1.4, std=0.1, n=5, bwa_time=1.0, margin=0.3
    )

    assert result["p_value"] > 0.05
    assert result["ci_upper_bound"] > 1.3
    assert result["non_inferior"] is False


def test_superior_when_ci_is_below_bwa():
    result = noninferiority_ttest_from_summary(
        mean=0.5, std=0.05, n=5, bwa_time=1.0, margin=0.3
    )

    assert result["ci_upper_bound"] < 1.0
    assert result["superior"] is True
