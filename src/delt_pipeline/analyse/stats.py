"""Descriptive and inferential statistics with scipy."""

from typing import Any

import numpy as np
from scipy import stats

from .._numbers import significant
from ..io import FeatureTable

SHAPIRO_MIN, SHAPIRO_MAX = 3, 5000


def describe_features(table: FeatureTable) -> list[dict[str, Any]]:
    """Summarise every feature: location, spread, shape and a Shapiro-Wilk normality p-value."""
    matrix = table.matrix()
    summaries = []
    for index, name in enumerate(table.features):
        values = matrix[:, index]
        varies = bool(np.ptp(values) > 0)
        normal_p = None
        if varies and SHAPIRO_MIN <= values.size <= SHAPIRO_MAX:
            normal_p = significant(stats.shapiro(values).pvalue)
        summaries.append(
            {
                "feature": name,
                "mean": significant(values.mean()),
                "std": significant(values.std(ddof=1)) if values.size > 1 else None,
                "min": significant(values.min()),
                "median": significant(np.median(values)),
                "max": significant(values.max()),
                "skewness": significant(stats.skew(values)) if varies else None,
                "shapiro_p": normal_p,
            }
        )
    return summaries


def correlation(table: FeatureTable) -> dict[str, Any]:
    """Spearman rank correlation between all features."""
    matrix = table.matrix()
    if matrix.shape[1] == 1:
        rho = np.ones((1, 1))
    else:
        rho = np.atleast_2d(stats.spearmanr(matrix).statistic)
        if rho.shape == (1, 1):  # scipy returns a scalar for exactly two variables
            rho = np.array([[1.0, rho[0, 0]], [rho[0, 0], 1.0]])
    return {
        "method": "spearman",
        "features": table.features,
        "matrix": [[significant(value) for value in row] for row in rho],
    }


def association(table: FeatureTable) -> list[dict[str, Any]]:
    """Test each feature against the target.

    Classification targets use a Kruskal-Wallis H-test across classes; regression
    targets use Spearman correlation. Returns an empty list without a target.
    """
    y, _ = table.encoded_target()
    if y is None:
        return []
    matrix = table.matrix()
    results = []
    for index, name in enumerate(table.features):
        values = matrix[:, index]
        if table.task == "classification":
            groups = [values[y == code] for code in np.unique(y)]
            test = "kruskal"
            statistic, p_value = float("nan"), float("nan")
            if len(groups) > 1 and np.ptp(values) > 0:  # the H-test is undefined otherwise
                statistic, p_value = stats.kruskal(*groups)
        else:
            test = "spearman"
            result = stats.spearmanr(values, y)
            statistic, p_value = result.statistic, result.pvalue
        results.append(
            {"feature": name, "test": test, "statistic": significant(statistic), "p_value": significant(p_value)}
        )
    return results
