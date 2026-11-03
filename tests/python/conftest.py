"""Shared fixtures: small, deterministic feature tables."""

import numpy as np
import polars as pl
import pytest

from delt_pipeline.io import FeatureTable
from delt_pipeline.wrangle import wrangle

N_ROWS = 24


@pytest.fixture
def raw_frame() -> pl.DataFrame:
    """A messy table: odd column names, a gap, a constant and a mostly-missing column."""
    rng = np.random.default_rng(7)
    label = np.array(["active", "inactive"] * (N_ROWS // 2))
    signal = np.where(label == "active", 2.0, -2.0)
    alpha = signal + rng.normal(0, 0.5, N_ROWS)
    alpha[3] = np.nan
    return pl.DataFrame(
        {
            "Sample ID": [f"s{index:02d}" for index in range(N_ROWS)],
            "smiles": ["CCO"] * N_ROWS,
            "feat alpha": alpha,
            "feat-beta": rng.normal(10, 2, N_ROWS),
            "gamma": -signal + rng.normal(0, 1, N_ROWS),
            "constant": [1.0] * N_ROWS,
            "sparse": [1.0, 2.0] + [None] * (N_ROWS - 2),
            "label": label,
        }
    )


@pytest.fixture
def table(raw_frame: pl.DataFrame) -> FeatureTable:
    """Cleaned binary-classification table."""
    return wrangle(raw_frame, "Sample ID", "label").table


@pytest.fixture
def unlabelled(raw_frame: pl.DataFrame) -> FeatureTable:
    """Cleaned table without a target."""
    return wrangle(raw_frame.drop("label"), "Sample ID").table


@pytest.fixture
def regression() -> FeatureTable:
    """Table with a continuous target."""
    rng = np.random.default_rng(11)
    x = rng.normal(0, 1, (N_ROWS, 3))
    frame = pl.DataFrame(
        {
            "id": [f"r{index:02d}" for index in range(N_ROWS)],
            "x1": x[:, 0],
            "x2": x[:, 1],
            "x3": x[:, 2],
            "y": 3 * x[:, 0] + rng.normal(0, 0.1, N_ROWS),
        }
    )
    return FeatureTable(frame, "id", "y")
