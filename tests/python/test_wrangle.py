import math

import polars as pl
import pytest

from delt_pipeline.wrangle import normalise_columns, normalise_name, wrangle


@pytest.mark.parametrize(
    ("raw", "expected"),
    [("Sample ID", "Sample_ID"), (" TP53 ", "TP53"), ("logP (calc.)", "logP_calc"), ("%%", "column")],
)
def test_normalise_name(raw: str, expected: str) -> None:
    assert normalise_name(raw) == expected


def test_normalise_columns_deduplicates() -> None:
    frame = pl.DataFrame({"a b": [1], "a-b": [2], "a_b": [3]})
    assert normalise_columns(frame).columns == ["a_b", "a_b_2", "a_b_3"]


def test_wrangle_cleans_table(raw_frame: pl.DataFrame) -> None:
    result = wrangle(raw_frame, "Sample ID", "label")
    qc = result.qc
    frame = result.table.frame
    assert frame.columns == ["Sample_ID", "smiles", "feat_alpha", "feat_beta", "gamma", "label"]
    assert qc["dropped_features"] == {"constant": "constant", "sparse": "too_many_missing"}
    assert qc["missing_values"]["feat_alpha"] == 1
    assert frame.get_column("feat_alpha").null_count() == 0
    assert frame.get_column("Sample_ID").is_sorted()
    assert (qc["rows_in"], qc["rows_out"], qc["task"]) == (24, 24, "classification")


def test_wrangle_impute_strategies(raw_frame: pl.DataFrame) -> None:
    alpha = raw_frame.get_column("feat alpha").drop_nans()
    median = wrangle(raw_frame, "Sample ID", impute="median").table.frame.get_column("feat_alpha")[3]
    mean = wrangle(raw_frame, "Sample ID", impute="mean").table.frame.get_column("feat_alpha")[3]
    dropped = wrangle(raw_frame, "Sample ID", impute="drop")
    assert math.isclose(median, alpha.median(), abs_tol=1e-6)
    assert math.isclose(mean, alpha.mean(), abs_tol=1e-6)
    assert dropped.qc["rows_dropped_missing_features"] == 1
    assert dropped.table.frame.height == 23


def test_wrangle_drops_rows_without_id_or_target() -> None:
    frame = pl.DataFrame({"id": ["a", None, "c", "d"], "x": [1.0, 2.0, 3.0, 4.0], "y": ["p", "q", None, "q"]})
    result = wrangle(frame, "id", "y")
    assert result.table.ids == ["a", "d"]


@pytest.mark.parametrize(
    ("frame", "kwargs", "message"),
    [
        (pl.DataFrame({"id": ["a"], "x": [1.0]}), {"impute": "zero"}, "unknown impute strategy"),
        (pl.DataFrame({"id": ["a"], "x": [1.0]}), {"target_column": "y"}, "column\\(s\\) not found: y"),
        (pl.DataFrame({"id": ["a", "a"], "x": [1.0, 2.0]}), {}, "duplicate identifiers"),
        (pl.DataFrame({"id": [None], "x": [1.0]}), {}, "no rows left"),
        (pl.DataFrame({"id": ["a", "b"], "x": [1.0, 1.0]}), {}, "no usable numeric feature"),
        (
            pl.DataFrame({"id": ["a", "b", "c", "d"], "x": [1.0, None, 2.0, None], "z": [None, 2.0, None, 3.0]}),
            {"impute": "drop", "max_missing": 0.9},
            "no rows left after dropping",
        ),
    ],
)
def test_wrangle_errors(frame: pl.DataFrame, kwargs: dict, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        wrangle(frame, "id", **kwargs)
