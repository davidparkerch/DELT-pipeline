import zipfile
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from sklearn.datasets import load_svmlight_file

from delt_pipeline.io import EXPORT_FORMATS, FeatureTable, export_table, parse_formats, read_table


@pytest.mark.parametrize(
    ("suffix", "writer"),
    [
        (".csv", lambda frame, path: frame.write_csv(path)),
        (".tsv", lambda frame, path: frame.write_csv(path, separator="\t")),
        (".parquet", lambda frame, path: frame.write_parquet(path)),
        (".arrow", lambda frame, path: frame.write_ipc(path)),
        (".jsonl", lambda frame, path: frame.write_ndjson(path)),
    ],
)
def test_read_table_formats(tmp_path: Path, suffix: str, writer) -> None:
    frame = pl.DataFrame({"id": ["a", "b"], "x": [1.5, 2.5]})
    path = tmp_path / f"table{suffix}"
    writer(frame, path)
    assert_frame_equal(read_table(path), frame)


def test_read_table_treats_na_as_missing(tmp_path: Path) -> None:
    path = tmp_path / "table.csv"
    path.write_text("id,x\na,NA\nb,2\n")
    assert read_table(path).get_column("x").null_count() == 1


def test_read_table_rejects_unknown_format(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unsupported table format"):
        read_table(tmp_path / "table.xlsx")


def test_feature_table_columns(table: FeatureTable) -> None:
    assert table.features == ["feat_alpha", "feat_beta", "gamma"]
    assert table.annotations == ["smiles"]
    assert table.ids[:2] == ["s00", "s01"]
    assert table.matrix().shape == (24, 3)


def test_classification_target(table: FeatureTable) -> None:
    y, classes = table.encoded_target()
    assert table.task == "classification"
    assert classes == ["active", "inactive"]
    assert set(np.unique(y)) == {0, 1}


def test_regression_target(regression: FeatureTable) -> None:
    y, classes = regression.encoded_target()
    assert regression.task == "regression"
    assert classes is None
    assert y.dtype == np.float64


def test_no_target(unlabelled: FeatureTable) -> None:
    assert unlabelled.task == "none"
    assert unlabelled.encoded_target() == (None, None)


def test_parse_formats() -> None:
    assert parse_formats(" Parquet, tsv ,parquet,,npz") == ["parquet", "tsv", "npz"]


@pytest.mark.parametrize(("text", "message"), [(" , ", "at least one"), ("parquet,xlsx", "unknown export format")])
def test_parse_formats_errors(text: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_formats(text)


def test_export_every_format(tmp_path: Path, table: FeatureTable) -> None:
    paths = export_table(table, "ds", list(EXPORT_FORMATS), tmp_path)
    assert [path.name for path in paths] == [f"ds{fmt.suffix}" for fmt in EXPORT_FORMATS.values()]
    for name in ("parquet", "arrow", "csv", "tsv", "jsonl"):
        assert_frame_equal(read_table(tmp_path / f"ds{EXPORT_FORMATS[name].suffix}"), table.frame, check_dtypes=False)

    arrays = np.load(tmp_path / "ds.npz")
    assert arrays["X"].shape == (24, 3)
    assert list(arrays["feature_names"]) == table.features
    assert list(arrays["classes"]) == ["active", "inactive"]
    assert len(arrays["y"]) == len(arrays["ids"]) == 24
    with zipfile.ZipFile(tmp_path / "ds.npz") as archive:
        assert {info.date_time for info in archive.infolist()} == {(1980, 1, 1, 0, 0, 0)}

    X, y = load_svmlight_file(str(tmp_path / "ds.svmlight"), zero_based=True)  # noqa: N806
    assert X.shape == (24, 3)
    assert set(y) == {0.0, 1.0}


def test_export_is_reproducible(tmp_path: Path, table: FeatureTable) -> None:
    first = export_table(table, "a", list(EXPORT_FORMATS), tmp_path)
    second = export_table(table, "b", list(EXPORT_FORMATS), tmp_path)
    assert [path.read_bytes() for path in first] == [path.read_bytes() for path in second]


def test_export_without_target(tmp_path: Path, unlabelled: FeatureTable) -> None:
    export_table(unlabelled, "ds", ["npz", "svmlight"], tmp_path)
    assert set(np.load(tmp_path / "ds.npz").files) == {"X", "ids", "feature_names"}
    _, y = load_svmlight_file(str(tmp_path / "ds.svmlight"), zero_based=True)
    assert not y.any()


def test_export_regression_npz(tmp_path: Path, regression: FeatureTable) -> None:
    export_table(regression, "ds", ["npz"], tmp_path)
    assert set(np.load(tmp_path / "ds.npz").files) == {"X", "ids", "feature_names", "y"}
