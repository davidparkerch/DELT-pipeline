"""Feature tables: reading, typing and exporting to downstream-ready formats."""

import io
import zipfile
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.datasets import dump_svmlight_file

NULL_VALUES = ["", "NA", "N/A", "NaN", "nan", "null", "NULL"]
MAX_CLASSES = 10
"""Numeric targets with at most this many distinct values are treated as classes."""

# Zip entries get a fixed timestamp so exports are byte-for-byte reproducible.
_ZIP_EPOCH = (1980, 1, 1, 0, 0, 0)


def read_table(path: Path) -> pl.DataFrame:
    """Read a CSV, TSV, Parquet, Arrow IPC or JSON Lines table.

    Args:
        path: Table file; the format is chosen from the file extension.

    Returns:
        The table as a polars DataFrame.

    Raises:
        ValueError: If the extension is not a supported table format.
    """
    suffix = path.suffix.lower()
    match suffix:
        case ".csv":
            return pl.read_csv(path, null_values=NULL_VALUES, infer_schema_length=None)
        case ".tsv" | ".tab" | ".txt":
            return pl.read_csv(path, separator="\t", null_values=NULL_VALUES, infer_schema_length=None)
        case ".parquet":
            return pl.read_parquet(path)
        case ".arrow" | ".feather" | ".ipc":
            return pl.read_ipc(path)
        case ".jsonl" | ".ndjson":
            return pl.read_ndjson(path)
        case _:
            msg = f"unsupported table format '{suffix}' for {path}"
            raise ValueError(msg)


@dataclass(frozen=True)
class FeatureTable:
    """A cleaned table with one identifier column, numeric features and an optional target."""

    frame: pl.DataFrame
    id_column: str
    target_column: str | None = None

    @cached_property
    def features(self) -> list[str]:
        """Numeric feature columns, in table order."""
        reserved = {self.id_column, self.target_column}
        return [name for name, dtype in self.frame.schema.items() if name not in reserved and dtype.is_numeric()]

    @cached_property
    def annotations(self) -> list[str]:
        """Non-numeric columns that are carried along but not used as features."""
        reserved = {self.id_column, self.target_column, *self.features}
        return [name for name in self.frame.columns if name not in reserved]

    @property
    def ids(self) -> list[str]:
        """Row identifiers."""
        return self.frame.get_column(self.id_column).cast(pl.String).to_list()

    def matrix(self) -> np.ndarray:
        """Return the feature matrix as a float64 array (rows x features)."""
        return self.frame.select(self.features).cast(pl.Float64).to_numpy()

    @cached_property
    def task(self) -> str:
        """Learning task implied by the target: ``classification``, ``regression`` or ``none``."""
        if self.target_column is None:
            return "none"
        target = self.frame.get_column(self.target_column)
        if target.dtype.is_numeric() and target.n_unique() > MAX_CLASSES:
            return "regression"
        return "classification"

    def encoded_target(self) -> tuple[np.ndarray | None, list[str] | None]:
        """Encode the target for machine learning.

        Returns:
            ``(y, classes)``: class codes plus sorted class labels for classification,
            float values and ``None`` for regression, or ``(None, None)`` without a target.
        """
        if self.target_column is None:
            return None, None
        target = self.frame.get_column(self.target_column)
        if self.task == "regression":
            return target.cast(pl.Float64).to_numpy(), None
        labels = target.cast(pl.String)
        classes = sorted(labels.unique().to_list())
        codes = np.array([classes.index(label) for label in labels.to_list()], dtype=np.int64)
        return codes, classes


Writer = Callable[[FeatureTable, Path], None]


@dataclass(frozen=True)
class ExportFormat:
    """An export format: its file suffix and writer."""

    suffix: str
    writer: Writer


def _write_parquet(table: FeatureTable, path: Path) -> None:
    table.frame.write_parquet(path)


def _write_arrow(table: FeatureTable, path: Path) -> None:
    table.frame.write_ipc(path)


def _write_csv(table: FeatureTable, path: Path) -> None:
    table.frame.write_csv(path)


def _write_tsv(table: FeatureTable, path: Path) -> None:
    table.frame.write_csv(path, separator="\t")


def _write_jsonl(table: FeatureTable, path: Path) -> None:
    table.frame.write_ndjson(path)


def _write_npz(table: FeatureTable, path: Path) -> None:
    arrays = {
        "X": table.matrix(),
        "ids": np.asarray(table.ids, dtype=str),
        "feature_names": np.asarray(table.features, dtype=str),
    }
    y, classes = table.encoded_target()
    if y is not None:
        arrays["y"] = y
    if classes is not None:
        arrays["classes"] = np.asarray(classes, dtype=str)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, array in arrays.items():
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, np.ascontiguousarray(array), allow_pickle=False)
            entry = zipfile.ZipInfo(f"{name}.npy", date_time=_ZIP_EPOCH)
            archive.writestr(entry, buffer.getvalue(), compress_type=zipfile.ZIP_DEFLATED)


def _write_svmlight(table: FeatureTable, path: Path) -> None:
    y, _ = table.encoded_target()
    if y is None:
        y = np.zeros(table.frame.height)
    dump_svmlight_file(table.matrix(), y, str(path), zero_based=True)


EXPORT_FORMATS: dict[str, ExportFormat] = {
    "parquet": ExportFormat(".parquet", _write_parquet),
    "arrow": ExportFormat(".arrow", _write_arrow),
    "csv": ExportFormat(".csv", _write_csv),
    "tsv": ExportFormat(".tsv", _write_tsv),
    "jsonl": ExportFormat(".jsonl", _write_jsonl),
    "npz": ExportFormat(".npz", _write_npz),
    "svmlight": ExportFormat(".svmlight", _write_svmlight),
}


def parse_formats(text: str) -> list[str]:
    """Parse a comma-separated list of export formats.

    Args:
        text: For example ``"parquet,tsv,npz"``.

    Returns:
        Unique format names in the order given.

    Raises:
        ValueError: If the list is empty or names an unknown format.
    """
    formats = list(dict.fromkeys(part.strip().lower() for part in text.split(",") if part.strip()))
    if not formats:
        msg = "at least one export format is required"
        raise ValueError(msg)
    unknown = [name for name in formats if name not in EXPORT_FORMATS]
    if unknown:
        msg = f"unknown export format(s) {', '.join(unknown)}; choose from {', '.join(EXPORT_FORMATS)}"
        raise ValueError(msg)
    return formats


def export_table(table: FeatureTable, prefix: str, formats: Sequence[str], outdir: Path) -> list[Path]:
    """Write ``table`` once per format as ``<outdir>/<prefix><suffix>``.

    Args:
        table: Table to export.
        prefix: File name stem, usually the dataset id.
        formats: Names from :data:`EXPORT_FORMATS`.
        outdir: Existing output directory.

    Returns:
        The written paths, in the order of ``formats``.
    """
    paths = []
    for name in formats:
        export = EXPORT_FORMATS[name]
        path = outdir / f"{prefix}{export.suffix}"
        export.writer(table, path)
        paths.append(path)
    return paths
