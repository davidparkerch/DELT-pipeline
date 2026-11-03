"""Cleaning raw feature tables with polars."""

import re
from dataclasses import dataclass
from typing import Any

import polars as pl

from ..io import FeatureTable

IMPUTE_STRATEGIES = ("median", "mean", "drop")
DECIMALS = 6
"""Features are rounded so results are stable across platforms."""


def normalise_name(name: str) -> str:
    """Make a column name safe for every export format (letters, digits, underscores).

    Case is preserved so gene symbols and descriptor names stay recognisable.
    """
    cleaned = re.sub(r"[^0-9A-Za-z_]+", "_", name.strip()).strip("_")
    return cleaned or "column"


def normalise_columns(frame: pl.DataFrame) -> pl.DataFrame:
    """Rename all columns with :func:`normalise_name`, de-duplicating clashes with a numeric suffix."""
    seen: dict[str, int] = {}
    names = []
    for column in frame.columns:
        name = normalise_name(column)
        count = seen.get(name, 0)
        seen[name] = count + 1
        names.append(name if count == 0 else f"{name}_{count + 1}")
    return frame.rename(dict(zip(frame.columns, names, strict=True)))


@dataclass(frozen=True)
class WrangleResult:
    """A cleaned table and a JSON-serialisable QC summary."""

    table: FeatureTable
    qc: dict[str, Any]


def _require_columns(frame: pl.DataFrame, *columns: str | None) -> None:
    missing = [column for column in columns if column is not None and column not in frame.columns]
    if missing:
        msg = f"column(s) not found: {', '.join(missing)}; available: {', '.join(frame.columns)}"
        raise ValueError(msg)


def wrangle(
    frame: pl.DataFrame,
    id_column: str,
    target_column: str | None = None,
    impute: str = "median",
    max_missing: float = 0.5,
) -> WrangleResult:
    """Clean a raw feature table.

    Steps: normalise column names, drop rows without an identifier (or target), reject
    duplicate identifiers, drop numeric features that are mostly missing or constant,
    impute the remaining gaps, round features and sort rows by identifier.

    Args:
        frame: Raw table.
        id_column: Column with unique row identifiers (sample, compound, ...).
        target_column: Optional label/response column kept for supervised analysis.
        impute: ``median``, ``mean`` or ``drop`` (drop rows with missing features).
        max_missing: Drop features with a larger fraction of missing values.

    Returns:
        The cleaned :class:`FeatureTable` and its QC summary.

    Raises:
        ValueError: On unknown strategies, missing columns, duplicate ids or no usable features.
    """
    if impute not in IMPUTE_STRATEGIES:
        msg = f"unknown impute strategy '{impute}'; choose from {', '.join(IMPUTE_STRATEGIES)}"
        raise ValueError(msg)
    frame = normalise_columns(frame)
    id_column = normalise_name(id_column)
    target_column = normalise_name(target_column) if target_column else None
    _require_columns(frame, id_column, target_column)

    rows_in = frame.height
    frame = frame.filter(pl.col(id_column).is_not_null()).with_columns(pl.col(id_column).cast(pl.String))
    if target_column is not None:
        frame = frame.filter(pl.col(target_column).is_not_null())
    if frame.height == 0:
        msg = "no rows left after removing rows without an identifier or target"
        raise ValueError(msg)
    duplicated = sorted(frame.filter(pl.col(id_column).is_duplicated()).get_column(id_column).unique().to_list())
    if duplicated:
        msg = f"duplicate identifiers in '{id_column}': {', '.join(duplicated)}"
        raise ValueError(msg)

    # NaN and null both mean "missing" from here on.
    floats = [name for name, dtype in frame.schema.items() if dtype.is_float()]
    frame = frame.with_columns(pl.col(name).fill_nan(None) for name in floats)

    candidates = FeatureTable(frame, id_column, target_column).features
    dropped: dict[str, str] = {}
    for name in candidates:
        column = frame.get_column(name)
        if column.null_count() / frame.height > max_missing:
            dropped[name] = "too_many_missing"
        elif column.drop_nulls().n_unique() <= 1:
            dropped[name] = "constant"
    features = [name for name in candidates if name not in dropped]
    if not features:
        msg = "no usable numeric feature columns"
        raise ValueError(msg)
    missing = {name: frame.get_column(name).null_count() for name in features}

    rows_before_impute = frame.height
    if impute == "drop":
        frame = frame.drop_nulls(subset=features)
    else:
        fill = pl.col(features).median() if impute == "median" else pl.col(features).mean()
        frame = frame.with_columns(pl.col(features).fill_null(fill))
    if frame.height == 0:
        msg = "no rows left after dropping rows with missing features"
        raise ValueError(msg)

    annotations = FeatureTable(frame, id_column, target_column).annotations
    targets = [target_column] if target_column is not None else []
    frame = (
        frame.with_columns(pl.col(features).cast(pl.Float64).round(DECIMALS))
        .select([id_column, *annotations, *features, *targets])
        .sort(id_column)
    )
    table = FeatureTable(frame, id_column, target_column)
    qc = {
        "id_column": id_column,
        "target_column": target_column,
        "task": table.task,
        "impute": impute,
        "max_missing": max_missing,
        "rows_in": rows_in,
        "rows_out": frame.height,
        "rows_dropped_missing_features": rows_before_impute - frame.height,
        "features": features,
        "annotations": annotations,
        "dropped_features": dict(sorted(dropped.items())),
        "missing_values": missing,
    }
    return WrangleResult(table, qc)
