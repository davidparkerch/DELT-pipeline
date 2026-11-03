"""Synthon aggregation: sum compound counts per building block (monosynthon) and per pair (disynthon).

A single DEL compound has a handful of counts, so its own ratio is noisy. A disynthon (two building
blocks, any third one) pools the counts of up to hundreds of compounds: if the pair binds, every member
is enriched, and the pooled ratio is far better determined. This is the "sharing of strength" that
DEL data offer instead of DESeq2's gene-wise priors, and the lines in the cube plot.

Pooled counts are scored with the same exposure-mapped Poisson rate ratio as single compounds.
"""

from typing import Any

import numpy as np
import polars as pl

from .._numbers import significant
from ..library.schema import CONTROL, LOAD, TARGET
from .enrichment import Exposure, enrichment_ratio, poisson_ratio
from .score import finite_series

LEVELS: dict[str, tuple[str, ...]] = {
    "mono_a": ("bb_a",),
    "mono_b": ("bb_b",),
    "mono_c": ("bb_c",),
    "di_ab": ("bb_a", "bb_b"),
    "di_ac": ("bb_a", "bb_c"),
    "di_bc": ("bb_b", "bb_c"),
}
"""Aggregation levels and the building-block columns that define a group."""

DISYNTHON_LEVELS = ("di_ab", "di_ac", "di_bc")

KEYS = ("bb_a", "bb_b", "bb_c")

TOP_SYNTHONS = 10
"""Disynthons listed per level in the QC summary."""


def aggregate_synthons(
    scores: pl.LazyFrame, exposure: Exposure, *, alpha: float = 0.5, level: float = 0.95
) -> pl.DataFrame:
    """Pool counts per monosynthon and disynthon, and score each group.

    Args:
        scores: Compound table with ``bb_a/b/c``, replicate counts and ``load`` (lazily scanned).
        exposure: Whole-library depth of the target and control arms (as for compounds).
        alpha: Pseudocount of the enrichment ratio.
        level: Confidence level of the rate-ratio interval.

    Returns:
        One row per synthon. ``level`` names the grouping and the building-block columns that are
        not part of it are null. Sorted by level, then Poisson lower bound (descending).
    """
    frames = []
    for name, keys in LEVELS.items():
        grouped = (
            scores.group_by(keys)
            .agg(
                pl.len().cast(pl.Int32).alias("n_compounds"),
                (pl.sum_horizontal(TARGET) > 0).sum().cast(pl.Int32).alias("n_with_target"),
                *(pl.col(column).sum().cast(pl.Int64) for column in (*TARGET, *CONTROL, LOAD)),
            )
            .collect(engine="streaming")
        )
        frames.append(
            grouped.with_columns(
                level=pl.lit(name),
                **{key: pl.lit(None, dtype=pl.UInt16) for key in KEYS if key not in keys},
            )
        )
    frame = pl.concat(frames, how="diagonal_relaxed").with_columns(
        target_sum=pl.sum_horizontal(TARGET), control_sum=pl.sum_horizontal(CONTROL)
    )
    c_t = frame.get_column("target_sum").to_numpy()
    c_c = frame.get_column("control_sum").to_numpy()
    estimate, lower, upper = poisson_ratio(c_t, c_c, exposure, level)
    frame = frame.with_columns(
        er=finite_series(enrichment_ratio(c_t, c_c, exposure, alpha)),
        poisson_ratio=finite_series(estimate),
        poisson_lower=finite_series(lower),
        poisson_upper=finite_series(upper),
    )
    order = list(LEVELS)
    return (
        frame.with_columns(pl.col(*KEYS).cast(pl.UInt16))
        .sort(
            pl.col("level").replace_strict(order, list(range(len(order))), return_dtype=pl.Int8),
            pl.col("poisson_lower"),
            "target_sum",
            *KEYS,
            descending=[False, True, True, False, False, False],
            nulls_last=True,
        )
        .select(
            "level",
            *KEYS,
            "n_compounds",
            "n_with_target",
            *TARGET,
            *CONTROL,
            LOAD,
            "target_sum",
            "control_sum",
            "er",
            "poisson_ratio",
            "poisson_lower",
            "poisson_upper",
        )
    )


def _correlation(frame: pl.DataFrame) -> float | None:
    """Mean pairwise Pearson correlation of log1p target counts across replicates."""
    logged = frame.select(pl.col(TARGET).log1p()).to_numpy()
    if logged.shape[0] < 2 or (logged.std(axis=0) == 0).any():  # noqa: PLR2004 - correlation needs two points
        return None
    matrix = np.corrcoef(logged, rowvar=False)
    pairs = [matrix[i, j] for i in range(len(TARGET)) for j in range(i + 1, len(TARGET))]
    return significant(sum(pairs) / len(pairs))


def synthon_qc(synthons: pl.DataFrame) -> dict[str, Any]:
    """Size and replicate agreement of every level, and the top disynthons of each kind."""
    levels = {}
    for name in LEVELS:
        subset = synthons.filter(pl.col("level") == name)
        levels[name] = {
            "n_synthons": subset.height,
            "replicate_correlation": _correlation(subset.filter(pl.col("target_sum") > 0)),
            "n_enriched": int(subset.filter(pl.col("poisson_lower") > 1).height),
        }
    top = {
        name: synthons.filter(pl.col("level") == name)
        .head(TOP_SYNTHONS)
        .select(*LEVELS[name], "n_compounds", "n_with_target", "target_sum", "control_sum", "poisson_lower")
        .to_dicts()
        for name in DISYNTHON_LEVELS
    }
    return {"levels": levels, "top": top}
