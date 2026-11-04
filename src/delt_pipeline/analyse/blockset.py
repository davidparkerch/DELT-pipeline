"""Building-block set enrichment: is a block's set of compounds ranked above the rest?

This is the gene-set-enrichment idea carried over to DEL data. A building block is a "gene set": the
compounds that contain it. Ranking all compounds by a score, a block that drives binding has its
compounds concentrated at the top. The test is the Mann-Whitney U test of the block's ranks against all
other ranks, computed from rank sums so that thousands of blocks over millions of compounds need one
ranking and one group-by. The p-values are corrected for the number of blocks tested (Benjamini-Hochberg).

Two cautions, stated in the docs: compound sets overlap heavily (every compound is in three sets), so the
tests are not independent and the q-values are indicative; and the universe is the compounds that were
sequenced at all, because compounds absent from every column are mostly failed syntheses.
"""

import math
from typing import Any

import numpy as np
import polars as pl
from scipy import stats

from .._numbers import significant
from ..library.schema import POSITIONS

SCORE = "poisson_lower"
"""The score the compounds are ranked by."""
Q_LEVEL = 0.05
TOP_BLOCKS = 10
"""Blocks listed in the QC summary."""


def blockset_enrichment(scores: pl.LazyFrame, top: pl.DataFrame | None = None, score: str = SCORE) -> pl.DataFrame:
    """Test every building block of the scored compounds.

    Args:
        scores: Scored compounds with ``bb_a``, ``bb_b``, ``bb_c`` and ``score`` (``null`` counts as 0).
        top: The screen's top compounds, to count how many of them contain each block (``n_top``).
        score: Column to rank by (higher is better).

    Returns:
        One row per block with at least one but not all compounds: ``position``, ``bb_id``,
        ``n_compounds``, ``n_top``, ``auc`` (probability that a compound of the block outranks one that is
        not, ties count half), ``z``, one-sided ``p_value``, ``log10_p``, ``q_value`` (Benjamini-Hochberg over
        all blocks) and the block's ``median_score``. Sorted by ``q_value``, then ``auc`` (descending).
    """
    frame = (
        scores.select(*(f"bb_{position}" for position in POSITIONS), pl.col(score).fill_null(0.0).alias("value"))
        .collect()
        .with_columns(rank=pl.col("value").rank("average"))
    )
    total = frame.height
    ties = frame.get_column("value").value_counts().get_column("count").cast(pl.Float64)
    tie_term = float((ties**3 - ties).sum()) / (total * (total - 1)) if total > 1 else 0.0
    rows = []
    for position in POSITIONS:
        column = f"bb_{position}"
        groups = frame.group_by(column).agg(
            pl.len().alias("n"), pl.col("rank").sum().alias("rank_sum"), pl.col("value").median().alias("median")
        )
        in_top = (
            top.group_by(column)
            .len()
            .rename({column: "bb_id", "len": "n_top"})
            .with_columns(pl.col("bb_id").cast(pl.UInt16))
            if top is not None
            else pl.DataFrame(schema={"bb_id": pl.UInt16, "n_top": pl.UInt32})
        )
        rows.append(
            groups.select(
                pl.lit(position).alias("position"),
                pl.col(column).cast(pl.UInt16).alias("bb_id"),
                pl.col("n").alias("n_compounds"),
                "rank_sum",
                pl.col("median").alias("median_score"),
            ).join(in_top, on="bb_id", how="left")
        )
    table = pl.concat(rows).filter((pl.col("n_compounds") > 0) & (pl.col("n_compounds") < total))
    n = table.get_column("n_compounds").to_numpy().astype(float)
    rest = total - n
    u = table.get_column("rank_sum").to_numpy() - n * (n + 1) / 2
    variance = n * rest / 12 * ((total + 1) - tie_term)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(variance > 0, (u - n * rest / 2 - 0.5) / np.sqrt(variance), 0.0)
        auc = u / (n * rest)
    p_value = stats.norm.sf(z)
    q_value = stats.false_discovery_control(p_value, method="bh") if p_value.size else p_value
    result = table.drop("rank_sum").with_columns(
        n_top=pl.col("n_top").fill_null(0),
        auc=pl.Series(auc),
        z=pl.Series(z),
        p_value=pl.Series(p_value),
        log10_p=pl.Series(stats.norm.logsf(z) / math.log(10)),
        q_value=pl.Series(q_value),
    )
    result = result.select(
        "position", "bb_id", "n_compounds", "n_top", "auc", "z", "p_value", "log10_p", "q_value", "median_score"
    )
    return result.with_columns(
        pl.col("auc", "z", "p_value", "log10_p", "q_value", "median_score").round_sig_figs(6)
    ).sort(["q_value", "auc", "position", "bb_id"], descending=[False, True, False, False])


def blockset_qc(table: pl.DataFrame, n_compounds: int, score: str = SCORE) -> dict[str, Any]:
    """Summary of a set-enrichment table: how many blocks are significant, and the strongest ones."""
    significant_blocks = table.filter(pl.col("q_value") < Q_LEVEL)
    return {
        "score": score,
        "n_compounds": n_compounds,
        "n_tests": table.height,
        "q_level": Q_LEVEL,
        "n_significant": {
            position: int(significant_blocks.filter(pl.col("position") == position).height) for position in POSITIONS
        },
        "top": [
            {key: significant(value) if isinstance(value, float) else value for key, value in row.items()}
            for row in table.head(TOP_BLOCKS).to_dicts()
        ],
    }
