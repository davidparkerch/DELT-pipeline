"""Validate enrichment scores against measured binding affinities (KD) of resynthesised compounds.

KinDEL ships KD values for 30-50 compounds per target, measured on-DNA (fluorescence polarisation) and
off-DNA (SPR). A good score ranks tight binders (low KD) first, so we report the **negative** Spearman
correlation between score and KD, as in the KinDEL paper (higher is better).
"""

from typing import Any

import numpy as np
import polars as pl
from scipy import stats

from .._numbers import significant
from ..library.schema import ID
from .score import SCORE_LABELS, SCORES

PUBLISHED: dict[tuple[str, str], dict[str, float]] = {
    ("mapk14", "on"): {"count_diff": 0.778, "poisson": 0.737},
    ("mapk14", "off"): {"count_diff": 0.353, "poisson": 0.166},
    ("ddr1", "on"): {"count_diff": 0.695, "poisson": 0.779},
    ("ddr1", "off"): {"count_diff": 0.355, "poisson": 0.441},
}
"""KinDEL baselines (Chen et al., ICML 2025, Tables 1-2: in-library, -Spearman vs KD)."""

MIN_COMPOUNDS = 3


def read_heldout(frame: pl.DataFrame) -> pl.DataFrame:
    """In-library held-out measurements (rows with a ``molecule_hash``) and their KD.

    Every measurement is kept, including repeated measurements of one compound (two in DDR1 off-DNA),
    exactly as the KinDEL benchmark does, so the correlations are comparable with the published ones.
    """
    return frame.select(pl.col(ID).cast(pl.String), pl.col("kd").cast(pl.Float64)).drop_nulls()


def _negative_spearman(score: np.ndarray, kd: np.ndarray) -> float | None:
    if score.size < MIN_COMPOUNDS or np.ptp(score) == 0 or np.ptp(kd) == 0:
        return None
    return significant(-stats.spearmanr(score, kd).statistic, 3)


def validate_scores(
    scores: pl.DataFrame, heldout: pl.DataFrame, absent: dict[str, float | None], target: str, condition: str
) -> dict[str, Any]:
    """Correlate every score with KD for one target and assay condition.

    Held-out compounds absent from the screen table had no counts in any sample; they get the scores
    of an all-zero compound (``absent``) rather than being dropped, which would bias the comparison.

    Args:
        scores: Scored sub-library (must contain every in-library held-out compound that has counts).
        heldout: From :func:`read_heldout`.
        absent: Scores of a compound with no counts (``analyse.score.absent_scores``).
        target: Screen target, e.g. ``mapk14`` (selects the published baselines).
        condition: ``on`` or ``off`` (DNA).

    Returns:
        One row per score with the correlation and, where KinDEL published one, its baseline.
    """
    joined = heldout.join(scores.select(ID, *SCORES), on=ID, how="left")
    n_absent = joined.filter(pl.col(SCORES[0]).is_null()).height
    joined = joined.with_columns(pl.col(name).fill_null(absent[name]) for name in SCORES)
    published = PUBLISHED.get((target.lower(), condition), {})
    rows = []
    for name in SCORES:
        score = joined.get_column(name).cast(pl.Float64).to_numpy()
        reference = published.get("poisson" if name == "poisson_lower" else name)
        rows.append(
            {
                "score": name,
                "label": SCORE_LABELS[name],
                "spearman": _negative_spearman(score, joined.get_column("kd").to_numpy()),
                "published": reference,
            }
        )
    return {
        "target": target,
        "condition": condition,
        "n_compounds": joined.height,
        "n_absent": n_absent,
        "scores": rows,
    }
