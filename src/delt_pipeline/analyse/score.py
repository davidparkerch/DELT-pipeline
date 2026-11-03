"""Score every compound of a sub-library count table with all enrichment methods."""

from typing import Any

import numpy as np
import polars as pl

from .._numbers import significant
from ..library.schema import BB_COLUMNS, CONTROL, ID, LOAD, TARGET
from .enrichment import (
    BetaPrior,
    Exposure,
    count_difference,
    enrichment_ratio,
    faver_z,
    fit_beta_prior,
    poisson_ratio,
    shrunk_log2_ratio,
)

SCORES = ("count_diff", "er", "poisson_lower", "z_faver", "eb_log2_ratio")
"""Score columns, in the order they are reported. Higher means more enriched for all of them."""

SCORE_LABELS = {
    "count_diff": "Counts (mean target - mean control)",
    "er": "Enrichment ratio (pseudocount)",
    "poisson_lower": "Poisson rate ratio, lower bound",
    "z_faver": "Normalised z-score (Faver 2019)",
    "eb_log2_ratio": "Empirical-Bayes shrunken log2 ratio",
}

SIGNIFICANT = 6
"""Significant digits kept in score columns (byte-stable outputs)."""

TOP_COMPOUNDS = 10
"""Compounds listed in the QC summary."""


def exposure_from_qc(counts_qc: dict[str, Any]) -> Exposure:
    """Target and control exposure: summed whole-library depth of the replicates of each arm."""
    depth = counts_qc["library"]["depth"]
    return Exposure(
        target=float(sum(depth[name] for name in TARGET)), control=float(sum(depth[name] for name in CONTROL))
    )


def _pooled(counts: pl.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    c_t = counts.select(pl.sum_horizontal(TARGET)).to_series().to_numpy()
    c_c = counts.select(pl.sum_horizontal(CONTROL)).to_series().to_numpy()
    return c_t, c_c


def compute_scores(
    c_t: np.ndarray,
    c_c: np.ndarray,
    exposure: Exposure,
    prior: BetaPrior,
    *,
    diversity: int,
    alpha: float,
    level: float,
) -> dict[str, np.ndarray]:
    """All scores (and the rate-ratio interval) from pooled counts; also used for absent compounds."""
    replicates = len(TARGET)
    estimate, lower, upper = poisson_ratio(c_t, c_c, exposure, level)
    return {
        "count_diff": count_difference(np.asarray(c_t) / replicates, np.asarray(c_c) / replicates),
        "er": enrichment_ratio(c_t, c_c, exposure, alpha),
        "poisson_ratio": estimate,
        "poisson_lower": lower,
        "poisson_upper": upper,
        "z_faver": faver_z(c_t, exposure.target, diversity),
        "eb_log2_ratio": shrunk_log2_ratio(c_t, c_c, prior, exposure),
    }


def _pairwise_correlation(frame: pl.DataFrame) -> dict[str, float | None]:
    """Pearson correlation of log1p target counts for every pair of replicates."""
    logged = frame.select(pl.col(TARGET).log1p())
    pairs = {}
    for index, first in enumerate(TARGET):
        for second in TARGET[index + 1 :]:
            x, y = logged.get_column(first).to_numpy(), logged.get_column(second).to_numpy()
            constant = x.size < 2 or x.std() == 0 or y.std() == 0  # noqa: PLR2004 - correlation needs two points
            pairs[f"{first}~{second}"] = None if constant else significant(np.corrcoef(x, y)[0, 1])
    return pairs


def replicate_correlation(counts: pl.DataFrame) -> dict[str, dict[str, float | None]]:
    """How well target replicates agree, per compound and per disynthon (BB-A x BB-B).

    Single compounds have so few counts that Poisson noise dominates and replicates barely correlate.
    Summed over the third building block, the same data agree closely: the reason to read DEL
    results at the synthon level.
    """
    compounds = counts.filter(pl.sum_horizontal(TARGET) > 0)
    disynthons = counts.group_by("bb_a", "bb_b").agg(pl.col(TARGET).sum())
    return {"compound": _pairwise_correlation(compounds), "disynthon_ab": _pairwise_correlation(disynthons)}


def score_table(
    counts: pl.DataFrame,
    counts_qc: dict[str, Any],
    sublibrary: dict[str, Any],
    *,
    alpha: float = 0.5,
    level: float = 0.95,
) -> tuple[pl.DataFrame, dict[str, Any]]:
    """Score a compact count table.

    The Beta prior is fitted only on compounds whose BB-A was drawn at random: BB-As chosen because a
    held-out (validated) compound uses them are enriched for binders and would bias it.

    Args:
        counts: Compact count table (``library.subset.write_counts``).
        counts_qc: Its QC summary, holding the whole-library depth.
        sublibrary: The BB-A selection (``library.subset.select_sublibrary``).
        alpha: Pseudocount of the enrichment ratio.
        level: Confidence level of the rate-ratio interval.

    Returns:
        ``(scores, qc)``: the count table with score columns added, sorted by ``poisson_lower``
        (descending), and a JSON-serialisable summary.
    """
    exposure = exposure_from_qc(counts_qc)
    c_t, c_c = _pooled(counts)
    random_bb_a = counts.get_column("bb_a").is_in(sublibrary["random"]).to_numpy()
    fit_on = random_bb_a if random_bb_a.any() else np.ones(counts.height, dtype=bool)
    prior = fit_beta_prior(c_t[fit_on], c_c[fit_on])
    scores = compute_scores(c_t, c_c, exposure, prior, diversity=sublibrary["diversity"], alpha=alpha, level=level)
    frame = counts.with_columns(
        target_sum=pl.Series(c_t, dtype=pl.Int32),
        control_sum=pl.Series(c_c, dtype=pl.Int32),
        target_support=pl.sum_horizontal(pl.col(name) > 0 for name in TARGET).cast(pl.Int8),
        random_bb_a=pl.Series(random_bb_a),
        # Undefined (0/0) and infinite (no control counts) ratios are stored as null:
        # not every export format can represent NaN or infinity.
        **{
            name: pl.Series(np.where(np.isfinite(values), values, np.nan)).fill_nan(None).round_sig_figs(SIGNIFICANT)
            for name, values in scores.items()
        },
    )
    frame = frame.sort(["poisson_lower", "target_sum", ID], descending=[True, True, False], nulls_last=True)
    top = frame.head(TOP_COMPOUNDS).select(ID, *BB_COLUMNS, "target_sum", "control_sum", LOAD, *SCORES)
    qc = {
        "n_compounds": frame.height,
        "exposure": {"target": exposure.target, "control": exposure.control},
        "alpha": alpha,
        "level": level,
        "diversity": sublibrary["diversity"],
        "prior": {
            "a": significant(prior.a),
            "b": significant(prior.b),
            "mean": significant(prior.mean),
            "null_fraction": significant(exposure.null_fraction),
            "n_compounds": prior.n_compounds,
            "fitted_on": "random BB-A" if random_bb_a.any() else "all",
        },
        "replicate_correlation": replicate_correlation(counts),
        "fraction_with_target_counts": significant(float((c_t > 0).mean())) if counts.height else None,
        "top": top.to_dicts(),
    }
    return frame, qc


def absent_scores(score_qc: dict[str, Any]) -> dict[str, float | None]:
    """Scores of a compound with no counts at all (absent from KinDEL tables, which skip such rows)."""
    exposure = Exposure(score_qc["exposure"]["target"], score_qc["exposure"]["control"])
    prior = BetaPrior(score_qc["prior"]["a"], score_qc["prior"]["b"], score_qc["prior"]["n_compounds"])
    zero = np.zeros(1)
    scores = compute_scores(
        zero, zero, exposure, prior, diversity=score_qc["diversity"], alpha=score_qc["alpha"], level=score_qc["level"]
    )
    return {name: significant(scores[name][0]) for name in SCORES}
