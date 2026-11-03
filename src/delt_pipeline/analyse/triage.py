"""Hit triage: put each target screen's top compounds into tiers, with flags that explain why.

The evidence, in the order a DEL chemist reads it:

=====================  ======================================================================
Flag                   Meaning
=====================  ======================================================================
``replicate_support``  counts in at least two of three target replicates
``disynthon_support``  one of its disynthons (AB, AC or BC) is enriched without the compound itself
``counter_binder``     enriched in a counter-screen (an unrelated protein, e.g. BCA)
``matrix_binder``      enriched in the no-protein (bead) control over its synthesis load
``shared_target``      also a top hit of another *target* screen
=====================  ======================================================================

Only ``counter_binder`` and ``matrix_binder`` mark artefacts. A compound shared by two kinase
targets (``shared_target``) is a plausible hinge binder, not a promiscuous one, so it is reported but
not penalised.

Tiers: ``1`` = replicate and disynthon support, ``2`` = one of them, ``3`` = neither, and
``flagged`` = counter or matrix binder, whatever its support.

The thresholds are lower confidence bounds of rate ratios, and each asks a different question:

* **counter-screen**: any significant enrichment (bound > 1) is a warning;
* **disynthon** (other members only): at least 1.5-fold, because family ratios are overdispersed
  (synthesis and bead effects are shared by whole families) and a bound above 1 is common even for
  non-binders;
* **bead control**: at least fourfold, because matrix counts and load are two noisy measurements of
  the same synthesis yield and their ratio scatters widely.
"""

from dataclasses import dataclass, replace
from typing import Any

import polars as pl

from .._numbers import significant
from ..library.schema import CONTROL, ID, LOAD
from .enrichment import Exposure, poisson_ratio
from .score import SCORES, exposure_from_qc, finite_series
from .synthons import DISYNTHON_LEVELS, KEYS, LEVELS

MIN_REPLICATES = 2
"""Target replicates with counts needed for ``replicate_support``."""

DISYNTHON_LOWER = 1.5
"""Lower bound of a disynthon's rate ratio needed for ``disynthon_support``."""

DISYNTHON_MIN_COMPOUNDS = 3
"""Other members with target counts a disynthon needs (one bright compound is not a family)."""

COUNTER_LOWER = 1.0
"""Lower bound of the rate ratio in a counter-screen that makes a ``counter_binder``."""

MATRIX_LOWER = 4.0
"""Lower bound of the matrix-over-load rate ratio that makes a ``matrix_binder``."""

TIERS = ("1", "2", "3", "flagged")

HIT_ID = "hit_id"
"""Unique row id of the hits table (``<screen>:<molecule_hash>``): a compound can be a hit of several screens."""

FLAGS = ("replicate_support", "disynthon_support", "counter_binder", "matrix_binder", "shared_target")


@dataclass(frozen=True)
class ScreenInput:
    """What triage needs from one scored screen."""

    id: str
    role: str
    """``target`` or ``counter``."""
    top: pl.DataFrame
    """The screen's top compounds (all score columns)."""
    scores: pl.LazyFrame
    """Every compound of the sub-library, scanned lazily (only looked up by ``molecule_hash``)."""
    synthons: pl.DataFrame
    """From :func:`~delt_pipeline.analyse.synthons.aggregate_synthons`."""
    counts_qc: dict[str, Any]
    """Whole-library depth of every count column."""


def matrix_exposure(counts_qc: dict[str, Any]) -> Exposure:
    """Exposure of the matrix (bead) control against the pre-selection load, over the whole library."""
    depth = counts_qc["library"]["depth"]
    return Exposure(target=float(sum(depth[name] for name in CONTROL)), control=float(depth[LOAD]))


def _best_disynthon(top: pl.DataFrame, synthons: pl.DataFrame, exposure: Exposure, level: float) -> pl.DataFrame:
    """Best-supported disynthon of each compound, *leaving the compound itself out*.

    A family supports a compound only if its other members are enriched: otherwise a single bright
    compound would make its own disynthons look enriched. The compound's counts are subtracted from
    each of its three disynthons, which are then re-scored with the exposure-mapped Poisson bound.

    Returns:
        ``molecule_hash``, ``disynthon`` (``AB``, ``AC`` or ``BC``) and ``disynthon_lower``; null when no
        disynthon has enough other members with counts.
    """
    own = top.select(ID, *KEYS, own_target="target_sum", own_control="control_sum")
    candidates = []
    for name in DISYNTHON_LEVELS:
        keys = list(LEVELS[name])
        pairs = synthons.filter(pl.col("level") == name).select(*keys, "n_with_target", "target_sum", "control_sum")
        candidates.append(
            own.join(pairs, on=keys, how="inner").select(
                ID,
                pl.lit(name.removeprefix("di_").upper()).alias("disynthon"),
                (pl.col("n_with_target") - (pl.col("own_target") > 0).cast(pl.Int32)).alias("others"),
                (pl.col("target_sum") - pl.col("own_target")).alias("rest_target"),
                (pl.col("control_sum") - pl.col("own_control")).alias("rest_control"),
            )
        )
    pairs = pl.concat(candidates).filter(pl.col("others") >= DISYNTHON_MIN_COMPOUNDS)
    _, lower, _ = poisson_ratio(
        pairs.get_column("rest_target").to_numpy(), pairs.get_column("rest_control").to_numpy(), exposure, level
    )
    best = (
        pairs.with_columns(disynthon_lower=finite_series(lower))
        .sort(ID, "disynthon_lower", "disynthon", descending=[False, True, False], nulls_last=True)
        .unique(ID, keep="first", maintain_order=True)
        .select(ID, "disynthon", "disynthon_lower")
    )
    return top.select(ID).join(best, on=ID, how="left")


def _counter_lower(hashes: pl.Series, counters: list[ScreenInput]) -> pl.DataFrame:
    """Highest Poisson lower bound of each compound over all counter-screens (null if never observed)."""
    frame = pl.DataFrame({ID: hashes.unique().sort()})
    if not counters:
        return frame.with_columns(counter_lower=pl.lit(None, dtype=pl.Float64))
    found = pl.concat(
        [
            screen.scores.filter(pl.col(ID).is_in(hashes.implode())).select(ID, "poisson_lower").collect()
            for screen in counters
        ]
    )
    best = found.group_by(ID).agg(counter_lower=pl.col("poisson_lower").max())
    return frame.join(best, on=ID, how="left")


def _tier() -> pl.Expr:
    support = pl.col("replicate_support").cast(pl.Int8) + pl.col("disynthon_support").cast(pl.Int8)
    return (
        pl.when(pl.col("counter_binder") | pl.col("matrix_binder"))
        .then(pl.lit("flagged"))
        .when(support == 2)  # noqa: PLR2004 - both kinds of support
        .then(pl.lit("1"))
        .when(support == 1)
        .then(pl.lit("2"))
        .otherwise(pl.lit("3"))
        .alias("tier")
    )


def _reasons() -> pl.Expr:
    """Human-readable evidence, e.g. ``3/3 replicates; disynthon AB (LB 4.2); also a top hit of ddr1``."""
    parts = [
        pl.format("{}/3 replicates", pl.col("target_support")),
        pl.when(pl.col("disynthon_support")).then(
            pl.format("disynthon {} (LB {})", pl.col("disynthon"), pl.col("disynthon_lower").round(1))
        ),
        pl.when(pl.col("counter_binder")).then(
            pl.format("counter-screen binder (LB {})", pl.col("counter_lower").round(1))
        ),
        pl.when(pl.col("matrix_binder")).then(pl.format("bead binder (LB {})", pl.col("matrix_lower").round(1))),
        pl.when(pl.col("shared_target")).then(pl.format("also a top hit of {}", pl.col("shared_with"))),
    ]
    return pl.concat_str(parts, separator="; ", ignore_nulls=True).alias("reasons")


def triage(screens: list[ScreenInput], *, level: float = 0.95) -> pl.DataFrame:
    """Tier and flag the top compounds of every target screen.

    Args:
        screens: Every scored screen; ``target`` screens are triaged, ``counter`` screens are only
            looked up.
        level: Confidence level of every rate-ratio bound (disynthon and bead).

    Returns:
        One row per (target screen, compound), keyed by ``hit_id`` and sorted by screen and rank. All score
        columns of the screen's top table are kept, followed by the evidence columns, the flags, ``tier``
        and ``reasons``.
    """
    targets = sorted((screen for screen in screens if screen.role == "target"), key=lambda screen: screen.id)
    counters = sorted((screen for screen in screens if screen.role == "counter"), key=lambda screen: screen.id)
    tops = {screen.id: set(screen.top.get_column(ID).to_list()) for screen in targets}
    frames = []
    for screen in targets:
        top = screen.top.with_row_index("rank", offset=1).with_columns(pl.col("rank").cast(pl.Int32))
        counter = _counter_lower(top.get_column(ID), counters)
        _, matrix_lb, _ = poisson_ratio(
            top.get_column("control_sum").to_numpy(),
            top.get_column(LOAD).to_numpy(),
            matrix_exposure(screen.counts_qc),
            level,
        )
        shared_with = [
            ";".join(other for other in sorted(tops) if other != screen.id and molecule in tops[other])
            for molecule in top.get_column(ID).to_list()
        ]
        frame = (
            top.with_columns(
                pl.lit(screen.id).alias("screen"),
                matrix_lower=finite_series(matrix_lb),
                shared_with=pl.Series(shared_with, dtype=pl.String),
            )
            .join(
                _best_disynthon(top, screen.synthons, exposure_from_qc(screen.counts_qc), level),
                on=ID,
                how="left",
                maintain_order="left",
            )
            .join(counter, on=ID, how="left", maintain_order="left")
            .with_columns(
                replicate_support=pl.col("target_support") >= MIN_REPLICATES,
                disynthon_support=(pl.col("disynthon_lower") > DISYNTHON_LOWER).fill_null(value=False),
                counter_binder=(pl.col("counter_lower") > COUNTER_LOWER).fill_null(value=False),
                matrix_binder=(pl.col("matrix_lower") > MATRIX_LOWER).fill_null(value=False),
                shared_target=pl.col("shared_with") != "",
                shared_with=pl.when(pl.col("shared_with") != "").then(pl.col("shared_with")),
            )
        )
        frames.append(frame.with_columns(_tier(), _reasons()))
    if not frames:
        msg = "triage needs at least one screen with role 'target'"
        raise ValueError(msg)
    hits = pl.concat(frames, how="vertical").with_columns(hit_id=pl.concat_str("screen", pl.lit(":"), ID))
    leading = [HIT_ID, "screen", "rank", "tier", ID]
    return hits.select(*leading, *(name for name in hits.columns if name not in leading))


def _tier_counts(hits: pl.DataFrame) -> dict[str, int]:
    counts = dict(hits.get_column("tier").value_counts().iter_rows())
    return {tier: int(counts.get(tier, 0)) for tier in TIERS}


def top_by(scores: pl.LazyFrame, score: str, n: int) -> pl.DataFrame:
    """The ``n`` best compounds by one score (ties broken by target counts, then id), all columns.

    Only the ranking columns are sorted, then the chosen rows are read back, so millions of
    compounds are never sorted with every column.
    """
    chosen = (
        scores.select(ID, score, "target_sum")
        .sort(score, "target_sum", ID, descending=[True, True, False], nulls_last=True)
        .head(n)
        .with_row_index("order")
        .collect()
    )
    return (
        scores.filter(pl.col(ID).is_in(chosen.get_column(ID).implode()))
        .collect()
        .join(chosen.select(ID, "order"), on=ID)
        .sort("order")
        .drop("order")
    )


def tiers_by_score(screens: list[ScreenInput], n: int, *, level: float = 0.95) -> dict[str, dict[str, dict[str, int]]]:
    """Tier counts of the top ``n`` compounds of each target screen, ranked by each score in turn.

    The same triage, applied to different rankings, shows what each score lets through. A target-only
    score (Faver z) has no control, so the bead binders it promotes all end up ``flagged``.
    """
    result = {}
    for score in SCORES:
        reranked = [
            replace(screen, top=top_by(screen.scores, score, n)) if screen.role == "target" else screen
            for screen in screens
        ]
        hits = triage(reranked, level=level)
        result[score] = {
            str(screen): _tier_counts(group) for (screen,), group in hits.group_by("screen", maintain_order=True)
        }
    return result


def triage_qc(
    hits: pl.DataFrame, *, level: float = 0.95, by_score: dict[str, dict[str, dict[str, int]]] | None = None
) -> dict[str, Any]:
    """Tier and flag counts per target screen, the thresholds used and (optionally) tiers by score."""
    screens = {}
    for (screen,), group in hits.group_by("screen", maintain_order=True):
        screens[str(screen)] = {
            "n_compounds": group.height,
            "tiers": _tier_counts(group),
            "flags": {flag: int(group.get_column(flag).sum()) for flag in FLAGS},
            "top_tier1": group.filter(pl.col("tier") == "1")
            .head(20)
            .select("rank", ID, "bb_a", "bb_b", "bb_c", "target_sum", "control_sum", "poisson_lower", "reasons")
            .to_dicts(),
        }
    thresholds = {
        "min_replicates": MIN_REPLICATES,
        "disynthon_lower": significant(DISYNTHON_LOWER),
        "disynthon_min_compounds": DISYNTHON_MIN_COMPOUNDS,
        "counter_lower": significant(COUNTER_LOWER),
        "matrix_lower": significant(MATRIX_LOWER),
        "level": level,
    }
    return {"screens": screens, "thresholds": thresholds, "n_hits": hits.height, "by_score": by_score or {}}
