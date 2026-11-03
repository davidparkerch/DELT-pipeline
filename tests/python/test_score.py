from pathlib import Path

import polars as pl
import pytest

from delt_pipeline.analyse import SCORES, absent_scores, score_table
from delt_pipeline.library import SimulatedLibrary


def compound(frame: pl.DataFrame, a: int, b: int, c: int | None = None) -> pl.DataFrame:
    selected = frame.filter((pl.col("bb_a") == a) & (pl.col("bb_b") == b))
    return selected if c is None else selected.filter(pl.col("bb_c") == c)


def rank(frame: pl.DataFrame, score: str) -> pl.DataFrame:
    return frame.sort(score, descending=True, nulls_last=True).with_row_index("rank")


def test_planted_signal(mapk14_scores: tuple[pl.DataFrame, dict], library: SimulatedLibrary) -> None:
    """Naive ER is fooled by a low-count outlier; the Poisson lower bound and EB shrinkage are not."""
    scores, _ = mapk14_scores
    planted = library.planted
    outlier = compound(scores, *planted.outlier)
    assert outlier.select("target_sum", "control_sum", "target_support").row(0) == (6, 0, 1)

    # The binder family member with the most evidence (most target counts).
    family = compound(scores, *planted.binder).sort("target_sum", descending=True)
    best = family.row(0, named=True)
    assert outlier.item(0, "er") > best["er"]
    assert outlier.item(0, "poisson_lower") < best["poisson_lower"]
    assert outlier.item(0, "eb_log2_ratio") < best["eb_log2_ratio"]

    by_er = rank(scores, "er")
    by_lower = rank(scores, "poisson_lower")
    assert compound(by_er, *planted.outlier).item(0, "rank") < 20
    assert compound(by_lower, *planted.outlier).item(0, "rank") > 30
    # Planted families (the MAPK14 binder and the kinase-shared one) own the top of the ranking.
    families = {planted.binder, planted.shared}
    top = by_lower.head(5).select("bb_a", "bb_b").rows()
    assert all(pair in families for pair in top)


def test_scores_table(mapk14_scores: tuple[pl.DataFrame, dict], mapk14_counts: tuple[Path, dict]) -> None:
    scores, qc = mapk14_scores
    assert scores.height == pl.read_parquet(mapk14_counts[0]).height
    assert set(SCORES) <= set(scores.columns)
    assert scores.equals(
        scores.sort("poisson_lower", "target_sum", descending=True, nulls_last=True, maintain_order=True)
    )
    assert scores.get_column("target_support").max() == 3
    # Undefined ratios are null rather than NaN or infinity.
    assert scores.get_column("poisson_ratio").null_count() > 0
    assert scores.get_column("poisson_ratio").is_nan().sum() == 0
    assert qc["prior"]["fitted_on"] == "random BB-A"
    assert qc["prior"]["n_compounds"] < scores.filter(pl.col("target_sum") + pl.col("control_sum") > 0).height
    assert len(qc["top"]) == 10
    correlation = qc["replicate_correlation"]
    assert set(correlation) == {"compound", "disynthon_ab"}
    assert set(correlation["compound"]) == {"target_1~target_2", "target_1~target_3", "target_2~target_3"}
    # Replicates agree far better on disynthons than on single compounds.
    assert min(correlation["disynthon_ab"].values()) > max(correlation["compound"].values())


def test_prior_falls_back_to_all_compounds(mapk14_counts: tuple[Path, dict], sublibrary: dict) -> None:
    path, counts_qc = mapk14_counts
    only_heldout = {**sublibrary, "random": []}
    _, qc = score_table(pl.read_parquet(path), counts_qc, only_heldout)
    assert qc["prior"]["fitted_on"] == "all"


def test_replicate_correlation_of_degenerate_tables(mapk14_counts: tuple[Path, dict], sublibrary: dict) -> None:
    path, counts_qc = mapk14_counts
    counts = pl.read_parquet(path).head(1).with_columns(target_1=pl.lit(1, pl.Int32), matrix_1=pl.lit(1, pl.Int32))
    _, qc = score_table(counts, counts_qc, sublibrary)
    assert {value for level in qc["replicate_correlation"].values() for value in level.values()} == {None}


def test_absent_scores(mapk14_scores: tuple[pl.DataFrame, dict]) -> None:
    scores, qc = mapk14_scores
    absent = absent_scores(qc)
    assert absent["count_diff"] == 0.0
    assert absent["poisson_lower"] == 0.0
    assert absent["z_faver"] < 0
    # Same as a scored compound that happens to have no counts at all.
    empty = scores.filter((pl.col("target_sum") == 0) & (pl.col("control_sum") == 0))
    assert empty.height > 0
    assert empty.item(0, "er") == pytest.approx(absent["er"])
    # The QC stores the prior rounded to 6 significant digits (byte-stable outputs), hence the tolerance.
    assert empty.item(0, "eb_log2_ratio") == pytest.approx(absent["eb_log2_ratio"], rel=1e-3)
