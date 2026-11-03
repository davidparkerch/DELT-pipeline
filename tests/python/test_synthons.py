import polars as pl

from delt_pipeline.analyse import Exposure
from delt_pipeline.analyse.score import exposure_from_qc
from delt_pipeline.analyse.synthons import LEVELS, aggregate_synthons, synthon_qc
from delt_pipeline.library import SimulatedLibrary


def tiny() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "bb_a": [0, 0, 1],
            "bb_b": [0, 1, 1],
            "bb_c": [0, 0, 0],
            "target_1": [1, 2, 0],
            "target_2": [0, 1, 0],
            "target_3": [0, 0, 0],
            "matrix_1": [1, 0, 2],
            "matrix_2": [0, 0, 0],
            "matrix_3": [0, 0, 0],
            "load": [5, 6, 7],
        },
        schema_overrides=dict.fromkeys(("bb_a", "bb_b", "bb_c"), pl.UInt16),
    )


def test_hand_computed_sums() -> None:
    synthons = aggregate_synthons(tiny().lazy(), Exposure(10.0, 10.0))
    assert synthons.get_column("level").unique(maintain_order=True).to_list() == list(LEVELS)
    a0 = synthons.filter((pl.col("level") == "mono_a") & (pl.col("bb_a") == 0)).row(0, named=True)
    assert (a0["n_compounds"], a0["n_with_target"], a0["target_sum"], a0["control_sum"], a0["load"]) == (2, 2, 4, 1, 11)
    assert a0["bb_b"] is None
    assert a0["bb_c"] is None
    assert a0["poisson_ratio"] == 4.0
    # BB-A 1 has no target counts: lower bound 0, ratio 0.
    a1 = synthons.filter((pl.col("level") == "mono_a") & (pl.col("bb_a") == 1)).row(0, named=True)
    assert (a1["poisson_lower"], a1["poisson_ratio"]) == (0.0, 0.0)
    # Disynthon B1-C0 pools compounds from two BB-As; a group without control counts has no upper bound.
    bc = synthons.filter((pl.col("level") == "di_bc") & (pl.col("bb_b") == 1)).row(0, named=True)
    assert (bc["n_compounds"], bc["target_sum"], bc["control_sum"]) == (2, 3, 2)
    ab = synthons.filter((pl.col("level") == "di_ab") & (pl.col("bb_b") == 1) & (pl.col("bb_a") == 0)).row(
        0, named=True
    )
    assert ab["poisson_upper"] is None
    assert ab["poisson_ratio"] is None
    # Within a level, rows are sorted by the lower bound.
    lower = synthons.filter(pl.col("level") == "di_ab").get_column("poisson_lower").to_list()
    assert lower == sorted(lower, reverse=True)


def test_planted_disynthon(
    mapk14_scores: tuple[pl.DataFrame, dict], mapk14_counts: tuple, library: SimulatedLibrary
) -> None:
    """The planted binder pair tops its level; the single-count outlier's pair does not come close."""
    scores, _ = mapk14_scores
    synthons = aggregate_synthons(scores.lazy(), exposure_from_qc(mapk14_counts[1]))
    planted = library.planted
    pairs = synthons.filter(pl.col("level") == "di_ab")
    di_ab = pairs.select("bb_a", "bb_b").rows()
    assert set(di_ab[:2]) == {planted.binder, planted.shared}
    # Six lucky counts barely move a pair of ~15 compounds: its bound stays near 1, half the binder's.
    lower = dict(zip(di_ab, pairs.get_column("poisson_lower"), strict=True))
    assert lower[planted.outlier[:2]] < 0.5 * lower[planted.binder]
    qc = synthon_qc(synthons)
    assert set(qc["levels"]) == set(LEVELS)
    assert qc["levels"]["di_ab"]["n_synthons"] == len(di_ab)
    assert qc["levels"]["di_ab"]["replicate_correlation"] > 0.5
    assert (qc["top"]["di_ab"][0]["bb_a"], qc["top"]["di_ab"][0]["bb_b"]) in {planted.binder, planted.shared}
    assert set(qc["top"]) == {"di_ab", "di_ac", "di_bc"}


def test_qc_without_variation() -> None:
    synthons = aggregate_synthons(tiny().head(1).lazy(), Exposure(10.0, 10.0))
    assert synthon_qc(synthons)["levels"]["mono_a"]["replicate_correlation"] is None
