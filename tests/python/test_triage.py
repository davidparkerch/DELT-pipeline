import polars as pl
import pytest

from delt_pipeline.analyse.triage import FLAGS, TIERS, ScreenInput, tiers_by_score, top_by, triage, triage_qc
from delt_pipeline.library import SimulatedLibrary


def member(frame: pl.DataFrame, screen: str | None, a: int, b: int, c: int | None = None) -> pl.DataFrame:
    selected = frame.filter((pl.col("bb_a") == a) & (pl.col("bb_b") == b))
    if screen is not None:
        selected = selected.filter(pl.col("screen") == screen)
    return selected if c is None else selected.filter(pl.col("bb_c") == c)


@pytest.fixture(scope="module")
def hits(screen_inputs: dict[str, ScreenInput]) -> pl.DataFrame:
    return triage(list(screen_inputs.values()))


def test_layout(hits: pl.DataFrame) -> None:
    assert hits.get_column("screen").unique().sort().to_list() == ["ddr1", "mapk14"]
    assert hits.height == 600
    assert hits.columns[:5] == ["hit_id", "screen", "rank", "tier", "molecule_hash"]
    assert hits.get_column("hit_id").is_unique().all()
    assert hits.item(0, "hit_id") == f"ddr1:{hits.item(0, 'molecule_hash')}"
    assert set(FLAGS) <= set(hits.columns)
    assert set(hits.get_column("tier").unique()) <= set(TIERS)
    assert hits.filter(pl.col("screen") == "mapk14").get_column("rank").to_list() == list(range(1, 301))


def test_planted_binder_is_tier_one(
    hits: pl.DataFrame, screen_inputs: dict[str, ScreenInput], library: SimulatedLibrary
) -> None:
    family = member(hits, "mapk14", *library.planted.binder)
    best = family.row(0, named=True)
    assert best["tier"] == "1"
    assert best["disynthon"] == "AB"
    assert best["reasons"].startswith("3/3 replicates; disynthon AB")
    # Support comes from the other members: the family's bound without this compound is lower.
    synthons = screen_inputs["mapk14"].synthons.filter(pl.col("level") == "di_ab")
    whole = member(synthons, None, *library.planted.binder).item(0, "poisson_lower")
    assert 1.5 < best["disynthon_lower"] < whole


def test_outlier_is_not_tier_one(hits: pl.DataFrame, library: SimulatedLibrary) -> None:
    outlier = member(hits, "mapk14", *library.planted.outlier).row(0, named=True)
    assert outlier["target_support"] == 1
    assert not outlier["replicate_support"]
    assert outlier["tier"] != "1"


def test_counter_binder_is_flagged(hits: pl.DataFrame, library: SimulatedLibrary) -> None:
    """Enriched in the BCA counter-screen as well: an artefact, however good its kinase counts."""
    promiscuous = member(hits, "ddr1", *library.planted.promiscuous).row(0, named=True)
    assert promiscuous["rank"] <= 5
    assert promiscuous["counter_binder"]
    assert promiscuous["tier"] == "flagged"
    assert "counter-screen binder" in promiscuous["reasons"]


def test_shared_target_is_not_flagged(hits: pl.DataFrame, library: SimulatedLibrary) -> None:
    """A family that binds both kinases but not BCA is a plausible hinge binder, not promiscuous."""
    family = member(hits, "mapk14", *library.planted.shared)
    shared = family.filter(pl.col("shared_target"))
    assert shared.height > 0
    assert not shared.get_column("counter_binder").any()
    assert not shared.get_column("matrix_binder").any()
    assert "flagged" not in shared.get_column("tier").to_list()
    assert shared.get_column("shared_with").unique().to_list() == ["ddr1"]
    assert "also a top hit of ddr1" in shared.item(0, "reasons")


def test_bead_binder_is_flagged(screen_inputs: dict[str, ScreenInput], library: SimulatedLibrary) -> None:
    """The bead binder is enriched in target and control alike, so the lower bound never ranks it high;
    handed to triage anyway, the matrix-over-load ratio flags it."""
    mapk14 = screen_inputs["mapk14"]
    a, b, c = library.planted.bead
    scores = mapk14.scores.collect()
    bead = scores.filter((pl.col("bb_a") == a) & (pl.col("bb_b") == b) & (pl.col("bb_c") == c))
    assert bead.height == 1
    with_bead = ScreenInput("mapk14", "target", pl.concat([mapk14.top, bead]), mapk14.scores, mapk14.synthons,
                            mapk14.counts_qc)  # fmt: skip
    hits = triage([with_bead])
    row = member(hits, "mapk14", a, b, c).row(0, named=True)
    assert row["matrix_binder"]
    assert row["tier"] == "flagged"
    assert "bead binder" in row["reasons"]
    # Without counter-screens there is nothing to look up.
    assert hits.get_column("counter_lower").null_count() == hits.height
    assert not hits.get_column("shared_target").any()


def test_needs_a_target(screen_inputs: dict[str, ScreenInput]) -> None:
    with pytest.raises(ValueError, match="role 'target'"):
        triage([screen_inputs["bca"]])


def test_qc(hits: pl.DataFrame) -> None:
    qc = triage_qc(hits)
    assert qc["n_hits"] == 600
    assert set(qc["screens"]) == {"mapk14", "ddr1"}
    mapk14 = qc["screens"]["mapk14"]
    assert sum(mapk14["tiers"].values()) == mapk14["n_compounds"] == 300
    assert list(mapk14["tiers"]) == list(TIERS)
    assert mapk14["flags"]["counter_binder"] >= 1
    assert 0 < len(mapk14["top_tier1"]) <= 20
    assert all("replicates" in row["reasons"] for row in mapk14["top_tier1"])
    assert qc["thresholds"]["disynthon_lower"] == 1.5


def test_tiers_by_score(screen_inputs: dict[str, ScreenInput]) -> None:
    """The same triage applied to each score's top compounds."""
    by_score = tiers_by_score(list(screen_inputs.values()), 20)
    assert set(by_score) == {"count_diff", "er", "poisson_lower", "z_faver", "eb_log2_ratio"}
    assert all(sum(tiers.values()) == 20 for score in by_score.values() for tiers in score.values())
    # The ranking by Poisson lower bound is the one the pipeline triages.
    hits = triage(list(screen_inputs.values()))
    poisson = hits.filter(pl.col("screen") == "mapk14").head(20).get_column("tier").value_counts()
    assert dict(poisson.iter_rows()) == {tier: n for tier, n in by_score["poisson_lower"]["mapk14"].items() if n}
    top = top_by(screen_inputs["mapk14"].scores, "z_faver", 20)
    assert top.height == 20
    assert top.get_column("z_faver").is_sorted(descending=True)
    qc = triage_qc(hits, by_score=by_score)
    assert qc["by_score"] == by_score
