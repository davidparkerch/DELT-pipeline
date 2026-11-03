import polars as pl

from delt_pipeline.analyse import PUBLISHED, SCORES, absent_scores, read_heldout, validate_scores
from delt_pipeline.library import SimulatedLibrary


def test_read_heldout_keeps_in_library_compounds(library: SimulatedLibrary) -> None:
    raw = library.heldout["mapk14_ondna"]
    heldout = read_heldout(raw)
    assert heldout.columns == ["molecule_hash", "kd"]
    assert heldout.height == raw.height - 2  # the two external compounds have no hash
    repeated = read_heldout(pl.concat([raw, raw.with_columns(kd=pl.col("kd") * 3)]))
    assert repeated.height == 2 * heldout.height  # every measurement counts, as in the KinDEL benchmark


def test_validate_scores(mapk14_scores: tuple[pl.DataFrame, dict], library: SimulatedLibrary) -> None:
    scores, qc = mapk14_scores
    heldout = read_heldout(library.heldout["mapk14_ondna"])
    result = validate_scores(scores, heldout, absent_scores(qc), "MAPK14", "on")
    assert result["n_compounds"] == heldout.height
    # Held-out compounds with no counts anywhere are absent from the table (as in KinDEL), not dropped.
    present = heldout.join(scores.select("molecule_hash"), on="molecule_hash").height
    assert result["n_absent"] == heldout.height - present
    assert [row["score"] for row in result["scores"]] == list(SCORES)
    by_score = {row["score"]: row for row in result["scores"]}
    # The planted families have low KD, so evidence-aware scores rank them well.
    assert by_score["poisson_lower"]["spearman"] > 0.4
    assert by_score["poisson_lower"]["published"] == PUBLISHED[("mapk14", "on")]["poisson"]
    assert by_score["count_diff"]["published"] == PUBLISHED[("mapk14", "on")]["count_diff"]
    assert by_score["er"]["published"] is None


def test_absent_compounds_get_zero_count_scores(
    mapk14_scores: tuple[pl.DataFrame, dict], library: SimulatedLibrary
) -> None:
    scores, qc = mapk14_scores
    heldout = read_heldout(library.heldout["mapk14_offdna"])
    baseline = validate_scores(scores, heldout, absent_scores(qc), "x", "off")["n_absent"]
    missing = heldout.join(scores.select("molecule_hash"), on="molecule_hash").item(0, "molecule_hash")
    result = validate_scores(scores.filter(pl.col("molecule_hash") != missing), heldout, absent_scores(qc), "x", "off")
    assert result["n_absent"] == baseline + 1
    assert result["n_compounds"] == heldout.height
    assert all(row["published"] is None for row in result["scores"])


def test_undefined_correlations(mapk14_scores: tuple[pl.DataFrame, dict], library: SimulatedLibrary) -> None:
    scores, qc = mapk14_scores
    heldout = read_heldout(library.heldout["mapk14_ondna"])
    constant = scores.with_columns(pl.lit(1.0).alias(name) for name in SCORES)
    result = validate_scores(constant, heldout, dict.fromkeys(SCORES, 1.0), "mapk14", "on")
    assert {row["spearman"] for row in result["scores"]} == {None}
    tiny = validate_scores(scores, heldout.head(2), absent_scores(qc), "mapk14", "on")
    assert {row["spearman"] for row in tiny["scores"]} == {None}
