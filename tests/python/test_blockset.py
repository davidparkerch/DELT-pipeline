import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from scipy import stats

from delt_pipeline.analyse.blockset import Q_LEVEL, blockset_enrichment, blockset_qc


def _scores(seed: int = 3, n: int = 400) -> pl.DataFrame:
    rng = np.random.default_rng(seed)
    frame = pl.DataFrame(
        {
            "bb_a": rng.integers(0, 4, n),
            "bb_b": rng.integers(0, 5, n),
            "bb_c": rng.integers(0, 3, n),
            "poisson_lower": np.round(rng.gamma(1.0, 1.0, n), 1),  # rounded: plenty of ties
        }
    )
    # Block B3 is shifted upwards.
    return frame.with_columns(
        poisson_lower=pl.when(pl.col("bb_b") == 3)
        .then(pl.col("poisson_lower") + 1.5)
        .otherwise(pl.col("poisson_lower"))
    )


def test_matches_the_mann_whitney_test() -> None:
    scores = _scores()
    table = blockset_enrichment(scores.lazy())
    assert table.height == 4 + 5 + 3
    value = scores.get_column("poisson_lower").to_numpy()
    for row in table.iter_rows(named=True):
        member = (scores.get_column(f"bb_{row['position']}") == row["bb_id"]).to_numpy()
        reference = stats.mannwhitneyu(value[member], value[~member], alternative="greater", method="asymptotic")
        assert row["p_value"] == pytest.approx(reference.pvalue, rel=1e-4, abs=1e-12)
        assert row["auc"] == pytest.approx(reference.statistic / (member.sum() * (~member).sum()), rel=1e-5)
        assert row["n_compounds"] == member.sum()
        assert row["median_score"] == pytest.approx(np.median(value[member]), rel=1e-5)
    first = table.row(0, named=True)
    assert (first["position"], first["bb_id"]) == ("b", 3)
    assert first["q_value"] < Q_LEVEL
    assert first["auc"] > 0.6
    q = stats.false_discovery_control(table.get_column("p_value").to_numpy(), method="bh")
    assert table.get_column("q_value").to_numpy() == pytest.approx(q, rel=1e-4, abs=1e-12)


def test_counts_blocks_in_the_top_table() -> None:
    scores = _scores()
    top = scores.sort("poisson_lower", descending=True).head(50)
    table = blockset_enrichment(scores.lazy(), top)
    row = table.filter((pl.col("position") == "b") & (pl.col("bb_id") == 3)).row(0, named=True)
    assert row["n_top"] == top.filter(pl.col("bb_b") == 3).height > 0
    assert table.get_column("n_top").sum() == 3 * 50
    assert blockset_enrichment(scores.lazy()).get_column("n_top").sum() == 0


def test_no_signal_gives_no_enrichment() -> None:
    flat = _scores().with_columns(poisson_lower=pl.lit(1.0))
    table = blockset_enrichment(flat.lazy())
    assert (table.get_column("auc") == 0.5).all()
    assert (table.get_column("z") == 0).all()
    assert (table.get_column("p_value") == 0.5).all()
    assert blockset_enrichment(flat.head(1).lazy()).height == 0  # one compound: no block to compare
    nulls = _scores().with_columns(poisson_lower=pl.lit(None, dtype=pl.Float64))
    assert (blockset_enrichment(nulls.lazy()).get_column("auc") == 0.5).all()  # null counts as 0


def test_planted_families_lead(pipeline_dir: Path) -> None:
    table = pl.read_parquet(pipeline_dir / "mapk14.blockset.parquet")
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    first = table.row(0, named=True)
    # MAPK14 has two planted families: (A2, B5) and the kinase-shared (A7, B11). Their B blocks lead the ranking
    # and the best one is significant. The A blocks hold 16 family compounds among hundreds of others, so a
    # block-level test cannot see them: a set test finds blocks that dominate, the disynthon table finds lines.
    assert {tuple(row) for row in table.head(2).select("position", "bb_id").rows()} == {("b", 5), ("b", 11)}
    assert first["q_value"] < Q_LEVEL
    assert first["n_top"] == top.filter(pl.col(f"bb_{first['position']}") == first["bb_id"]).height
    # Control: the BCA family block (B20) is not called in MAPK14.
    bca_block = table.filter((pl.col("position") == "b") & (pl.col("bb_id") == 20)).row(0, named=True)
    assert bca_block["q_value"] > 0.5


def test_blockset_qc(pipeline_dir: Path) -> None:
    table = pl.read_parquet(pipeline_dir / "mapk14.blockset.parquet")
    qc = json.loads((pipeline_dir / "mapk14.blockset_qc.json").read_text())
    assert qc == blockset_qc(table, qc["n_compounds"])
    assert qc["score"] == "poisson_lower"
    assert qc["n_tests"] == table.height
    assert sum(qc["n_significant"].values()) == table.filter(pl.col("q_value") < Q_LEVEL).height
    assert len(qc["top"]) == 10
    assert qc["top"][0]["position"] == table.row(0, named=True)["position"]
