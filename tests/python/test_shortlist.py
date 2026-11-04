import json
from pathlib import Path

import polars as pl
import pytest
from rdkit import DataStructs

from delt_pipeline.analyse.shortlist import (
    COLUMNS,
    cluster,
    fingerprints,
    mean_similarity,
    select_picks,
    shortlist,
)

ANALOGUES = ["c1ccccc1CC(=O)N", "c1ccccc1CCC(=O)N", "CCOC(=O)C1CCNCC1"]


def test_fingerprints() -> None:
    prints = fingerprints(ANALOGUES)
    assert [print_.GetNumBits() for print_ in prints] == [2048] * 3
    assert prints[0].ToBitString() == fingerprints([ANALOGUES[0]])[0].ToBitString()
    with pytest.raises(ValueError, match="cannot parse"):
        fingerprints(["not a molecule"])


def test_cluster_by_similarity() -> None:
    prints = fingerprints(ANALOGUES)
    close = DataStructs.TanimotoSimilarity(prints[0], prints[1])
    far = DataStructs.TanimotoSimilarity(prints[0], prints[2])
    assert close > far
    threshold = (close + far) / 2
    assert cluster(prints, threshold) == [0, 0, 1]
    assert cluster(prints, close + 0.01) == [0, 1, 2]  # nothing is that alike
    # Clusters are numbered by their best-ranked member, whatever the order of the analogues.
    assert cluster([prints[2], prints[0], prints[1]], threshold) == [0, 1, 1]
    assert cluster(prints[:1], threshold) == [0]
    assert cluster([], threshold) == []
    for bad in (0.0, 1.0):
        with pytest.raises(ValueError, match="similarity"):
            cluster(prints, bad)


def test_mean_similarity() -> None:
    prints = fingerprints(ANALOGUES)
    pairs = [(0, 1), (0, 2), (1, 2)]
    expected = sum(DataStructs.TanimotoSimilarity(prints[i], prints[j]) for i, j in pairs) / 3
    assert mean_similarity(prints) == pytest.approx(expected)
    assert mean_similarity([prints[0], prints[0]]) == 1.0
    assert mean_similarity(prints[:1]) is None


def test_select_picks_by_hand() -> None:
    clusters = [0, 0, 1, 1, 2]
    keys = [None, "x", "x", "x", None]
    # Round 1: the best of every cluster (0, 2, 4). Round 2: compound 1; compound 3 would be the third of "x".
    assert select_picks(clusters, keys, n_picks=10, max_per_cluster=2, max_per_disynthon=2) == [
        (0, 1),
        (2, 1),
        (4, 1),
        (1, 2),
    ]
    assert select_picks(clusters, keys, n_picks=2) == [(0, 1), (2, 1)]
    # Without a disynthon cap the fourth compound is taken as well.
    assert select_picks(clusters, keys, n_picks=10, max_per_disynthon=9)[-1] == (3, 2)
    # A third round finds nobody left to pick.
    assert len(select_picks(clusters, keys, n_picks=10, max_per_cluster=3, max_per_disynthon=2)) == 4
    assert select_picks([], [], n_picks=3) == []
    for bad in ({"n_picks": 0}, {"max_per_cluster": 0}, {"max_per_disynthon": 0}):
        with pytest.raises(ValueError, match="at least 1"):
            select_picks(clusters, keys, **bad)


def test_shortlist_is_more_diverse_than_the_ranking(pipeline_dir: Path) -> None:
    hits = pl.read_parquet(pipeline_dir / "hits.parquet")
    picks, qc = shortlist(hits, n_picks=6, similarity=0.5)
    assert picks.columns == list(COLUMNS)
    assert picks.get_column("pick_rank").to_list() == list(range(1, picks.height + 1))
    assert set(picks.get_column("tier")) == {"1"}
    assert set(picks.get_column("screen")) == {"mapk14"}
    assert picks.get_column("hit_id").is_unique().all()
    summary = qc["screens"]["mapk14"]
    assert summary["n_picks"] == picks.height <= 6
    assert summary["n_eligible"] == hits.filter(pl.col("tier") == "1").height
    assert summary["picks"]["clusters"] >= summary["naive_top"]["clusters"]
    assert 0 < summary["picks"]["mean_similarity"] <= 1
    assert summary["picks"]["max_rank"] == picks.get_column("rank").max()
    assert qc["parameters"] == {
        "n_picks": 6,
        "similarity": 0.5,
        "max_per_cluster": 2,
        "max_per_disynthon": 3,
        "tiers": ["1"],
    }
    # Deterministic.
    again, _ = shortlist(hits, n_picks=6, similarity=0.5)
    assert again.equals(picks)


def test_shortlist_caps_picks_per_cluster_and_disynthon(pipeline_dir: Path) -> None:
    hits = pl.read_parquet(pipeline_dir / "hits.parquet")
    picks, _ = shortlist(hits, n_picks=40, similarity=0.9, max_per_cluster=1, max_per_disynthon=1)
    assert picks.get_column("cluster_id").is_unique().all()
    keys = picks.get_column("disynthon_key").drop_nulls()
    assert keys.is_unique().all()
    assert picks.get_column("cluster_size").min() >= 1
    wide, _ = shortlist(hits, n_picks=40, similarity=0.9, max_per_cluster=3, max_per_disynthon=40)
    assert set(wide.get_column("pick_round")) <= {1, 2, 3}


def test_shortlist_without_eligible_hits(pipeline_dir: Path) -> None:
    hits = pl.read_parquet(pipeline_dir / "hits.parquet")
    picks, qc = shortlist(hits, tiers=("nothing",))
    assert picks.height == 0
    assert qc["screens"]["mapk14"]["n_picks"] == 0
    assert qc["screens"]["mapk14"]["picks"] == {
        "clusters": 0,
        "disynthons": 0,
        "mean_similarity": None,
        "max_rank": None,
    }
    empty, qc = shortlist(hits.head(0))
    assert empty.height == 0
    assert qc["screens"] == {}


def test_shortlist_cli_output(pipeline_dir: Path) -> None:
    picks = pl.read_parquet(pipeline_dir / "picks.parquet")
    qc = json.loads((pipeline_dir / "picks_qc.json").read_text())
    assert picks.columns == list(COLUMNS)
    assert qc["parameters"]["n_picks"] == 6
    assert qc["screens"]["mapk14"]["n_picks"] == picks.height
