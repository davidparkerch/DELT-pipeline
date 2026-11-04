import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from delt_pipeline.analyse.model import (
    CONTROL_PSEUDOCOUNT,
    VARIANTS,
    fingerprint_matrix,
    heldout_blocks,
    model_report,
    published_ml,
    sharing_a_disynthon,
    training_sample,
)
from delt_pipeline.library import SimulatedLibrary, read_building_blocks

SMALL = {"n_top": 150, "n_random": 150}


def test_published_ml_table() -> None:
    for target in ("mapk14", "ddr1"):
        rows = published_ml(target)
        assert len(rows) == 16  # 8 models x 2 splits
        assert {row["split"] for row in rows} == {"random", "disynthon"}
        assert {row["model"] for row in rows} >= {"RF", "XGBoost", "Chemprop (M)", "DEL-Compose (S)"}
        assert all(0 < row["on_dna_sd"] < 0.2 for row in rows)
    # Two values checked by eye against Table 1 of the paper.
    mapk14 = {(row["model"], row["split"]): row for row in published_ml("MAPK14")}
    assert mapk14[("RF", "random")]["on_dna"] == 0.694
    assert mapk14[("DEL-Compose (S)", "disynthon")]["off_dna"] == 0.716
    assert published_ml("bca") == []


def test_training_sample(mapk14_scores: tuple[pl.DataFrame, dict]) -> None:
    scores, _qc = mapk14_scores
    sample = training_sample(scores.lazy(), **SMALL)
    assert sample.columns == ["molecule_hash", "smiles", "bb_a", "bb_b", "bb_c", "target_sum", "control_sum"]
    assert sample.get_column("molecule_hash").is_unique().all()
    assert sample.equals(training_sample(scores.lazy(), **SMALL))  # deterministic
    assert SMALL["n_top"] < sample.height < SMALL["n_top"] + 2 * SMALL["n_random"]
    # The top part holds every compound with more counts than its weakest member.
    top_floor = scores.sort("target_sum", descending=True).item(SMALL["n_top"] - 1, "target_sum")
    clearly_above = scores.filter(pl.col("target_sum") > top_floor).get_column("molecule_hash")
    assert clearly_above.is_in(sample.get_column("molecule_hash").implode()).all()
    assert training_sample(scores.lazy(), 5, 0).height >= 5


def test_fingerprint_matrix() -> None:
    matrix = fingerprint_matrix(["c1ccccc1CC(=O)N", "CCO"])
    assert matrix.shape == (2, 2048)
    assert matrix.dtype == np.float32
    assert set(np.unique(matrix)) == {0.0, 1.0}
    assert matrix[0].sum() > matrix[1].sum() > 0
    with pytest.raises(ValueError, match="cannot parse"):
        fingerprint_matrix(["not a molecule"])


def test_heldout_blocks(library: SimulatedLibrary, library_dir: Path, mapk14_scores: tuple[pl.DataFrame, dict]) -> None:
    scores, _qc = mapk14_scores
    blocks = read_building_blocks(library_dir / "synthons")
    held = library.heldout["mapk14_ondna"]
    result = heldout_blocks(held, blocks)
    in_library = held.filter(pl.col("molecule_hash").is_not_null())
    assert result.height == in_library.height
    # The ids agree with the ones in the scored table.
    joined = result.join(scores.select("molecule_hash", "bb_a", "bb_b", "bb_c"), on="molecule_hash", suffix="_scored")
    assert joined.height > 0
    for position in "abc":
        assert (joined.get_column(f"bb_{position}") == joined.get_column(f"bb_{position}_scored")).all()
    unknown = pl.DataFrame({"molecule_hash": ["x"], "smiles_a": ["nope"], "smiles_b": ["nope"], "smiles_c": ["nope"]})
    assert heldout_blocks(unknown, blocks).height == 0


def test_sharing_a_disynthon() -> None:
    schema = {"bb_a": pl.UInt16, "bb_b": pl.UInt16, "bb_c": pl.UInt16}
    sample = pl.DataFrame(
        {
            "molecule_hash": ["ab", "ab2", "ac", "bc", "none"],
            "bb_a": [0, 0, 0, 7, 9],
            "bb_b": [0, 0, 5, 1, 9],
            "bb_c": [1, 2, 3, 3, 9],
        },
        schema={"molecule_hash": pl.String, **schema},
    )
    held = pl.DataFrame({"bb_a": [0], "bb_b": [1], "bb_c": [3]}, schema=schema)
    # Held-out (0, 1, 3): AB = (0, 1), AC = (0, 3), BC = (1, 3). "ab"/"ab2" share none (their B is 0)...
    assert sorted(sharing_a_disynthon(sample, held).to_list()) == ["ac", "bc"]
    other = pl.DataFrame({"bb_a": [0], "bb_b": [0], "bb_c": [7]}, schema=schema)
    assert sorted(sharing_a_disynthon(sample, other).to_list()) == ["ab", "ab2"]


def test_model_report(mapk14_scores: tuple[pl.DataFrame, dict], library: SimulatedLibrary, library_dir: Path) -> None:
    scores, _qc = mapk14_scores
    blocks = read_building_blocks(library_dir / "synthons")
    heldout = {"on": library.heldout["mapk14_ondna"], "off": library.heldout["mapk14_offdna"]}
    report = model_report(scores.lazy(), blocks, heldout, "mapk14", **SMALL)
    assert [(r["variant"], r["split"], r["condition"]) for r in report["results"]] == [
        (variant, split, condition)
        for split in ("random", "disynthon")
        for variant in VARIANTS
        for condition in ("on", "off")
    ]
    assert all(r["spearman"] is None or -1 <= r["spearman"] <= 1 for r in report["results"])
    in_library = {c: t.filter(pl.col("molecule_hash").is_not_null()).height for c, t in heldout.items()}
    assert {(r["condition"], r["n_compounds"]) for r in report["results"]} == set(in_library.items())
    training = report["training"]
    assert training["disynthon"]["n_train"] < training["random"]["n_train"]
    assert training["n_excluded_disynthon"] == training["random"]["n_train"] - training["disynthon"]["n_train"] > 0
    assert training["sample"]["n_top"] == SMALL["n_top"]
    assert report["published"]["on"]["poisson"] == 0.737
    assert len(report["published_ml"]) == 16
    assert report["model"]["loss"] == "poisson"
    # The planted binder families have low KD, so a model that learns them ranks the validated compounds well.
    random_on = {
        r["variant"]: r["spearman"] for r in report["results"] if r["split"] == "random" and r["condition"] == "on"
    }
    assert random_on["control_weighted"] > 0.2
    assert model_report(scores.lazy(), blocks, heldout, "mapk14", **SMALL) == report  # deterministic
    assert CONTROL_PSEUDOCOUNT > 0


def test_model_without_heldout_data(mapk14_scores: tuple[pl.DataFrame, dict], library_dir: Path) -> None:
    scores, _qc = mapk14_scores
    report = model_report(scores.lazy(), read_building_blocks(library_dir / "synthons"), {}, "bca", **SMALL)
    assert report["results"] == []
    assert report["training"] == {}
    assert report["published_ml"] == []


def test_model_cli_output(pipeline_dir: Path) -> None:
    mapk14 = json.loads((pipeline_dir / "mapk14.model.json").read_text())
    assert len(mapk14["results"]) == 8
    assert json.loads((pipeline_dir / "bca.model.json").read_text())["results"] == []
