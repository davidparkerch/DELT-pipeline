from pathlib import Path

import polars as pl
import pytest

from delt_pipeline.library.blocks import (
    DESCRIPTORS,
    LOW_LOAD_RATIO,
    building_block_table,
    describe_block,
    describe_blocks,
    load_statistics,
    read_reacted,
)
from delt_pipeline.library.schema import COUNTS
from delt_pipeline.library.simulate import N_A, N_B, N_C


def test_describe_block_by_hand() -> None:
    # Glycine with the attachment point on the amine: the fragment is capped with a hydrogen.
    glycine = describe_block("N([Fr])CC(=O)O")
    assert glycine["canonical_smiles"] == "NCC(=O)O"
    assert glycine["fragment_smiles"] == "*NCC(=O)O"
    assert glycine["inchikey"] == "DHMQDGOQFOQNFH-UHFFFAOYSA-N"  # glycine
    assert glycine["attachments"] == 1
    assert glycine["heavy_atoms"] == 5
    assert glycine["mw"] == 75.067
    assert glycine["rings"] == 0
    assert glycine["hbd"] == 2
    assert glycine["fraction_csp3"] == 0.5
    assert glycine["murcko_scaffold"] is None
    assert set(DESCRIPTORS) <= set(glycine)


def test_describe_block_scaffold_stereo_and_two_attachments() -> None:
    block = describe_block("[At]NC1=CC=C(C[Fr])C=C1")
    assert block["attachments"] == 2
    assert block["rings"] == 1
    assert block["murcko_scaffold"] == "c1ccccc1"
    assert block["canonical_smiles"] == "Cc1ccc(N)cc1"
    # Stereochemistry survives capping.
    assert "@" in describe_block("CC[C@H](N[Fr])C(=O)NC")["canonical_smiles"]
    assert "@" not in describe_block("CC(N[Fr])C(=O)NC")["canonical_smiles"]
    assert describe_block("CCO")["attachments"] == 0


def test_describe_block_rejects_bad_smiles() -> None:
    with pytest.raises(ValueError, match="cannot parse"):
        describe_block("not a molecule")


def test_read_reacted(library_dir: Path, tmp_path: Path) -> None:
    reacted = read_reacted(library_dir / "synthons")
    assert reacted.columns == ["position", "bb_id", "reagent", "reacted"]
    assert reacted.height == N_A + N_B + N_C
    assert reacted.filter(pl.col("position") == "b").get_column("bb_id").to_list() == list(range(N_B))
    for position in "abc":
        (tmp_path / f"reacted_{position}.csv").write_text(f"smiles_{position}\nC\n")
    with pytest.raises(ValueError, match="no 'reacted_a' column"):
        read_reacted(tmp_path)


def test_blocks_sharing_a_fragment(library_dir: Path) -> None:
    blocks = describe_blocks(read_reacted(library_dir / "synthons"))
    shared = blocks.filter(pl.col("n_same_fragment") > 1)
    assert shared.select("position", "bb_id", "same_fragment_as").rows() == [
        ("c", N_C - 2, N_C - 2),
        ("c", N_C - 1, N_C - 2),
    ]
    assert blocks.filter(pl.col("block_id") == "A0000").height == 1
    assert blocks.filter(pl.col("n_same_fragment") == 1).height == N_A + N_B + N_C - 2


def test_load_statistics_by_hand() -> None:
    counts = pl.DataFrame(
        {"bb_a": [0, 0, 1], "bb_b": [0, 1, 1], "bb_c": [0, 0, 0], "load": [4, 0, 2]},
        schema={"bb_a": pl.UInt16, "bb_b": pl.UInt16, "bb_c": pl.UInt16, "load": pl.Int32},
    )
    sublibrary = {"heldout": [0], "random": [1]}  # two BB-As
    statistics = load_statistics(counts, sublibrary, {"a": 2, "b": 2, "c": 1})
    row = statistics.filter((pl.col("position") == "b") & (pl.col("bb_id") == 1)).row(0, named=True)
    # B1 could be in 2 (BB-A) x 1 (C) = 2 compounds; 1 of them was sequenced before selection, with load 2.
    assert (row["n_possible"], row["n_with_load"]) == (2, 1)
    assert row["mean_load"] == 1.0
    assert row["fraction_with_load"] == 0.5
    row = statistics.filter((pl.col("position") == "a") & (pl.col("bb_id") == 0)).row(0, named=True)
    assert (row["n_possible"], row["mean_load"]) == (2, 2.0)  # 2 B x 1 C possible; load 4 + 0


def test_building_block_table(pipeline_dir: Path, sublibrary: dict, tmp_path: Path) -> None:
    table = pl.read_parquet(pipeline_dir / "blocks.parquet")
    assert table.height == N_A + N_B + N_C
    assert table.get_column("block_id").is_unique().all()
    assert table.equals(table.sort("position", "bb_id"))
    sources = table.filter(pl.col("position") == "a").get_column("sublibrary_source").value_counts()
    assert dict(sources.iter_rows()) == {
        "heldout": len(sublibrary["heldout"]),
        "random": len(sublibrary["random"]),
        "excluded": N_A - 8,
    }
    assert set(table.filter(pl.col("position") != "a").get_column("sublibrary_source")) == {"all"}
    excluded = table.filter(pl.col("sublibrary_source") == "excluded")
    assert excluded.get_column("mean_load").null_count() == excluded.height
    # The load of a block is its share of the sub-library's total load.
    counts = pl.read_parquet(pipeline_dir / "mapk14.counts.parquet")
    block = table.filter((pl.col("position") == "b") & (pl.col("bb_id") == 5)).row(0, named=True)
    expected = counts.filter(pl.col("bb_b") == 5).get_column("load").sum() / block["n_possible"]
    assert block["mean_load"] == pytest.approx(expected, rel=1e-4)
    assert set(table.get_column("low_load").drop_nulls()) <= {True, False}


def test_building_block_qc(pipeline_dir: Path, library_dir: Path, sublibrary: dict, tmp_path: Path) -> None:
    import json  # noqa: PLC0415

    qc = json.loads((pipeline_dir / "blocks_qc.json").read_text())
    assert qc["n_blocks"] == {"a": N_A, "b": N_B, "c": N_C}
    assert qc["n_in_sublibrary"] == 8 + N_B + N_C
    assert qc["n_shared_fragments"] == 2
    assert qc["load_shared_across_screens"] is True
    assert qc["low_load_ratio"] == LOW_LOAD_RATIO
    # A screen with a different pre-selection load is reported, not hidden.
    other = pl.read_parquet(pipeline_dir / "bca.counts.parquet").with_columns(load=pl.col("load") + 1)
    other.write_parquet(tmp_path / "other.counts.parquet")
    _table, qc = building_block_table(
        read_reacted(library_dir / "synthons"),
        [pipeline_dir / "mapk14.counts.parquet", tmp_path / "other.counts.parquet"],
        sublibrary,
    )
    assert qc["load_shared_across_screens"] is False
    with pytest.raises(ValueError, match="at least one count table"):
        building_block_table(read_reacted(library_dir / "synthons"), [], sublibrary)
    assert set(COUNTS) >= {"load"}
