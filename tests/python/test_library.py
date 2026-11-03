from pathlib import Path

import numpy as np
import polars as pl
import pytest

from delt_pipeline.library import (
    SimulatedLibrary,
    heldout_bb_a,
    library_depth,
    library_diversity,
    read_building_blocks,
    select_sublibrary,
    simulate_library,
    write_counts,
)
from delt_pipeline.library.schema import BB_COLUMNS, COUNTS, KINDEL_COUNTS
from delt_pipeline.library.simulate import N_A, N_B, N_C, molecule_hash


def test_simulation_is_deterministic(library: SimulatedLibrary) -> None:
    again = simulate_library()
    for screen, frame in library.screens.items():
        assert frame.equals(again.screens[screen])
    assert not library.screens["mapk14"].equals(simulate_library(seed=8).screens["mapk14"])


def test_simulation_mirrors_kindel(library: SimulatedLibrary) -> None:
    assert set(library.screens) == {"mapk14", "ddr1", "bca"}
    frame = library.screens["mapk14"]
    assert set(KINDEL_COUNTS) <= set(frame.columns)
    assert frame.height < N_A * N_B * N_C  # all-zero compounds are absent
    assert frame.select(pl.sum_horizontal(KINDEL_COUNTS).min()).item() > 0
    first = frame.row(0, named=True)
    assert first["molecule_hash"] == molecule_hash(first["smiles"])
    heldout = library.heldout["mapk14_ondna"]
    assert heldout.get_column("molecule_hash").null_count() == 2  # external compounds
    assert (heldout.get_column("kd") > 0).all()


def test_read_building_blocks(library_dir: Path, tmp_path: Path) -> None:
    building_blocks = read_building_blocks(library_dir / "synthons")
    assert [building_blocks[p].height for p in "abc"] == [N_A, N_B, N_C]
    assert building_blocks["a"].get_column("bb_a").to_list() == list(range(N_A))
    assert library_diversity(building_blocks) == N_A * N_B * N_C

    for position in "bc":
        (tmp_path / f"reacted_{position}.csv").write_text(f"smiles_{position}\nC\nCC\n")
    (tmp_path / "reacted_a.csv").write_text("smiles\nC\n")
    with pytest.raises(ValueError, match="no 'smiles_a' column"):
        read_building_blocks(tmp_path)
    (tmp_path / "reacted_a.csv").write_text("smiles_a\nC\nC\n")
    with pytest.raises(ValueError, match="more than once"):
        read_building_blocks(tmp_path)


def test_heldout_bb_a(library: SimulatedLibrary, library_dir: Path) -> None:
    building_blocks = read_building_blocks(library_dir / "synthons")
    ids, unmatched = heldout_bb_a(building_blocks, list(library.heldout.values()))
    planted = library.planted
    assert {planted.binder[0], planted.ddr1_binder[0], planted.shared[0]} <= set(ids)
    assert unmatched == 0
    unknown = pl.DataFrame({"smiles_a": ["not-a-building-block"]})
    assert heldout_bb_a(building_blocks, [unknown]) == ([], 1)
    assert heldout_bb_a(building_blocks, []) == ([], 0)


def test_select_sublibrary(library: SimulatedLibrary, library_dir: Path, sublibrary: dict) -> None:
    building_blocks = read_building_blocks(library_dir / "synthons")
    assert len(sublibrary["heldout"]) + len(sublibrary["random"]) == 8
    assert not set(sublibrary["heldout"]) & set(sublibrary["random"])
    assert sublibrary["diversity"] == N_A * N_B * N_C
    # Fewer BB-As than held-out ones: keep the held-out ones, add none.
    small = select_sublibrary(building_blocks, list(library.heldout.values()), n_bb_a=1, seed=1)
    assert small["heldout"] == sublibrary["heldout"]
    assert small["random"] == []
    # More than the library has: everything.
    everything = select_sublibrary(building_blocks, [], n_bb_a=1000, seed=1)
    assert len(everything["random"]) == N_A
    with pytest.raises(ValueError, match="at least 1"):
        select_sublibrary(building_blocks, [], n_bb_a=0, seed=1)


def test_depth_comes_from_the_whole_library(library_dir: Path, mapk14_counts: tuple[Path, dict]) -> None:
    _path, qc = mapk14_counts
    full = pl.read_parquet(library_dir / "mapk14.parquet")
    assert qc["library"]["n_compounds"] == full.height
    assert qc["library"]["depth"]["target_1"] == full.get_column("seq_target_1").sum()
    # The sub-library is smaller, and its own totals must not be used for normalisation.
    assert qc["sublibrary"]["depth"]["target_1"] < qc["library"]["depth"]["target_1"]
    assert library_depth(library_dir / "mapk14.parquet") == qc["library"]


def test_count_table(mapk14_counts: tuple[Path, dict], sublibrary: dict) -> None:
    path, qc = mapk14_counts
    counts = pl.read_parquet(path)
    assert counts.columns == ["molecule_hash", "smiles", *BB_COLUMNS, *COUNTS]
    assert set(counts.get_column("bb_a").unique()) <= set(sublibrary["heldout"] + sublibrary["random"])
    assert counts.equals(counts.sort(*BB_COLUMNS, "molecule_hash"))
    assert qc["sublibrary"]["n_compounds"] == counts.height
    assert qc["sublibrary"]["n_unmapped"] == 0
    assert 0 < qc["sublibrary"]["fraction_zero"]["target_1"] < 1


def test_unmapped_building_blocks_are_counted(library_dir: Path, tmp_path: Path, sublibrary: dict) -> None:
    building_blocks = read_building_blocks(library_dir / "synthons")
    building_blocks["c"] = building_blocks["c"].head(N_C - 1)  # forget the last C
    qc = write_counts(library_dir / "mapk14.parquet", building_blocks, sublibrary, tmp_path / "c.parquet")
    counts = pl.read_parquet(tmp_path / "c.parquet")
    assert qc["sublibrary"]["n_unmapped"] > 0
    assert counts.get_column("bb_c").max() == N_C - 2
    assert np.isclose(qc["sublibrary"]["n_compounds"], counts.height)
