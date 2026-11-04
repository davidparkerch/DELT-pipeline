"""Shared fixtures: small deterministic feature tables and the simulated DEL library."""

import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from delt_pipeline.analyse import score_table
from delt_pipeline.io import FeatureTable
from delt_pipeline.library import (
    SimulatedLibrary,
    read_building_blocks,
    select_sublibrary,
    simulate_library,
    write_counts,
)

N_ROWS = 24


def _features() -> pl.DataFrame:
    rng = np.random.default_rng(7)
    label = np.array(["active", "inactive"] * (N_ROWS // 2))
    signal = np.where(label == "active", 2.0, -2.0)
    return pl.DataFrame(
        {
            "Sample_ID": [f"s{index:02d}" for index in range(N_ROWS)],
            "smiles": ["CCO"] * N_ROWS,
            "feat_alpha": signal + rng.normal(0, 0.5, N_ROWS),
            "feat_beta": rng.normal(10, 2, N_ROWS),
            "gamma": -signal + rng.normal(0, 1, N_ROWS),
            "label": label,
        }
    )


@pytest.fixture
def table() -> FeatureTable:
    """Clean binary-classification table."""
    return FeatureTable(_features(), "Sample_ID", "label")


@pytest.fixture
def unlabelled() -> FeatureTable:
    """Clean table without a target."""
    return FeatureTable(_features().drop("label"), "Sample_ID")


@pytest.fixture
def regression() -> FeatureTable:
    """Table with a continuous target."""
    rng = np.random.default_rng(11)
    x = rng.normal(0, 1, (N_ROWS, 3))
    frame = pl.DataFrame(
        {
            "id": [f"r{index:02d}" for index in range(N_ROWS)],
            "x1": x[:, 0],
            "x2": x[:, 1],
            "x3": x[:, 2],
            "y": 3 * x[:, 0] + rng.normal(0, 0.1, N_ROWS),
        }
    )
    return FeatureTable(frame, "id", "y")


@pytest.fixture(scope="session")
def library() -> SimulatedLibrary:
    """The simulated library (same seed as assets/example)."""
    return simulate_library()


@pytest.fixture(scope="session")
def library_dir(tmp_path_factory: pytest.TempPathFactory, library: SimulatedLibrary) -> Path:
    """The simulated library written in the KinDEL file layout."""
    root = tmp_path_factory.mktemp("library")
    (root / "heldout").mkdir()
    (root / "synthons").mkdir()
    for screen, frame in library.screens.items():
        frame.write_parquet(root / f"{screen}.parquet")
    for name, frame in library.heldout.items():
        frame.write_csv(root / "heldout" / f"{name}.csv")
    for position, frame in library.building_blocks.items():
        frame.write_csv(root / "synthons" / f"reacted_{position}.csv")
    return root


@pytest.fixture(scope="session")
def sublibrary(library: SimulatedLibrary, library_dir: Path) -> dict:
    """Eight BB-As: every held-out one plus random ones."""
    building_blocks = read_building_blocks(library_dir / "synthons")
    return select_sublibrary(building_blocks, list(library.heldout.values()), n_bb_a=8, seed=42)


@pytest.fixture(scope="session")
def mapk14_counts(tmp_path_factory: pytest.TempPathFactory, library_dir: Path, sublibrary: dict) -> tuple[Path, dict]:
    """The MAPK14 count table of the sub-library and its QC."""
    output = tmp_path_factory.mktemp("counts") / "mapk14.counts.parquet"
    qc = write_counts(
        library_dir / "mapk14.parquet", read_building_blocks(library_dir / "synthons"), sublibrary, output
    )
    return output, json.loads(json.dumps(qc))


@pytest.fixture(scope="session")
def mapk14_scores(mapk14_counts: tuple[Path, dict], sublibrary: dict) -> tuple[pl.DataFrame, dict]:
    """Scored MAPK14 sub-library and its QC."""
    path, qc = mapk14_counts
    return score_table(pl.read_parquet(path), qc, sublibrary)


def run_cli(*args: str | Path) -> int:
    """Run the CLI in-process with stringified arguments."""
    from delt_pipeline.cli import main  # noqa: PLC0415 - keep the CLI import local to the helper

    return main([str(arg) for arg in args])


@pytest.fixture(scope="session")
def pipeline_dir(tmp_path_factory: pytest.TempPathFactory, library_dir: Path) -> Path:
    """Every CLI step for two screens (MAPK14 with held-out KD, BCA without), like the pipeline runs them."""
    out = tmp_path_factory.mktemp("pipeline")
    synthons = library_dir / "synthons"
    heldout = sorted((library_dir / "heldout").glob("*.csv"))
    sub = out / "sublibrary.json"
    assert (
        run_cli("select-sublibrary", "--synthons", synthons, "--heldout", *heldout, "--n-bb-a", 8, "--output", sub) == 0
    )
    for screen, target in (("mapk14", "MAPK14"), ("bca", "BCA")):
        counts, qc = out / f"{screen}.counts.parquet", out / f"{screen}.counts_qc.json"
        assert run_cli("subset", "--input", library_dir / f"{screen}.parquet", "--synthons", synthons,
                       "--sublibrary", sub, "--output", counts, "--qc", qc) == 0  # fmt: skip
        scores, top, score_qc = (
            out / f"{screen}.scores.parquet",
            out / f"{screen}.top.parquet",
            out / f"{screen}.score_qc.json",
        )
        assert run_cli("score", "--input", counts, "--counts-qc", qc, "--sublibrary", sub, "--top-n", 300,
                       "--output", scores, "--top", top, "--qc", score_qc) == 0  # fmt: skip
        pooled, synthon_qc = out / f"{screen}.synthons.parquet", out / f"{screen}.synthon_qc.json"
        assert run_cli("synthons", "--scores", scores, "--counts-qc", qc, "--output", pooled,
                       "--qc", synthon_qc) == 0  # fmt: skip
        assert run_cli("blockset", "--scores", scores, "--top", top, "--output", out / f"{screen}.blockset.parquet",
                       "--qc", out / f"{screen}.blockset_qc.json") == 0  # fmt: skip
        validation = out / f"{screen}.validation.json"
        assays = []
        if screen == "mapk14":
            assays = [
                "--ondna",
                library_dir / "heldout/mapk14_ondna.csv",
                "--offdna",
                library_dir / "heldout/mapk14_offdna.csv",
            ]
        assert run_cli("validate", "--scores", scores, "--score-qc", score_qc, "--target", target.lower(), *assays,
                       "--output", validation) == 0  # fmt: skip
        exports = out / "exports" / screen
        assert run_cli("export", "--input", top, "--id-column", "molecule_hash", "--formats", "parquet,tsv",
                       "--prefix", f"{screen}_hits", "--outdir", exports) == 0  # fmt: skip
        assert run_cli("describe", "--input", top, "--id-column", "molecule_hash", "--id", screen,
                       "--exports", *sorted(exports.iterdir()), "--exports-base", f"../exports/{screen}/",
                       "--dictionary", Path(__file__).parents[2] / "assets/hits_dictionary.csv",
                       "--license", "https://github.com/insitro/kindel/blob/main/LICENSE.md",
                       "--source", "10.48550/arXiv.2410.08938", "--creator-name", "Ada",
                       "--creator-orcid", "0000-0002-1825-0097",
                       "--description", f"Simulated {target} screen", "--keywords", "DEL; kinase",
                       "--pipeline", "demo", "--pipeline-version", "1.0", "--homepage", "https://example.org/demo",
                       "--jsonld", out / f"{screen}.metadata.jsonld", "--datapackage", exports / "datapackage.json",
                       "--identifier", out / f"{screen}.identifier.txt") == 0  # fmt: skip
    assert run_cli("blocks", "--synthons", synthons, "--sublibrary", sub,
                   "--counts", out / "mapk14.counts.parquet", out / "bca.counts.parquet",
                   "--output", out / "blocks.parquet",
                   "--qc", out / "blocks_qc.json") == 0  # fmt: skip
    roles = ["--screen", "mapk14", "target", "--screen", "bca", "counter"]
    assert run_cli("triage", *roles, "--input-dir", out, "--compare-top-n", 50,
                   "--output", out / "hits.parquet", "--qc", out / "triage_qc.json") == 0  # fmt: skip
    assert run_cli("shortlist", "--hits", out / "hits.parquet", "--n-picks", 6, "--similarity", 0.5,
                   "--output", out / "picks.parquet", "--qc", out / "picks_qc.json") == 0  # fmt: skip
    assert run_cli("explore", "--hits", out / "hits.parquet", "--picks", out / "picks.parquet",
                   "--blocks", out / "blocks.parquet", "--blockset", "mapk14", out / "mapk14.blockset.parquet",
                   "--blockset", "bca", out / "bca.blockset.parquet", "--metadata", out / "mapk14.metadata.jsonld",
                   "--title", "Test explorer", "--report-url", "index.html", "--pipeline", "demo",
                   "--pipeline-version", "1.0", "--homepage", "https://example.org/demo",
                   "--output", out / "explorer.html") == 0  # fmt: skip
    for screen, target in (("mapk14", "MAPK14"), ("bca", "BCA")):
        assert run_cli("visualise", "--top", out / f"{screen}.top.parquet", "--validation",
                       out / f"{screen}.validation.json", "--screen", target,
                       "--scores", out / f"{screen}.scores.parquet", "--synthons", out / f"{screen}.synthons.parquet",
                       "--hits", out / "hits.parquet", "--id", screen,
                       "--counts-qc", out / f"{screen}.counts_qc.json",
                       "--output", out / f"{screen}.figures.json") == 0  # fmt: skip
    return out


@pytest.fixture(scope="session")
def screen_inputs(tmp_path_factory: pytest.TempPathFactory, library: SimulatedLibrary, library_dir: Path) -> dict:
    """Every simulated screen on the whole library (15 BB-As), scored and aggregated, ready for triage."""
    from delt_pipeline.analyse.score import exposure_from_qc  # noqa: PLC0415
    from delt_pipeline.analyse.synthons import aggregate_synthons  # noqa: PLC0415
    from delt_pipeline.analyse.triage import ScreenInput  # noqa: PLC0415

    building_blocks = read_building_blocks(library_dir / "synthons")
    whole = select_sublibrary(building_blocks, list(library.heldout.values()), n_bb_a=15, seed=42)
    out = tmp_path_factory.mktemp("triage")
    inputs = {}
    for screen, role in (("mapk14", "target"), ("ddr1", "target"), ("bca", "counter")):
        counts_qc = json.loads(
            json.dumps(
                write_counts(library_dir / f"{screen}.parquet", building_blocks, whole, out / f"{screen}.parquet")
            )
        )
        scores, _ = score_table(pl.read_parquet(out / f"{screen}.parquet"), counts_qc, whole)
        synthons = aggregate_synthons(scores.lazy(), exposure_from_qc(counts_qc))
        inputs[screen] = ScreenInput(screen, role, scores.head(300), scores.lazy(), synthons, counts_qc)
    return inputs
