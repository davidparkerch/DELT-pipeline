import json
import re
import runpy
import sys
from pathlib import Path

import polars as pl
import pytest
from conftest import run_cli

from delt_pipeline import __version__
from delt_pipeline.cli import main


def test_simulate(tmp_path: Path) -> None:
    assert run_cli("simulate", "--outdir", tmp_path, "--seed", 3) == 0
    assert {path.name for path in tmp_path.glob("*.parquet")} == {"mapk14.parquet", "ddr1.parquet", "bca.parquet"}
    assert len(list((tmp_path / "heldout").glob("*.csv"))) == 4
    assert len(list((tmp_path / "synthons").glob("reacted_*.csv"))) == 3


def test_pipeline_steps(pipeline_dir: Path) -> None:
    sublibrary = json.loads((pipeline_dir / "sublibrary.json").read_text())
    assert len(sublibrary["heldout"]) + len(sublibrary["random"]) == 8
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    assert top.height == 300
    assert json.loads((pipeline_dir / "bca.validation.json").read_text()) == []
    validation = json.loads((pipeline_dir / "mapk14.validation.json").read_text())
    assert [result["condition"] for result in validation] == ["on", "off"]
    assert set(json.loads((pipeline_dir / "mapk14.figures.json").read_text())) == {
        "cube",
        "ma",
        "replicates",
        "saturation",
        "validation",
    }
    assert set(json.loads((pipeline_dir / "bca.figures.json").read_text())) == {
        "cube",
        "ma",
        "replicates",
        "saturation",
    }
    assert {path.name for path in (pipeline_dir / "exports/mapk14").iterdir()} == {
        "mapk14_hits.parquet",
        "mapk14_hits.tsv",
        "datapackage.json",
    }
    metadata = json.loads((pipeline_dir / "mapk14.metadata.jsonld").read_text())
    assert metadata["keywords"] == ["DEL", "kinase"]
    assert metadata["license"] == "https://github.com/insitro/kindel/blob/main/LICENSE.md"
    assert (pipeline_dir / "mapk14.identifier.txt").read_text().strip() == metadata["@id"]
    variables = {variable["name"]: variable for variable in metadata["variableMeasured"]}
    assert variables["poisson_lower"]["description"].startswith("Lower bound")


def test_synthons_and_triage(pipeline_dir: Path) -> None:
    synthons = pl.read_parquet(pipeline_dir / "mapk14.synthons.parquet")
    assert synthons.get_column("level").unique().len() == 6
    assert set(json.loads((pipeline_dir / "mapk14.synthon_qc.json").read_text())) == {"levels", "top"}
    hits = pl.read_parquet(pipeline_dir / "hits.parquet")
    assert hits.get_column("screen").unique().to_list() == ["mapk14"]
    assert hits.height == 300
    qc = json.loads((pipeline_dir / "triage_qc.json").read_text())
    assert set(qc["screens"]) == {"mapk14"}
    assert set(qc["by_score"]["z_faver"]) == {"mapk14"}


def test_triage_rejects_unknown_role(pipeline_dir: Path, tmp_path: Path) -> None:
    args = ["triage", "--screen", "mapk14", "decoy", "--input-dir", pipeline_dir]
    assert run_cli(*args, "--output", tmp_path / "h.parquet", "--qc", tmp_path / "q.json") == 1


def test_validate_without_heldout_skips_reading_scores(tmp_path: Path) -> None:
    output = tmp_path / "v.json"
    missing = tmp_path / "missing"
    assert run_cli("validate", "--scores", missing, "--score-qc", missing, "--target", "x", "--output", output) == 0
    assert json.loads(output.read_text()) == []


def test_validate_one_assay_and_describe_without_identifier(
    pipeline_dir: Path, library_dir: Path, tmp_path: Path
) -> None:
    output = tmp_path / "v.json"
    scored = ["--scores", pipeline_dir / "mapk14.scores.parquet", "--score-qc", pipeline_dir / "mapk14.score_qc.json"]
    offdna = library_dir / "heldout/mapk14_offdna.csv"
    assert run_cli("validate", *scored, "--target", "mapk14", "--offdna", offdna, "--output", output) == 0
    assert [result["condition"] for result in json.loads(output.read_text())] == ["off"]
    jsonld, package = tmp_path / "m.jsonld", tmp_path / "m.json"
    describe = ["describe", "--input", pipeline_dir / "mapk14.top.parquet", "--id-column", "molecule_hash", "--id", "m"]
    describe += ["--pipeline", "demo", "--pipeline-version", "1", "--jsonld", jsonld, "--datapackage", package]
    assert run_cli(*describe) == 0
    assert json.loads(package.read_text())["resources"] == []


def test_report(pipeline_dir: Path, tmp_path: Path) -> None:
    versions, report = tmp_path / "versions.yml", tmp_path / "report.html"
    versions.write_text("demo: 1.0\n")
    args = ["report", "--input-dir", pipeline_dir, "--versions", versions, "--title", "T"]
    assert run_cli(*args, "--pipeline", "demo", "--pipeline-version", "1.0", "--output", report) == 0
    assert "<title>T</title>" in report.read_text()


@pytest.mark.parametrize(
    ("args", "message"),
    [
        (["export", "--input", "{top}", "--id-column", "molecule_hash", "--formats", "xlsx", "--prefix", "x"], "xlsx"),
        (
            [
                "report",
                "--input-dir",
                "{empty}",
                "--versions",
                "v.yml",
                "--pipeline",
                "p",
                "--pipeline-version",
                "1",
                "--output",
                "r.html",
            ],
            "no '\\*.metadata.jsonld'",
        ),
        (
            [
                "score",
                "--input",
                "{counts}",
                "--counts-qc",
                "{counts_qc}",
                "--sublibrary",
                "{sub}",
                "--top-n",
                "0",
                "--output",
                "s",
                "--top",
                "t",
                "--qc",
                "q",
            ],
            "--top-n must be at least 1",
        ),
        (
            [
                "subset",
                "--input",
                "{missing}",
                "--synthons",
                "{missing}",
                "--sublibrary",
                "{sub}",
                "--output",
                "o",
                "--qc",
                "q",
            ],
            "missing",
        ),
    ],
)
def test_errors_exit_nonzero(tmp_path: Path, pipeline_dir: Path, args, message, capsys) -> None:
    (tmp_path / "empty").mkdir()
    paths = {
        "top": pipeline_dir / "mapk14.top.parquet",
        "counts": pipeline_dir / "mapk14.counts.parquet",
        "counts_qc": pipeline_dir / "mapk14.counts_qc.json",
        "sub": pipeline_dir / "sublibrary.json",
        "empty": tmp_path / "empty",
        "missing": tmp_path / "missing",
    }
    assert main([arg.format(**paths) for arg in args]) == 1
    assert re.search(message, capsys.readouterr().err)


def test_version(capsys) -> None:
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == __version__


def test_module_entry_point(monkeypatch, capsys) -> None:
    monkeypatch.setattr(sys, "argv", ["delt_pipeline", "--version"])
    with pytest.raises(SystemExit) as exit_info:
        runpy.run_module("delt_pipeline", run_name="__main__")
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == __version__


def test_dictionary_covers_published_tables(pipeline_dir: Path) -> None:
    """Every column of the published tables (hits per screen and triaged, building blocks, shortlist) is described."""
    dictionary = set(pl.read_csv(Path(__file__).parents[2] / "assets/hits_dictionary.csv").get_column("column"))
    for table in ("mapk14.top.parquet", "hits.parquet", "blocks.parquet", "picks.parquet"):
        assert set(pl.read_parquet(pipeline_dir / table).columns) <= dictionary, table
