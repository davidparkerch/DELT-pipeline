import json
import re
import runpy
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from delt_pipeline import __version__
from delt_pipeline.cli import main


@pytest.fixture
def raw_csv(tmp_path: Path, raw_frame: pl.DataFrame) -> Path:
    path = tmp_path / "raw.csv"
    raw_frame.write_csv(path)
    return path


def run(*args: str | Path) -> int:
    return main([str(arg) for arg in args])


def test_full_pipeline(tmp_path: Path, raw_csv: Path) -> None:
    table = ["--id-column", "Sample ID", "--target-column", "label"]
    clean, qc = tmp_path / "ds.clean.parquet", tmp_path / "ds.qc.json"
    embedding, stats = tmp_path / "ds.embedding.parquet", tmp_path / "ds.stats.json"
    figures, exports = tmp_path / "ds.figures.json", tmp_path / "exports"

    assert run("wrangle", "--input", raw_csv, *table, "--output", clean, "--qc", qc) == 0
    assert json.loads(qc.read_text())["rows_out"] == 24
    assert run("analyse", "--input", clean, *table, "--cv-folds", "3", "--embedding", embedding, "--stats", stats) == 0
    assert (
        run("visualise", "--input", clean, *table, "--embedding", embedding, "--stats", stats, "--output", figures) == 0
    )
    assert "pca" in json.loads(figures.read_text())
    assert run("export", "--input", clean, *table, "--formats", "csv,npz", "--prefix", "ds", "--outdir", exports) == 0
    assert np.load(exports / "ds.npz")["X"].shape == (24, 3)

    sums = tmp_path / "SHA256SUMS"
    sums.write_text(f"{'0' * 64}  ds.csv\n")
    dictionary = tmp_path / "dictionary.csv"
    dictionary.write_text("column,description\nfeat alpha,Alpha\n")
    jsonld, package, identifier = tmp_path / "ds.metadata.jsonld", tmp_path / "datapackage.json", tmp_path / "id.txt"
    describe = ["describe", "--input", clean, *table, "--id", "ds", "--pipeline", "demo", "--pipeline-version", "1.0"]
    describe += ["--exports", exports / "ds.csv", exports / "ds.npz", "--sha256sums", sums, "--stats", stats]
    describe += ["--dictionary", dictionary, "--keywords", "chem; ml,", "--exports-base", "../exports/"]
    describe += ["--jsonld", jsonld, "--datapackage", package, "--identifier", identifier]
    assert run(*describe) == 0
    metadata = json.loads(jsonld.read_text())
    assert metadata["keywords"] == ["chem", "ml"]
    assert identifier.read_text().strip() == metadata["@id"]
    assert {item["contentUrl"] for item in metadata["distribution"]} >= {"../exports/ds.csv", "ds.clean.parquet"}
    assert json.loads(package.read_text())["resources"][0]["hash"] == f"sha256:{'0' * 64}"

    versions, report = tmp_path / "versions.yml", tmp_path / "report.html"
    versions.write_text("demo: 1.0\n")
    assert (
        run(
            "report",
            "--input-dir",
            tmp_path,
            "--versions",
            versions,
            "--title",
            "T",
            "--pipeline",
            "demo",
            "--pipeline-version",
            "1.0",
            "--output",
            report,
        )
        == 0
    )
    assert "<title>T</title>" in report.read_text()


def test_describe_minimal(tmp_path: Path, raw_csv: Path) -> None:
    clean, qc = tmp_path / "u.clean.parquet", tmp_path / "u.qc.json"
    run("wrangle", "--input", raw_csv, "--id-column", "Sample ID", "--output", clean, "--qc", qc)
    jsonld, package = tmp_path / "u.jsonld", tmp_path / "u.json"
    args = ["describe", "--input", clean, "--id-column", "Sample ID", "--id", "u", "--pipeline", "demo"]
    assert run(*args, "--pipeline-version", "1", "--jsonld", jsonld, "--datapackage", package) == 0
    assert json.loads(package.read_text())["resources"] == []


@pytest.mark.parametrize(
    ("args", "message"),
    [
        (["export", "--input", "{clean}", "--id-column", "Sample ID", "--formats", "xlsx", "--prefix", "x"], "xlsx"),
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
        (["wrangle", "--input", "{missing}", "--id-column", "x", "--output", "o", "--qc", "q"], "missing"),
    ],
)
def test_errors_exit_nonzero(tmp_path: Path, raw_csv, args, message, capsys) -> None:
    clean = tmp_path / "c.parquet"
    run("wrangle", "--input", raw_csv, "--id-column", "Sample ID", "--output", clean, "--qc", tmp_path / "q.json")
    (tmp_path / "empty").mkdir()
    paths = {"clean": clean, "empty": tmp_path / "empty", "missing": tmp_path / "missing.csv"}
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
