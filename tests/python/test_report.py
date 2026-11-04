from pathlib import Path

from delt_pipeline.report import discover_ids, load_bundle, parse_versions, render_report

PIPELINE = {"name": "demo-pipeline", "version": "1.2.3", "homepage": "https://example.org/demo"}


def test_discover_and_load(pipeline_dir: Path) -> None:
    assert discover_ids(pipeline_dir) == ["bca", "mapk14"]
    bundle = load_bundle(pipeline_dir, "mapk14")
    assert bundle.id == "mapk14"
    assert set(bundle.figures) == {"cube", "ma", "replicates", "saturation", "validation"}
    assert [result["condition"] for result in bundle.validation] == ["on", "off"]
    assert load_bundle(pipeline_dir, "bca").validation == []


def test_parse_versions(tmp_path: Path) -> None:
    path = tmp_path / "versions.yml"
    path.write_text("# comment\npython: 3.14\ngzip: '1.13'\n  nested: x\nnot a pair\n")
    assert parse_versions(path) == [("gzip", "1.13"), ("python", "3.14")]


def test_render_report(pipeline_dir: Path) -> None:
    bundles = [load_bundle(pipeline_dir, screen) for screen in discover_ids(pipeline_dir)]
    html = render_report("Demo <report>", bundles, [("polars", "2.1")], PIPELINE)
    assert "<title>Demo &lt;report&gt;</title>" in html
    assert 'id="ds-mapk14"' in html
    assert 'id="ds-bca"' in html
    assert 'id="fig-mapk14-validation"' in html
    assert 'id="fig-bca-cube"' in html
    assert "Plotly.newPlot" in html
    assert "<script src" not in html  # fully offline
    assert "(0.737)" in html  # KinDEL's published Poisson baseline next to ours
    assert "https://github.com/insitro/kindel/blob/main/LICENSE.md" in html
    assert "https://orcid.org/0000-0002-1825-0097" in html
    assert 'href="https://example.org/demo"' in html
    assert "2 DEL selections" in html
    assert "Data license:" in html  # both screens share the KinDEL license
    assert 'id="triage"' not in html
    assert "Replicate r" in html


def test_render_report_with_triage(pipeline_dir: Path, tmp_path: Path) -> None:
    from conftest import run_cli  # noqa: PLC0415

    versions, report = tmp_path / "versions.yml", tmp_path / "report.html"
    versions.write_text("demo: 1.0\n")
    args = ["report", "--input-dir", pipeline_dir, "--versions", versions, "--pipeline", "demo"]
    assert run_cli(*args, "--pipeline-version", "1.0", "--output", report) == 0
    html = report.read_text()
    assert 'id="triage"' in html
    assert "mapk14: best tier-1 hits" in html
    assert "replicates; disynthon" in html
    assert "What each score lets through" in html
    assert "Normalised z-score (Faver 2019)" in html


def test_render_report_without_validation(pipeline_dir: Path) -> None:
    bundle = load_bundle(pipeline_dir, "bca")
    html = render_report("Single", [bundle], [], {"name": "p", "version": "0", "homepage": ""})
    assert "No screen has held-out KD data." in html
    assert "1 DEL selection scored" in html
    flags = {"counter_binder": 0, "matrix_binder": 0, "shared_target": 0}
    screen = {"n_compounds": 1, "tiers": dict.fromkeys(["1", "2", "3", "flagged"], 0), "flags": flags, "top_tier1": []}
    thresholds = {"min_replicates": 2, "disynthon_lower": 1.5, "counter_lower": 1, "matrix_lower": 4}
    empty = {"screens": {"x": screen}, "thresholds": thresholds, "by_score": {}}
    pipeline = {"name": "p", "version": "0", "homepage": ""}
    assert "No tier-1 hits." in render_report("T", [bundle], [], pipeline, empty)
