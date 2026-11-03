import json
from pathlib import Path

import pytest

from delt_pipeline.analyse import analyse
from delt_pipeline.fair import DatasetInfo, Provenance, describe_dataset
from delt_pipeline.io import FeatureTable
from delt_pipeline.report import discover_ids, load_bundle, parse_versions, render_report
from delt_pipeline.visualise import build_figures, figures_to_json
from delt_pipeline.wrangle import wrangle

PIPELINE = {"name": "demo-pipeline", "version": "1.2.3", "homepage": "https://example.org/demo"}


def write_bundle(directory: Path, dataset_id: str, table: FeatureTable, qc: dict, info: DatasetInfo) -> None:
    result = analyse(table, cv_folds=3)
    clean = directory / f"{dataset_id}.clean.parquet"
    table.frame.write_parquet(clean)
    description = describe_dataset(table, clean, [], info, Provenance("demo-pipeline", "1.2.3"))
    (directory / f"{dataset_id}.qc.json").write_text(json.dumps(qc))
    (directory / f"{dataset_id}.stats.json").write_text(json.dumps(result.stats))
    (directory / f"{dataset_id}.figures.json").write_text(
        figures_to_json(build_figures(table, result.embedding, result.stats))
    )
    (directory / f"{dataset_id}.metadata.jsonld").write_text(json.dumps(description.jsonld))


@pytest.fixture
def bundle_dir(tmp_path: Path, raw_frame) -> Path:
    labelled = wrangle(raw_frame, "Sample ID", "label")
    unlabelled = wrangle(raw_frame.drop("label", "constant", "sparse"), "Sample ID")
    rich = DatasetInfo(id="labelled", source="doi:10.1000/x", creator_name="Ada", creator_orcid="0000-0002-1825-0097")
    write_bundle(tmp_path, "labelled", labelled.table, labelled.qc, rich)
    write_bundle(
        tmp_path, "unlabelled", unlabelled.table, unlabelled.qc, DatasetInfo(id="unlabelled", creator_name="Bo")
    )
    return tmp_path


def test_discover_and_load(bundle_dir: Path) -> None:
    assert discover_ids(bundle_dir) == ["labelled", "unlabelled"]
    bundle = load_bundle(bundle_dir, "labelled")
    assert bundle.id == "labelled"
    assert set(bundle.figures) == {"pca", "correlation", "distributions", "importance"}


def test_parse_versions(tmp_path: Path) -> None:
    path = tmp_path / "versions.yml"
    path.write_text("# comment\npython: 3.14\ngzip: '1.13'\n  nested: x\nnot a pair\n")
    assert parse_versions(path) == [("gzip", "1.13"), ("python", "3.14")]


def test_render_report(bundle_dir: Path, tmp_path: Path) -> None:
    bundles = [load_bundle(bundle_dir, dataset_id) for dataset_id in discover_ids(bundle_dir)]
    html = render_report("Demo <report>", bundles, [("gzip", "1.13")], PIPELINE)
    assert "<title>Demo &lt;report&gt;</title>" in html
    assert 'id="ds-labelled"' in html
    assert 'id="ds-unlabelled"' in html
    assert 'id="fig-labelled-importance"' in html
    assert "Plotly.newPlot" in html
    assert "<script src" not in html  # fully offline
    assert "https://doi.org/10.1000/x" in html
    assert "https://orcid.org/0000-0002-1825-0097" in html
    assert "No features dropped." in html
    assert 'href="https://example.org/demo"' in html


def test_render_report_edge_cases(bundle_dir: Path) -> None:
    bundle = load_bundle(bundle_dir, "labelled")
    bundle.stats["model"] = {"skipped": "too few samples"}
    bundle.stats["embedding"]["silhouette"] = None
    bundle.stats["features"][0]["shapiro_p"] = None
    bundle.stats["association"][0]["p_value"] = None
    html = render_report("Single", [bundle], [], {"name": "p", "version": "0", "homepage": ""})
    assert "too few samples" in html
    assert '<div class="label">Silhouette</div>' not in html
    assert "1 dataset profiled" in html
    bundle.stats["model"] = None
    assert "Supervised baseline" not in render_report("x", [bundle], [], PIPELINE)
