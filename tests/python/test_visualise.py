import json

from delt_pipeline.analyse import analyse
from delt_pipeline.io import FeatureTable
from delt_pipeline.visualise import build_figures, figures_to_json


def test_figures_with_target_and_model(table: FeatureTable) -> None:
    result = analyse(table, n_components=2, n_clusters=2, cv_folds=3)
    figures = build_figures(table, result.embedding, result.stats)
    assert list(figures) == ["pca", "correlation", "distributions", "importance"]
    assert {trace.name for trace in figures["pca"].data} == {"label=active", "label=inactive"}


def test_figures_without_target(unlabelled: FeatureTable) -> None:
    result = analyse(unlabelled, n_components=1, n_clusters=2)
    figures = build_figures(unlabelled, result.embedding, result.stats)
    assert "importance" not in figures
    assert figures["pca"].layout.yaxis.title.text == "PC2 (n/a)"
    assert all(trace.name.startswith("cluster=") for trace in figures["pca"].data)


def test_figures_json_is_deterministic(table: FeatureTable) -> None:
    result = analyse(table, cv_folds=3)
    first = figures_to_json(build_figures(table, result.embedding, result.stats))
    second = figures_to_json(build_figures(table, result.embedding, result.stats))
    assert first == second
    assert list(json.loads(first)) == sorted(json.loads(first))
