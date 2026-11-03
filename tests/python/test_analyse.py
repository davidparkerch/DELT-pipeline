import polars as pl
import pytest

from delt_pipeline.analyse import analyse, association, correlation, describe_features, embed, evaluate_model
from delt_pipeline.io import FeatureTable


def test_describe_features(table: FeatureTable) -> None:
    summaries = describe_features(table)
    assert [item["feature"] for item in summaries] == table.features
    assert all(item["shapiro_p"] is not None for item in summaries)


def test_describe_degenerate_features() -> None:
    single = FeatureTable(pl.DataFrame({"id": ["a"], "x": [1.0]}), "id")
    flat = FeatureTable(pl.DataFrame({"id": ["a", "b", "c"], "x": [2.0, 2.0, 2.0]}), "id")
    assert describe_features(single)[0]["std"] is None
    assert describe_features(flat)[0]["skewness"] is None
    assert describe_features(flat)[0]["shapiro_p"] is None


@pytest.mark.parametrize("n_features", [1, 2, 3])
def test_correlation_shapes(n_features: int) -> None:
    columns = {f"x{index}": [1.0, 2.0, 3.0, float(4 + index)] for index in range(n_features)}
    table = FeatureTable(pl.DataFrame({"id": ["a", "b", "c", "d"], **columns}), "id")
    result = correlation(table)
    assert len(result["matrix"]) == n_features
    assert all(len(row) == n_features for row in result["matrix"])
    assert result["matrix"][0][0] == 1.0


def test_association(table: FeatureTable, regression: FeatureTable, unlabelled: FeatureTable) -> None:
    assert {item["test"] for item in association(table)} == {"kruskal"}
    assert {item["test"] for item in association(regression)} == {"spearman"}
    assert association(unlabelled) == []


def test_association_identical_values() -> None:
    frame = pl.DataFrame({"id": ["a", "b", "c", "d"], "x": [1.0, 1.0, 1.0, 1.0], "y": ["p", "q", "p", "q"]})
    result = association(FeatureTable(frame, "id", "y"))
    assert result[0]["p_value"] is None


def test_embed_clips_parameters(table: FeatureTable) -> None:
    embedding, summary = embed(table, n_components=10, n_clusters=2, seed=1)
    assert embedding.columns == ["Sample_ID", "PC1", "PC2", "PC3", "cluster", "label"]
    assert summary["n_components"] == 3
    assert sum(summary["cluster_sizes"].values()) == 24
    assert summary["silhouette"] is not None


def test_embed_single_cluster(unlabelled: FeatureTable) -> None:
    embedding, summary = embed(unlabelled, n_components=2, n_clusters=1, seed=1)
    assert "label" not in embedding.columns
    assert summary["silhouette"] is None


def test_model_binary_classification(table: FeatureTable) -> None:
    result = evaluate_model(table, cv_folds=3, seed=1)
    assert result["estimator"] == "RandomForestClassifier"
    assert set(result["metrics"]) == {"accuracy", "balanced_accuracy", "f1_macro", "roc_auc"}
    assert result["metrics"]["accuracy"]["mean"] > 0.8
    assert len(result["feature_importance"]) == 3


def test_model_multiclass() -> None:
    labels = ["a", "b", "c"] * 4
    frame = pl.DataFrame({"id": [str(i) for i in range(12)], "x": [float("abc".index(v)) for v in labels], "y": labels})
    result = evaluate_model(FeatureTable(frame, "id", "y"), cv_folds=3, seed=1)
    assert "roc_auc" not in result["metrics"]


def test_model_regression(regression: FeatureTable) -> None:
    result = evaluate_model(regression, cv_folds=4, seed=1)
    assert set(result["metrics"]) == {"r2", "mae", "rmse"}
    assert result["metrics"]["mae"]["mean"] > 0
    assert result["feature_importance"][0]["feature"] == "x1"


def test_model_skips_without_enough_data() -> None:
    frame = pl.DataFrame({"id": ["a", "b", "c"], "x": [1.0, 2.0, 3.0], "y": ["p", "p", "q"]})
    result = evaluate_model(FeatureTable(frame, "id", "y"), cv_folds=5, seed=1)
    assert "skipped" in result
    assert "metrics" not in result


def test_model_without_target(unlabelled: FeatureTable) -> None:
    assert evaluate_model(unlabelled, cv_folds=3, seed=1) is None


def test_analyse(table: FeatureTable) -> None:
    result = analyse(table, n_components=2, n_clusters=2, cv_folds=3, seed=42)
    assert result.embedding.height == 24
    assert result.stats["task"] == "classification"
    assert set(result.stats) == {
        "n_samples",
        "n_features",
        "task",
        "features",
        "correlation",
        "association",
        "embedding",
        "model",
    }
