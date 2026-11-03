"""Data analysis with scipy and scikit-learn."""

from dataclasses import dataclass
from typing import Any

import polars as pl

from ..io import FeatureTable
from .embed import embed
from .model import evaluate_model
from .stats import association, correlation, describe_features


@dataclass(frozen=True)
class AnalysisResult:
    """Per-sample embedding table and JSON-serialisable statistics."""

    embedding: pl.DataFrame
    stats: dict[str, Any]


def analyse(
    table: FeatureTable,
    *,
    n_components: int = 2,
    n_clusters: int = 3,
    cv_folds: int = 5,
    seed: int = 42,
    threads: int = 1,
) -> AnalysisResult:
    """Run descriptive statistics, the PCA/k-means embedding and the supervised baseline."""
    embedding, embedding_summary = embed(table, n_components, n_clusters, seed)
    stats = {
        "n_samples": table.frame.height,
        "n_features": len(table.features),
        "task": table.task,
        "features": describe_features(table),
        "correlation": correlation(table),
        "association": association(table),
        "embedding": embedding_summary,
        "model": evaluate_model(table, cv_folds, seed, threads),
    }
    return AnalysisResult(embedding, stats)


__all__ = ["AnalysisResult", "analyse", "association", "correlation", "describe_features", "embed", "evaluate_model"]
