"""Interactive figures with plotly, serialised as deterministic JSON."""

import json
from typing import Any

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import polars as pl

from .._numbers import significant
from ..io import FeatureTable

TEMPLATE = "plotly_white"


def _layout(title: str, **kwargs: object) -> dict[str, object]:
    return {"title": {"text": title}, "template": TEMPLATE, "margin": {"l": 60, "r": 20, "t": 60, "b": 60}, **kwargs}


def pca_scatter(table: FeatureTable, embedding: pl.DataFrame, stats: dict[str, Any]) -> go.Figure:
    """Samples on the first two principal components, coloured by target (or cluster)."""
    ratios = stats["embedding"]["explained_variance_ratio"]
    x = embedding.get_column("PC1").to_list()
    has_pc2 = "PC2" in embedding.columns
    y = embedding.get_column("PC2").to_list() if has_pc2 else [0.0] * len(x)
    colour_by = table.target_column if table.target_column is not None else "cluster"
    groups = embedding.get_column(colour_by).cast(pl.String).to_list()
    ids = embedding.get_column(table.id_column).to_list()
    figure = go.Figure()
    for group in sorted(set(groups)):
        members = [index for index, value in enumerate(groups) if value == group]
        figure.add_trace(
            go.Scatter(
                x=[x[i] for i in members],
                y=[y[i] for i in members],
                text=[ids[i] for i in members],
                mode="markers",
                name=f"{colour_by}={group}",
                hovertemplate="%{text}<br>PC1=%{x}<br>PC2=%{y}<extra></extra>",
            )
        )
    pc2_title = f"PC2 ({ratios[1]:.1%})" if has_pc2 else "PC2 (n/a)"
    figure.update_layout(
        _layout("PCA of standardised features", xaxis={"title": f"PC1 ({ratios[0]:.1%})"}, yaxis={"title": pc2_title})
    )
    return figure


def correlation_heatmap(stats: dict[str, Any]) -> go.Figure:
    """Spearman correlation matrix between features."""
    correlation = stats["correlation"]
    figure = go.Figure(
        go.Heatmap(
            z=correlation["matrix"],
            x=correlation["features"],
            y=correlation["features"],
            zmin=-1,
            zmax=1,
            colorscale="RdBu",
            reversescale=True,
            colorbar={"title": {"text": "rho"}},
        )
    )
    figure.update_layout(_layout("Spearman correlation", yaxis={"autorange": "reversed"}))
    return figure


def feature_distributions(table: FeatureTable) -> go.Figure:
    """Box plots of z-scored features so different units share one axis."""
    matrix = table.matrix()
    spread = matrix.std(axis=0)
    scaled = (matrix - matrix.mean(axis=0)) / np.where(spread > 0, spread, 1.0)
    figure = go.Figure()
    for index, name in enumerate(table.features):
        figure.add_trace(go.Box(y=[significant(value) for value in scaled[:, index]], name=name, boxmean=True))
    figure.update_layout(_layout("Feature distributions (z-scores)", showlegend=False, yaxis={"title": "z-score"}))
    return figure


def importance_bar(model: dict[str, Any]) -> go.Figure:
    """Random-forest feature importances, most important first."""
    ranked = model["feature_importance"]
    figure = go.Figure(
        go.Bar(x=[item["importance"] for item in ranked], y=[item["feature"] for item in ranked], orientation="h")
    )
    figure.update_layout(
        _layout(
            f"Feature importance ({model['estimator']})", yaxis={"autorange": "reversed"}, xaxis={"title": "importance"}
        )
    )
    return figure


def build_figures(table: FeatureTable, embedding: pl.DataFrame, stats: dict[str, Any]) -> dict[str, go.Figure]:
    """Build every figure that the data supports, keyed by a stable name."""
    figures = {
        "pca": pca_scatter(table, embedding, stats),
        "correlation": correlation_heatmap(stats),
        "distributions": feature_distributions(table),
    }
    model = stats.get("model")
    if model and model.get("feature_importance"):
        figures["importance"] = importance_bar(model)
    return figures


def figures_to_json(figures: dict[str, go.Figure]) -> str:
    """Serialise figures to one JSON document with sorted keys (byte-stable)."""
    payload = {name: json.loads(pio.to_json(figure, engine="json")) for name, figure in figures.items()}
    return json.dumps(payload, sort_keys=True, indent=1) + "\n"
