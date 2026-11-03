"""Unsupervised structure: PCA embedding and k-means clustering with scikit-learn."""

from typing import Any

import polars as pl
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

from .._numbers import significant
from ..io import FeatureTable


def embed(table: FeatureTable, n_components: int, n_clusters: int, seed: int) -> tuple[pl.DataFrame, dict[str, Any]]:
    """Project standardised features with PCA and cluster samples with k-means.

    ``n_components`` and ``n_clusters`` are clipped to what the data supports.

    Returns:
        A table with the identifier, ``PC1..PCn``, ``cluster`` and (if present) the target,
        plus a summary with explained variance and the silhouette score.
    """
    scaled = StandardScaler().fit_transform(table.matrix())
    n_samples, n_features = scaled.shape
    components = max(1, min(n_components, n_features, n_samples))
    clusters = max(1, min(n_clusters, n_samples))

    pca = PCA(n_components=components, random_state=seed)
    coordinates = pca.fit_transform(scaled)
    labels = KMeans(n_clusters=clusters, n_init=10, random_state=seed).fit_predict(scaled)

    silhouette = None
    if 1 < len(set(labels)) < n_samples:
        silhouette = significant(silhouette_score(scaled, labels))

    columns: dict[str, Any] = {table.id_column: table.ids}
    for index in range(components):
        columns[f"PC{index + 1}"] = [significant(value) for value in coordinates[:, index]]
    columns["cluster"] = pl.Series(labels, dtype=pl.Int32)
    if table.target_column is not None:
        columns[table.target_column] = table.frame.get_column(table.target_column)
    summary = {
        "n_components": components,
        "explained_variance_ratio": [significant(value) for value in pca.explained_variance_ratio_],
        "n_clusters": clusters,
        "cluster_sizes": {str(label): int((labels == label).sum()) for label in sorted(set(labels))},
        "silhouette": silhouette,
    }
    return pl.DataFrame(columns), summary
