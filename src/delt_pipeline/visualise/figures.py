"""Interactive figures with plotly, serialised as deterministic JSON."""

import json
from typing import Any

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import polars as pl

from ..library.schema import BB_COLUMNS, ID

TEMPLATE = "plotly_white"
COLOURSCALE = "Viridis"
ON_DNA, OFF_DNA, PUBLISHED = "#2f6fde", "#e07b28", "#5c6577"


def _layout(title: str, **kwargs: object) -> dict[str, object]:
    return {"title": {"text": title}, "template": TEMPLATE, "margin": {"l": 60, "r": 20, "t": 60, "b": 60}, **kwargs}


def cube(top: pl.DataFrame, screen: str) -> go.Figure:
    """The DEL "cube": top compounds placed by their three building blocks.

    A line of points is a disynthon (two building blocks bind, the third varies); a single bright point
    with no neighbours is weaker evidence than a whole line.
    """
    score = np.log2(np.clip(top.get_column("poisson_lower").fill_null(0.0).to_numpy(), 1e-3, None))
    text = [
        f"{molecule[:12]}…<br>target {target} / control {control}"
        for molecule, target, control in top.select(ID, "target_sum", "control_sum").iter_rows()
    ]
    trace = go.Scatter3d(
        x=top.get_column("bb_a").to_list(),
        y=top.get_column("bb_b").to_list(),
        z=top.get_column("bb_c").to_list(),
        mode="markers",
        text=text,
        marker={
            "size": 3,
            "color": score.round(3).tolist(),
            "colorscale": COLOURSCALE,
            "colorbar": {"title": {"text": "log2 lower bound"}},
            "opacity": 0.8,
        },
        hovertemplate="A %{x} · B %{y} · C %{z}<br>%{text}<extra></extra>",
    )
    axes = {
        f"{axis}axis": {"title": {"text": f"Building block {position.upper()}"}}
        for axis, position in zip("xyz", (name[-1] for name in BB_COLUMNS), strict=True)
    }
    return go.Figure(
        trace,
        layout=_layout(f"{screen}: top {top.height} compounds by Poisson lower bound", scene=axes, height=560),
    )


SHORT_LABELS = {
    "count_diff": "Counts",
    "er": "Enrichment ratio",
    "poisson_lower": "Poisson lower bound",
    "z_faver": "Faver z",
    "eb_log2_ratio": "EB shrunken ratio",
}


def validation_bars(validations: list[dict[str, Any]], screen: str) -> go.Figure:
    """-Spearman (score vs KD) of every score, on- and off-DNA, next to KinDEL's published baselines."""
    figure = go.Figure()
    for result in validations:
        condition = result["condition"]
        rows = result["scores"]
        colour = ON_DNA if condition == "on" else OFF_DNA
        labels = [SHORT_LABELS[row["score"]] for row in rows]
        figure.add_trace(
            go.Bar(
                x=labels,
                y=[row["spearman"] for row in rows],
                name=f"{condition}-DNA, this pipeline (n={result['n_compounds']})",
                marker_color=colour,
            )
        )
        figure.add_trace(
            go.Bar(
                x=labels,
                y=[row["published"] for row in rows],
                name=f"{condition}-DNA, KinDEL published",
                marker={"color": "rgba(0,0,0,0)", "line": {"color": colour, "width": 2}, "pattern": {"shape": "/"}},
                marker_pattern_fgcolor=colour,
            )
        )
    figure.update_layout(
        _layout(
            f"{screen}: -Spearman of score vs measured KD (higher is better)",
            barmode="group",
            yaxis={"title": {"text": "-Spearman rho"}, "range": [-0.2, 1], "zeroline": True},
            legend={"orientation": "h", "y": -0.18},
            height=480,
        )
    )
    return figure


def build_figures(top: pl.DataFrame, validations: list[dict[str, Any]], screen: str) -> dict[str, go.Figure]:
    """Every figure the data supports, keyed by a stable name."""
    figures = {"cube": cube(top, screen)}
    if validations:
        figures["validation"] = validation_bars(validations, screen)
    return figures


def figures_to_json(figures: dict[str, go.Figure]) -> str:
    """Serialise figures to one JSON document with sorted keys (byte-stable)."""
    payload = {name: json.loads(pio.to_json(figure, engine="json")) for name, figure in figures.items()}
    return json.dumps(payload, sort_keys=True, indent=1) + "\n"
