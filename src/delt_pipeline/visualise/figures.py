"""Interactive figures with plotly, serialised as deterministic JSON."""

import json
from typing import Any

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import polars as pl
from plotly.subplots import make_subplots

from .._numbers import significant
from ..library.depth import RECOMMENDED_OVERSAMPLING
from ..library.schema import BB_COLUMNS, ID, TARGET

TEMPLATE = "plotly_white"
COLOURSCALE = "Viridis"
ON_DNA, OFF_DNA, PUBLISHED = "#2f6fde", "#e07b28", "#5c6577"

TIER_COLOURS = {"1": "#2a78d6", "2": "#1baf7a", "3": "#9a9890", "flagged": "#eb6834"}
"""Triage tiers on the cube: three categorical hues plus a muted grey for tier 3 (all pairs distinguishable)."""
TIER_LABELS = {"1": "Tier 1", "2": "Tier 2", "3": "Tier 3", "flagged": "Flagged (counter or bead)"}
POINTS, HIGHLIGHT = "#2a78d6", "#eb6834"

SAMPLE = 5000
"""Points drawn in the replicate and MA plots (a seeded sample keeps the report small)."""


def _layout(title: str, **kwargs: object) -> dict[str, object]:
    return {"title": {"text": title}, "template": TEMPLATE, "margin": {"l": 60, "r": 20, "t": 60, "b": 60}, **kwargs}


def _cube_trace(top: pl.DataFrame, **kwargs: object) -> go.Scatter3d:
    text = [
        f"{molecule[:12]}…<br>target {target} / control {control}"
        for molecule, target, control in top.select(ID, "target_sum", "control_sum").iter_rows()
    ]
    return go.Scatter3d(
        x=top.get_column("bb_a").to_list(),
        y=top.get_column("bb_b").to_list(),
        z=top.get_column("bb_c").to_list(),
        mode="markers",
        text=text,
        hovertemplate="A %{x} · B %{y} · C %{z}<br>%{text}<extra></extra>",
        **kwargs,
    )


def cube(top: pl.DataFrame, screen: str, tiers: pl.DataFrame | None = None) -> go.Figure:
    """The DEL "cube": top compounds placed by their three building blocks.

    A line of points is a disynthon (two building blocks bind, the third varies); a single bright point
    with no neighbours is weaker evidence than a whole line. With ``tiers`` (``molecule_hash``, ``tier``)
    points are coloured by triage tier, otherwise by the Poisson lower bound.
    """
    if tiers is None or tiers.is_empty():
        score = np.log2(np.clip(top.get_column("poisson_lower").fill_null(0.0).to_numpy(), 1e-3, None))
        marker = {
            "size": 3,
            "color": score.round(3).tolist(),
            "colorscale": COLOURSCALE,
            "colorbar": {"title": {"text": "log2 lower bound"}},
            "opacity": 0.8,
        }
        traces = [_cube_trace(top, marker=marker, showlegend=False)]
    else:
        tiered = top.join(tiers.select(ID, "tier"), on=ID, how="left", maintain_order="left")
        traces = [
            _cube_trace(
                tiered.filter(pl.col("tier") == tier),
                name=label,
                marker={"size": 3, "color": TIER_COLOURS[tier], "opacity": 0.85},
            )
            for tier, label in TIER_LABELS.items()
        ]
    axes = {
        f"{axis}axis": {"title": {"text": f"Building block {position.upper()}"}}
        for axis, position in zip("xyz", (name[-1] for name in BB_COLUMNS), strict=True)
    }
    return go.Figure(
        traces,
        layout=_layout(
            f"{screen}: top {top.height} compounds by Poisson lower bound",
            scene=axes,
            height=560,
            legend={"itemsizing": "constant"},
        ),
    )


def sample_rows(frame: pl.DataFrame, n: int = SAMPLE, seed: int = 0) -> pl.DataFrame:
    """At most ``n`` rows, drawn with a fixed seed and kept in their original order (byte-stable)."""
    if frame.height <= n:
        return frame
    keep = np.sort(np.random.default_rng(seed).choice(frame.height, size=n, replace=False))
    return frame[keep]


def _log_counts(frame: pl.DataFrame) -> tuple[list[float], list[float], float | None]:
    x = np.log1p(frame.get_column(TARGET[0]).to_numpy().astype(float))
    y = np.log1p(frame.get_column(TARGET[1]).to_numpy().astype(float))
    constant = x.size < 2 or x.std() == 0 or y.std() == 0  # noqa: PLR2004 - correlation needs two points
    r = None if constant else significant(np.corrcoef(x, y)[0, 1], 2)
    return x.round(3).tolist(), y.round(3).tolist(), r


def replicate_agreement(compounds: pl.DataFrame, disynthons: pl.DataFrame, screen: str) -> go.Figure:
    """Target replicate 1 vs 2 per compound and per A-B disynthon: noise at one level, signal at the other.

    Both panels show log(1 + count) of entities with target counts; the correlation is computed on the
    points shown.
    """
    panels = []
    for name, frame in (("compound", compounds), ("A-B disynthon", disynthons)):
        x, y, r = _log_counts(frame.filter(pl.sum_horizontal(TARGET) > 0))
        panels.append((name, x, y, r))
    figure = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=[
            f"Per {name} (r = {r if r is not None else 'n/a'}, n = {len(x):,})" for name, x, _, r in panels
        ],
    )
    for column, (name, x, y, _) in enumerate(panels, start=1):
        figure.add_trace(
            go.Scattergl(
                x=x,
                y=y,
                mode="markers",
                name=name,
                showlegend=False,
                marker={"size": 4, "color": POINTS, "opacity": 0.35},
                hovertemplate="replicate 1 %{x:.2f} · replicate 2 %{y:.2f}<extra></extra>",
            ),
            row=1,
            col=column,
        )
        figure.update_xaxes(title_text="log(1 + count), target replicate 1", row=1, col=column)
        figure.update_yaxes(title_text="log(1 + count), target replicate 2", row=1, col=column)
    figure.update_layout(_layout(f"{screen}: replicates agree on synthons, not on single compounds", height=460))
    return figure


def ma_plot(compounds: pl.DataFrame, top: pl.DataFrame, screen: str) -> go.Figure:
    """Enrichment ratio against mean count (an MA plot): extreme ratios live at low counts.

    The top compounds by Poisson lower bound are highlighted: they combine a high ratio with enough
    counts to trust it.
    """
    figure = go.Figure()
    for name, frame, colour, size in (
        ("Sampled compounds with counts", compounds, POINTS, 4),
        (f"Top {top.height} by Poisson lower bound", top, HIGHLIGHT, 6),
    ):
        mean = (frame.get_column("target_sum") + frame.get_column("control_sum")).to_numpy() / 2.0
        ratio = np.log2(frame.get_column("er").to_numpy())
        figure.add_trace(
            go.Scattergl(
                x=np.log10(mean).round(3).tolist(),
                y=ratio.round(3).tolist(),
                mode="markers",
                name=name,
                marker={"size": size, "color": colour, "opacity": 0.5},
                hovertemplate="log10 mean count %{x:.2f} · log2 ER %{y:.2f}<extra></extra>",
            )
        )
    figure.update_layout(
        _layout(
            f"{screen}: enrichment ratio against mean count",
            xaxis={"title": {"text": "log10 mean count per arm (target, control)"}},
            yaxis={"title": {"text": "log2 enrichment ratio"}},
            legend={"orientation": "h", "y": -0.2},
            height=460,
        )
    )
    return figure


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


SATURATION_COLUMNS = (
    ("target", "Target (pooled)", POINTS),
    ("matrix", "Matrix control (pooled)", PUBLISHED),
    ("load", "Pre-selection load", HIGHLIGHT),
)
"""Count columns drawn in the saturation plot: key in the QC file, label and colour."""


def saturation(counts_qc: dict[str, Any], screen: str) -> go.Figure:
    """How much of the library the counts sampled: rarefaction curves and the target count histogram.

    Left, the fraction of the library observed against counts per library member, for the data as
    sequenced (each curve ends at the real depth), a perfectly uniform library, and the tenfold
    oversampling the protocol literature recommends. Right, how many compounds have 0, 1, 2, ... counts.
    """
    columns = counts_qc["library"]["saturation"]
    figure = make_subplots(
        rows=1, cols=2, subplot_titles=("Library coverage vs. depth", "Counts per compound (target)")
    )
    for key, label, colour in SATURATION_COLUMNS:
        column = columns[key]
        depth = column["oversampling"]
        figure.add_trace(
            go.Scatter(
                x=[significant(f * depth) for f in column["rarefaction"]["fraction"]],
                y=column["rarefaction"]["observed"],
                mode="lines+markers",
                name=f"{label}: {depth:.2g} counts/compound",
                line={"color": colour},
            ),
            row=1,
            col=1,
        )
    grid = np.logspace(-3, np.log10(2 * RECOMMENDED_OVERSAMPLING), 60)
    figure.add_trace(
        go.Scatter(
            x=[significant(x) for x in grid],
            y=[significant(1 - np.exp(-x)) for x in grid],
            mode="lines",
            name="Uniform library (Poisson)",
            line={"color": "#9a9890", "dash": "dot"},
        ),
        row=1,
        col=1,
    )
    figure.add_vline(x=RECOMMENDED_OVERSAMPLING, line={"color": "#9a9890", "dash": "dash"}, row=1, col=1)
    figure.add_annotation(  # domain coordinates: annotations on log axes would need log10 units
        xref="x domain",
        yref="y domain",
        x=0.95,
        y=0.22,
        xanchor="right",
        text="10x oversampling<br>recommended",
        showarrow=False,
        font={"size": 11, "color": "#5c6577"},
    )
    histogram = columns["target"]["histogram"]
    labels = [str(k) for k in range(len(histogram) - 1)] + [f"{len(histogram) - 1}+"]
    figure.add_trace(
        go.Bar(x=labels, y=histogram, marker={"color": POINTS}, name="Compounds", showlegend=False), row=1, col=2
    )
    # Log axes need an explicit range (in log10 units) and decade ticks, or plotly stretches them to fit the shapes.
    figure.update_xaxes(
        type="log",
        range=[-3, float(np.log10(2 * RECOMMENDED_OVERSAMPLING))],
        dtick=1,
        title_text="Counts per library member",
        row=1,
        col=1,
    )
    figure.update_yaxes(title_text="Fraction of library observed", row=1, col=1)
    figure.update_xaxes(title_text="Pooled target counts", row=1, col=2)
    figure.update_yaxes(type="log", dtick=1, title_text="Compounds", row=1, col=2)
    figure.update_layout(
        **_layout(f"{screen}: sequencing depth and library coverage", legend={"orientation": "h", "x": 0, "y": -0.3})
    )
    return figure


def build_figures(
    top: pl.DataFrame,
    validations: list[dict[str, Any]],
    screen: str,
    *,
    compounds: pl.DataFrame | None = None,
    synthons: pl.DataFrame | None = None,
    tiers: pl.DataFrame | None = None,
    counts_qc: dict[str, Any] | None = None,
) -> dict[str, go.Figure]:
    """Every figure the data supports, keyed by a stable name.

    Args:
        top: The screen's top compounds.
        validations: From :func:`~delt_pipeline.analyse.validate.validate_scores` (may be empty).
        screen: Name shown in titles.
        compounds: A sample of compounds with counts (replicate and MA plots).
        synthons: The screen's synthon table (disynthon replicate agreement).
        tiers: ``molecule_hash`` and ``tier`` of the screen's triaged hits (cube colours).
        counts_qc: The screen's count QC (saturation plot).
    """
    figures = {"cube": cube(top, screen, tiers)}
    if compounds is not None:
        figures["ma"] = ma_plot(compounds, top, screen)
        if synthons is not None:
            disynthons = sample_rows(synthons.filter(pl.col("level") == "di_ab"))
            figures["replicates"] = replicate_agreement(compounds, disynthons, screen)
    if counts_qc is not None:
        figures["saturation"] = saturation(counts_qc, screen)
    if validations:
        figures["validation"] = validation_bars(validations, screen)
    return figures


def figures_to_json(figures: dict[str, go.Figure]) -> str:
    """Serialise figures to one JSON document with sorted keys (byte-stable)."""
    payload = {name: json.loads(pio.to_json(figure, engine="json")) for name, figure in figures.items()}
    return json.dumps(payload, sort_keys=True, indent=1) + "\n"
