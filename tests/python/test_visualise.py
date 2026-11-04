import json
from pathlib import Path

import polars as pl
import pytest

from delt_pipeline.visualise import (
    build_figures,
    cube,
    figures_to_json,
    ma_plot,
    replicate_agreement,
    sample_rows,
    saturation,
    validation_bars,
)


def test_figures(pipeline_dir: Path) -> None:
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    validation = json.loads((pipeline_dir / "mapk14.validation.json").read_text())
    figures = build_figures(top, validation, "MAPK14")
    assert set(figures) == {"cube", "validation"}
    assert len(figures["cube"].data[0].x) == top.height
    assert len(validation_bars(validation, "MAPK14").data) == 4  # bars + published markers, on and off
    assert set(build_figures(top, [], "BCA")) == {"cube"}
    assert cube(top, "MAPK14").layout.scene.xaxis.title.text == "Building block A"


def test_cube_by_tier(pipeline_dir: Path) -> None:
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    hits = pl.read_parquet(pipeline_dir / "hits.parquet")
    figure = cube(top, "MAPK14", hits.select("molecule_hash", "tier"))
    assert [trace.name for trace in figure.data] == ["Tier 1", "Tier 2", "Tier 3", "Flagged (counter or bead)"]
    assert sum(len(trace.x) for trace in figure.data) == top.height
    # An empty tier table (a counter-screen has no triaged hits) falls back to the score colours.
    assert cube(top, "BCA", hits.head(0)).data[0].marker.colorscale is not None


def test_replicates_and_ma(pipeline_dir: Path) -> None:
    scores = pl.read_parquet(pipeline_dir / "mapk14.scores.parquet")
    synthons = pl.read_parquet(pipeline_dir / "mapk14.synthons.parquet")
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    compounds = scores.filter(pl.col("target_sum") + pl.col("control_sum") > 0)
    figure = replicate_agreement(compounds, synthons.filter(pl.col("level") == "di_ab"), "MAPK14")
    titles = [annotation.text for annotation in figure.layout.annotations]
    assert titles[0].startswith("Per compound (r = ")
    assert titles[1].startswith("Per A-B disynthon (r = ")
    # Constant counts have no correlation.
    flat = compounds.head(3).with_columns(target_1=pl.lit(1), target_2=pl.lit(1))
    assert "r = n/a" in replicate_agreement(flat, flat, "x").layout.annotations[0].text
    assert len(ma_plot(compounds, top, "MAPK14").data) == 2


def test_sample_rows_is_seeded() -> None:
    frame = pl.DataFrame({"x": range(100)})
    assert sample_rows(frame, 200).equals(frame)
    first = sample_rows(frame, 10)
    assert first.height == 10
    assert first.equals(sample_rows(frame, 10))
    assert first.get_column("x").is_sorted()


def test_figures_json_is_stable(pipeline_dir: Path) -> None:
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    first = figures_to_json(build_figures(top, [], "MAPK14"))
    assert first == figures_to_json(build_figures(top, [], "MAPK14"))
    assert list(json.loads(first)) == ["cube"]


def test_optional_inputs(pipeline_dir: Path, tmp_path: Path) -> None:
    from conftest import run_cli  # noqa: PLC0415

    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    scores = pl.read_parquet(pipeline_dir / "mapk14.scores.parquet").filter(pl.col("target_sum") > 0)
    assert set(build_figures(top, [], "MAPK14", compounds=scores)) == {"cube", "ma"}
    output = tmp_path / "figures.json"
    args = [
        "visualise",
        "--top",
        pipeline_dir / "mapk14.top.parquet",
        "--validation",
        pipeline_dir / "bca.validation.json",
    ]
    assert run_cli(*args, "--screen", "MAPK14", "--output", output) == 0
    assert list(json.loads(output.read_text())) == ["cube"]


def test_saturation(pipeline_dir: Path) -> None:
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    counts_qc = json.loads((pipeline_dir / "mapk14.counts_qc.json").read_text())
    figure = saturation(counts_qc, "MAPK14")
    names = [trace.name for trace in figure.data]
    assert names[:3] == [
        f"Target (pooled): {counts_qc['library']['saturation']['target']['oversampling']:.2g} counts/compound",
        f"Matrix control (pooled): {counts_qc['library']['saturation']['matrix']['oversampling']:.2g} counts/compound",
        f"Pre-selection load: {counts_qc['library']['saturation']['load']['oversampling']:.2g} counts/compound",
    ]
    assert "Uniform library (Poisson)" in names
    assert figure.data[-1].x[-1] == "10+"
    assert tuple(figure.layout.xaxis.range) == (-3, pytest.approx(1.30103))  # the axis is not left to autorange
    assert "saturation" in build_figures(top, [], "MAPK14", counts_qc=counts_qc)
