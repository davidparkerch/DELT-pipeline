import json
from pathlib import Path

import polars as pl

from delt_pipeline.visualise import build_figures, cube, figures_to_json, validation_bars


def test_figures(pipeline_dir: Path) -> None:
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    validation = json.loads((pipeline_dir / "mapk14.validation.json").read_text())
    figures = build_figures(top, validation, "MAPK14")
    assert set(figures) == {"cube", "validation"}
    assert len(figures["cube"].data[0].x) == top.height
    assert len(validation_bars(validation, "MAPK14").data) == 4  # bars + published markers, on and off
    assert set(build_figures(top, [], "BCA")) == {"cube"}
    assert cube(top, "MAPK14").layout.scene.xaxis.title.text == "Building block A"


def test_figures_json_is_stable(pipeline_dir: Path) -> None:
    top = pl.read_parquet(pipeline_dir / "mapk14.top.parquet")
    first = figures_to_json(build_figures(top, [], "MAPK14"))
    assert first == figures_to_json(build_figures(top, [], "MAPK14"))
    assert list(json.loads(first)) == ["cube"]
