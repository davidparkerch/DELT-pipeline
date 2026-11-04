import base64
import json
import re
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from rdkit.Chem import MolFromSmiles, rdFingerprintGenerator

from delt_pipeline.explore import BLOCK_COLUMNS, HIT_COLUMNS, build_payload, fingerprint_bytes, render_explorer
from delt_pipeline.explore.payload import blockset_map, hit_rows

META = {"title": "Test", "license": "CC0", "license_url": "https://example.org/license", "report_url": None}


def _inputs(pipeline_dir: Path) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, dict[str, pl.DataFrame]]:
    return (
        pl.read_parquet(pipeline_dir / "hits.parquet"),
        pl.read_parquet(pipeline_dir / "picks.parquet"),
        pl.read_parquet(pipeline_dir / "blocks.parquet"),
        {s: pl.read_parquet(pipeline_dir / f"{s}.blockset.parquet") for s in ("mapk14", "bca")},
    )


def test_fingerprint_bytes_are_the_morgan_bits() -> None:
    smiles = ["c1ccccc1CC(=O)N", "CCO"]
    packed = fingerprint_bytes(smiles)
    assert [len(item) for item in packed] == [256, 256]
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for text, item in zip(smiles, packed, strict=True):
        bits = np.unpackbits(np.frombuffer(item, dtype=np.uint8))
        assert np.flatnonzero(bits).tolist() == sorted(generator.GetFingerprint(MolFromSmiles(text)).GetOnBits())
    with pytest.raises(ValueError, match="cannot parse"):
        fingerprint_bytes(["not a molecule"])


def test_hit_rows(pipeline_dir: Path) -> None:
    hits, picks, _blocks, _sets = _inputs(pipeline_dir)
    rows = hit_rows(hits, picks)
    assert len(rows) == hits.height
    column = {name: index for index, name in enumerate(HIT_COLUMNS)}
    assert all(len(row) == len(HIT_COLUMNS) for row in rows)
    assert [row[column["rank"]] for row in rows] == sorted(row[column["rank"]] for row in rows)
    picked = {row[column["hit_id"]]: row[column["pick_rank"]] for row in rows if row[column["pick_rank"]]}
    assert picked == dict(picks.select("hit_id", "pick_rank").iter_rows())
    assert {row[column["replicate_support"]] for row in rows} <= {0, 1}  # booleans become integers
    assert all(not isinstance(value, float) or np.isfinite(value) for row in rows for value in row)  # no inf or NaN


def test_payload_and_page(pipeline_dir: Path) -> None:
    hits, picks, blocks, sets = _inputs(pipeline_dir)
    payload = build_payload(hits, picks, blocks, sets, META)
    assert payload["screens"] == ["mapk14"]
    assert len(payload["hits"]) == len(payload["fingerprints"]) == hits.height
    assert payload["hit_columns"] == list(HIT_COLUMNS)
    assert payload["block_columns"] == list(BLOCK_COLUMNS)
    assert len(payload["blocks"]) == blocks.height
    assert set(payload["blocksets"]) == {"bca", "mapk14"}
    assert len(base64.b64decode(payload["fingerprints"][0])) == 256
    key = "B0005"
    assert payload["blocksets"]["mapk14"][key][0] == pytest.approx(
        sets["mapk14"].filter((pl.col("position") == "b") & (pl.col("bb_id") == 5)).item(0, "auc"), rel=1e-3
    )
    assert blockset_map({}) == {}

    page = render_explorer(payload)
    assert render_explorer(payload) == page  # byte-stable
    match = re.search(r'<script id="payload" type="application/json">(.*?)</script>', page, re.DOTALL)
    assert match
    assert json.loads(match.group(1)) == payload
    assert "plotly.js v" in page
    assert "<title>Test</title>" in page
    assert "https://example.org/license" in page


def test_page_cannot_be_broken_out_of(pipeline_dir: Path) -> None:
    hits, picks, blocks, sets = _inputs(pipeline_dir)
    hostile = hits.with_columns(reasons=pl.lit("</script><script>alert(1)</script>"))
    page = render_explorer(build_payload(hostile, picks, blocks, sets, {**META, "title": "<b>x</b>"}))
    assert "<script>alert(1)" not in page
    assert "<title>&lt;b&gt;x&lt;/b&gt;</title>" in page
    match = re.search(r'<script id="payload" type="application/json">(.*?)</script>', page, re.DOTALL)
    assert match
    assert json.loads(match.group(1))["hits"][0][HIT_COLUMNS.index("reasons")] == "</script><script>alert(1)</script>"


def test_explore_cli_output(pipeline_dir: Path) -> None:
    page = (pipeline_dir / "explorer.html").read_text(encoding="utf-8")
    assert "<title>Test explorer</title>" in page
    assert 'href="index.html"' in page
    assert "https://example.org/demo" in page
