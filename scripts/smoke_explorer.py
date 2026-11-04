"""Drive the DELT Explorer in a headless browser and check it against values computed independently in Python.

The page's JavaScript is thin glue, so it has no unit tests of its own. This script is its test: it opens a
built page, filters, searches, marks and exports through the real interface and compares what the page shows
with what RDKit and numpy compute from the same payload. It needs a browser and (for structures and search)
network access to the RDKit.js CDN, so it is not part of ``just ci``.

Run: ``just smoke-explorer results/explorer.html`` (installs playwright on the fly with ``uv run --with``).
``--frames DIR`` saves screenshots (the README animation is made from them); ``--chromium PATH`` uses a
browser you already have.
"""

import argparse
import base64
import csv
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from playwright.sync_api import Page, sync_playwright
from rdkit import Chem

PAYLOAD = re.compile(r'<script id="payload" type="application/json">(.*?)</script>', re.DOTALL)
SULFONAMIDE = "S(=O)(=O)N"
MAX_PAGE_BYTES = 12_000_000
MIN_LOWER_BOUND = 5.0
MIN_SIMILARITY = 0.5
ROWS_PER_PAGE = 100
FINGERPRINTS_CHECKED = 5
WEBGL_FLAGS = ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"]


@dataclass
class Context:
    """The page under test, the payload it embeds and where screenshots go."""

    page: Page
    payload: dict[str, Any]
    frames: Path | None
    screen: str

    def col(self, name: str) -> int:
        """Index of a hit column."""
        return self.payload["hit_columns"].index(name)

    @property
    def in_screen(self) -> list[list[Any]]:
        """The hit rows of the screen shown."""
        return [row for row in self.payload["hits"] if row[self.col("screen")] == self.screen]

    def ids(self, rows: list[list[Any]]) -> set[str]:
        """The hit ids of some rows."""
        return {row[self.col("hit_id")] for row in rows}

    def visible(self) -> set[str]:
        """The hit ids the page currently shows."""
        return set(self.page.evaluate("window.DELT.visible()"))

    def marked(self) -> int:
        """How many hits are marked in the page."""
        return self.page.evaluate("window.DELT.state.marked.size")

    def shot(self, name: str) -> None:
        """Save a screenshot when ``--frames`` was given."""
        if self.frames:
            self.page.screenshot(path=str(self.frames / f"{name}.png"))


def check(condition: bool, message: str) -> None:
    """Print the result and stop at the first failure."""
    print(("ok    " if condition else "FAIL  ") + message)
    if not condition:
        sys.exit(1)


def unpack(text: str) -> np.ndarray:
    """A base64 fingerprint as a 0/1 vector."""
    return np.unpackbits(np.frombuffer(base64.b64decode(text), dtype=np.uint8))


def tanimoto(a: np.ndarray, b: np.ndarray) -> float:
    """Tanimoto similarity of two bit vectors."""
    common = int(np.sum(a & b))
    total = int(np.sum(a) + np.sum(b)) - common
    return common / total if total else 0.0


def set_tiers(ctx: Context, tiers: set[str]) -> None:
    """Tick exactly these tier boxes."""
    script = (
        "(tiers) => document.querySelectorAll('#tiers input').forEach(b => { b.checked = tiers.includes(b.value); })"
    )
    ctx.page.evaluate(script, sorted(tiers))
    ctx.page.locator("#tiers").dispatch_event("change")


def check_start(ctx: Context) -> None:
    """The page starts with every hit of the first screen, both plots and a working RDKit.js."""
    check(ctx.page.evaluate("window.DELT.rdkitState") == "ready", "RDKit.js loaded")
    check(ctx.visible() == ctx.ids(ctx.in_screen), f"{ctx.screen}: all {len(ctx.in_screen)} hits shown")
    check(ctx.page.locator("#scatter canvas").count() > 0, "scatter plot drawn")
    check(ctx.page.locator("#cube canvas").count() > 0, "cube drawn")
    ctx.shot("01-overview")
    pairs = zip(ctx.payload["hits"], ctx.payload["fingerprints"], strict=True)
    for row, packed in list(pairs)[:FINGERPRINTS_CHECKED]:
        browser_bits = ctx.page.evaluate("(smiles) => window.DELT.fingerprintOf(smiles)", row[ctx.col("smiles")])
        check(
            bytes(browser_bits) == base64.b64decode(packed),
            f"RDKit.js fingerprint equals Python for rank {row[ctx.col('rank')]}",
        )


def check_filters(ctx: Context) -> None:
    """Tier, shortlist and score filters show what Python computes."""
    rows = ctx.in_screen
    set_tiers(ctx, {"1"})
    tier_one = [row for row in rows if row[ctx.col("tier")] == "1"]
    check(ctx.visible() == ctx.ids(tier_one), f"tier filter: {len(tier_one)} tier-1 hits")
    ctx.page.locator("#picks-only").check()
    picks = [row for row in tier_one if row[ctx.col("pick_rank")]]
    check(bool(picks) and ctx.visible() == ctx.ids(picks), f"shortlist filter: {len(picks)} picks")
    ctx.page.locator("#picks-only").uncheck()
    set_tiers(ctx, {"1", "2", "3", "flagged"})
    ctx.page.locator("#min-lb").fill(str(MIN_LOWER_BOUND))
    strong = [row for row in rows if row[ctx.col("poisson_lower")] >= MIN_LOWER_BOUND]
    check(ctx.visible() == ctx.ids(strong), f"lower-bound filter: {len(strong)} hits with LB >= {MIN_LOWER_BOUND}")
    ctx.page.locator("#min-lb").fill("0")


def check_substructure(ctx: Context) -> None:
    """Substructure search finds the blocks RDKit finds in Python, and filters the hits to them."""
    block_col = {name: i for i, name in enumerate(ctx.payload["block_columns"])}
    query = Chem.MolFromSmarts(SULFONAMIDE)
    matched = {
        row[block_col["block_id"]]
        for row in ctx.payload["blocks"]
        if Chem.MolFromSmiles(row[block_col["fragment_smiles"]]).HasSubstructMatch(query)
    }
    ctx.page.locator("#query").fill(SULFONAMIDE)
    ctx.page.locator("#search").click()
    shown = set(ctx.page.locator("#blocks .card b").all_inner_texts())
    check(
        ctx.page.locator("#blocks-panel").is_visible() and bool(matched),
        f"substructure search: {len(matched)} blocks match {SULFONAMIDE}",
    )
    check(bool(shown) and shown <= matched, "the blocks shown are among the ones RDKit matches in Python")
    containing = [
        row
        for row in ctx.in_screen
        if matched & {f"A{row[ctx.col('bb_a')]:04d}", f"B{row[ctx.col('bb_b')]:04d}", f"C{row[ctx.col('bb_c')]:04d}"}
    ]
    check(ctx.visible() == ctx.ids(containing), f"hits filtered to those with a matching block: {len(containing)}")
    ctx.shot("02-substructure")
    ctx.page.locator("#clear-search").click()


def check_similarity(ctx: Context) -> None:
    """Similarity search keeps the hits whose Tanimoto similarity (computed here from the same bits) is high enough."""
    hits = ctx.payload["hits"]
    bits = [unpack(text) for text in ctx.payload["fingerprints"]]
    reference = ctx.in_screen[0]
    index = hits.index(reference)
    ctx.page.locator("#mode").select_option("sim")
    ctx.page.locator("#min-sim").fill(str(MIN_SIMILARITY))
    ctx.page.locator("#query").fill(reference[ctx.col("smiles")])
    ctx.page.locator("#search").click()
    similar = [
        row
        for i, row in enumerate(hits)
        if row[ctx.col("screen")] == ctx.screen and tanimoto(bits[index], bits[i]) >= MIN_SIMILARITY
    ]
    check(
        ctx.visible() == ctx.ids(similar), f"similarity search: {len(similar)} hits with Tanimoto >= {MIN_SIMILARITY}"
    )
    check(reference[ctx.col("hit_id")] in ctx.visible(), "the query compound finds itself")
    ctx.page.locator("#clear-search").click()
    ctx.page.locator("#mode").select_option("sub")
    ctx.page.locator("#query").fill("[[")
    ctx.page.locator("#search").click()
    check("Cannot search" in ctx.page.locator("#status").inner_text(), "an invalid query is reported, not thrown")
    ctx.page.locator("#clear-search").click()


def lasso(ctx: Context) -> None:
    """Draw a lasso over the upper right of the scatter plot."""
    box = ctx.page.locator("#scatter .nsewdrag").bounding_box()
    x0, y0, width, height = box["x"], box["y"], box["width"], box["height"]
    ctx.page.mouse.move(x0 + width * 0.55, y0 + height * 0.1)
    ctx.page.mouse.down()
    for fx, fy in ((0.95, 0.1), (0.95, 0.9), (0.55, 0.9), (0.55, 0.1)):
        ctx.page.mouse.move(x0 + width * fx, y0 + height * fy, steps=6)
    ctx.page.mouse.up()
    ctx.page.wait_for_timeout(500)


def check_marking(ctx: Context) -> None:
    """Marks from the scatter, the table and the cube agree everywhere and export to CSV."""
    lasso(ctx)
    marked = set(ctx.page.evaluate("[...window.DELT.state.marked]"))
    check(bool(marked), f"lasso in the scatter marks {len(marked)} hits")
    first_page = ctx.ids(sorted(ctx.in_screen, key=lambda row: row[ctx.col("rank")])[:ROWS_PER_PAGE])
    check(
        ctx.page.locator("#hits tr.marked").count() == len(first_page & marked),
        "the same hits are marked in the table (first page)",
    )
    check(f"{len(marked)} marked" in ctx.page.locator("#status").inner_text(), "the status line counts the marked hits")
    ctx.shot("03-marked")
    ctx.page.locator("#hits tbody tr:not(.marked) input.mark").first.click()
    check(ctx.marked() == len(marked) + 1, "a table checkbox adds a mark")
    ctx.page.evaluate(
        "document.getElementById('cube').emit('plotly_click', {points: [{curveNumber: 0, pointNumber: 0}]})"
    )
    check(ctx.marked() in {len(marked), len(marked) + 2}, "a click in the cube toggles a mark")
    with ctx.page.expect_download() as download:
        ctx.page.locator("#export").click()
    rows = list(csv.reader(Path(download.value.path()).read_text(encoding="utf-8").splitlines()))
    check(rows[0] == ctx.payload["hit_columns"], "exported CSV has the hit columns")
    check(len(rows) - 1 == ctx.marked(), f"exported CSV has the {len(rows) - 1} marked hits")
    ctx.page.locator("#clear-marks").click()
    check(ctx.marked() == 0, "marks can be cleared")


def check_detail_and_screens(ctx: Context) -> None:
    """A row opens the compound with its structure, and the screen can be switched."""
    ctx.page.locator("#hits tbody tr.hit").first.click()
    check(
        ctx.page.locator("#detail .mol svg").count() >= 1, "clicking a row draws the compound and its building blocks"
    )
    ctx.shot("04-detail")
    for other in (s for s in ctx.payload["screens"] if s != ctx.screen):
        ctx.page.locator("#screen").select_option(other)
        check(ctx.page.evaluate("window.DELT.state.screen") == other and bool(ctx.visible()), f"switching to {other}")
        break


def main() -> None:
    """Open the page and run every check."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("page", type=Path)
    parser.add_argument("--frames", type=Path, help="save screenshots here")
    parser.add_argument("--height", type=int, default=1100, help="browser window height in pixels")
    parser.add_argument("--chromium", type=Path, help="Chromium executable (default: the one playwright installed)")
    args = parser.parse_args()
    if args.frames:
        args.frames.mkdir(parents=True, exist_ok=True)
    payload = json.loads(PAYLOAD.search(args.page.read_text(encoding="utf-8")).group(1))
    check(args.page.stat().st_size < MAX_PAGE_BYTES, f"page is {args.page.stat().st_size / 1e6:.1f} MB (< 12 MB)")

    errors: list[str] = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=str(args.chromium) if args.chromium else None,
            args=WEBGL_FLAGS,  # headless browsers have no GPU: draw WebGL in software
        )
        page = browser.new_page(viewport={"width": 1500, "height": args.height}, accept_downloads=True)
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(args.page.resolve().as_uri())
        page.wait_for_selector("#hits tbody tr")
        page.evaluate("window.DELT.ready")
        ctx = Context(page, payload, args.frames, page.evaluate("window.DELT.state.screen"))
        for step in (
            check_start,
            check_filters,
            check_substructure,
            check_similarity,
            check_marking,
            check_detail_and_screens,
        ):
            step(ctx)
        browser.close()
    check(not errors, f"no console errors{': ' + '; '.join(errors) if errors else ''}")
    print("explorer smoke test passed")


if __name__ == "__main__":
    main()
