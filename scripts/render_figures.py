"""Render a run's plotly figures to PNG for the README (needs kaleido, which downloads Chrome once).

Usage: python scripts/render_figures.py RESULTS_DIR OUTPUT_DIR
"""

import json
import sys
from pathlib import Path

import kaleido
import plotly.io as pio


def main(results: Path, output: Path) -> None:
    """Write ``<screen>_<figure>.png`` for every ``screens/*/*.figures.json`` in ``results``."""
    kaleido.get_chrome_sync()
    output.mkdir(parents=True, exist_ok=True)
    for path in sorted(results.glob("screens/*/*.figures.json")):
        screen = path.name.removesuffix(".figures.json")
        for name, figure in json.loads(path.read_text()).items():
            target = output / f"{screen}_{name}.png"
            pio.write_image(pio.from_json(json.dumps(figure)), target, width=900, height=600, scale=2)
            print(f"wrote {target}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
