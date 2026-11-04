"""Render the explorer page: one self-contained HTML file (plotly.js, styles, script and data inline)."""

import json
from importlib import resources
from typing import Any

from jinja2 import Environment, PackageLoader, StrictUndefined, select_autoescape
from plotly.offline import get_plotlyjs

RDKIT_JS_VERSION = "2026.9.1"
"""RDKit.js is loaded from a CDN (its WebAssembly build is 7 MB); the page works without it, minus structures."""


def render_explorer(payload: dict[str, Any]) -> str:
    """Render the page.

    Args:
        payload: From :func:`~delt_pipeline.explore.payload.build_payload`; ``payload["meta"]`` must hold
            ``title``, ``license`` and ``license_url``.

    Returns:
        The HTML document. The output depends only on the payload, so it is byte-stable.
    """
    environment = Environment(
        loader=PackageLoader(__package__, "templates"),
        autoescape=select_autoescape(default=True),
        undefined=StrictUndefined,
        keep_trailing_newline=True,
    )
    # No "<" may survive in the script element ("</script>", "<!--"); \u003c is the same character in JSON.
    data = json.dumps(payload, sort_keys=True, separators=(",", ":")).replace("<", "\\u003c")
    return environment.get_template("explorer.html.j2").render(
        meta=payload["meta"],
        data=data,
        plotly_js=get_plotlyjs(),
        app_js=resources.files(__package__).joinpath("static/app.js").read_text(encoding="utf-8"),
        rdkit_version=RDKIT_JS_VERSION,
    )
