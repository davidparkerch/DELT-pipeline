"""Standalone HTML report rendered with jinja2."""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from jinja2 import Environment, PackageLoader, StrictUndefined, select_autoescape
from plotly.offline import get_plotlyjs

from ..analyse.score import SCORE_LABELS

METADATA_SUFFIX = ".metadata.jsonld"
PARTS = {
    "counts_qc": ".counts_qc.json",
    "score_qc": ".score_qc.json",
    "synthon_qc": ".synthon_qc.json",
    "validation": ".validation.json",
    "model": ".model.json",
    "figures": ".figures.json",
    "metadata": METADATA_SUFFIX,
}


@dataclass(frozen=True)
class DatasetBundle:
    """Everything the report shows for one screen."""

    id: str
    counts_qc: dict[str, Any]
    score_qc: dict[str, Any]
    synthon_qc: dict[str, Any]
    validation: list[dict[str, Any]]
    model: dict[str, Any]
    figures: dict[str, Any]
    metadata: dict[str, Any]


VARIANT_LABELS = {"counts": "Counts only", "control_weighted": "Weighted by the control"}
SPLIT_LABELS = {"random": "Random", "disynthon": "Disynthon"}


def model_rows(model: dict[str, Any]) -> list[dict[str, Any]]:
    """The fingerprint-model results as table rows: one per variant and split, with KinDEL's ML range per split.

    Each row holds the label, ``on`` and ``off`` (-Spearman vs KD) and, for the split, the lowest and highest
    value among KinDEL's published machine-learning models (``None`` when the target has none).
    """
    published = model.get("published_ml", [])
    rows = []
    order = {key: index for index, key in enumerate([*VARIANT_LABELS, *SPLIT_LABELS])}
    for variant, split in sorted(
        {(r["variant"], r["split"]) for r in model.get("results", [])},
        key=lambda pair: (order[pair[0]], order[pair[1]]),
    ):
        by_condition = {
            r["condition"]: r["spearman"] for r in model["results"] if (r["variant"], r["split"]) == (variant, split)
        }
        same_split = [row for row in published if row["split"] == split]
        rows.append(
            {
                "label": f"{VARIANT_LABELS[variant]}, {SPLIT_LABELS[split].lower()} split",
                "on": by_condition.get("on"),
                "off": by_condition.get("off"),
                "published_on": (min(r["on_dna"] for r in same_split), max(r["on_dna"] for r in same_split))
                if same_split
                else None,
                "published_off": (min(r["off_dna"] for r in same_split), max(r["off_dna"] for r in same_split))
                if same_split
                else None,
            }
        )
    return rows


def discover_ids(directory: Path) -> list[str]:
    """Dataset ids in ``directory``, found from their ``<id>.metadata.jsonld`` files."""
    return sorted(path.name.removesuffix(METADATA_SUFFIX) for path in directory.glob(f"*{METADATA_SUFFIX}"))


def load_bundle(directory: Path, dataset_id: str) -> DatasetBundle:
    """Load the QC, score, validation, figure and metadata files of one screen (``<id>.<part>``)."""
    parts = {key: json.loads((directory / f"{dataset_id}{suffix}").read_text()) for key, suffix in PARTS.items()}
    return DatasetBundle(id=dataset_id, **parts)


def parse_versions(path: Path) -> list[tuple[str, str]]:
    """Parse ``tool: version`` lines (the pipeline's ``versions.yml``), sorted by tool."""
    versions = []
    for line in path.read_text().splitlines():
        tool, separator, version = line.partition(":")
        if separator and tool.strip() and not line.startswith((" ", "#")):
            versions.append((tool.strip(), version.strip().strip("'\"")))
    return sorted(versions)


def render_report(
    title: str,
    bundles: list[DatasetBundle],
    versions: list[tuple[str, str]],
    pipeline: dict[str, str],
    triage: dict[str, Any] | None = None,
) -> str:
    """Render the self-contained report (plotly.js, CSS and data are all inline).

    Args:
        title: Report title.
        bundles: One bundle per dataset.
        versions: ``(tool, version)`` pairs.
        pipeline: ``name``, ``version`` and optional ``homepage`` of the generating pipeline.
        triage: The triage summary (``analyse.triage.triage_qc``), if the run triaged hits.

    Returns:
        The HTML document.
    """
    environment = Environment(
        loader=PackageLoader(__package__, "templates"),
        autoescape=select_autoescape(default=True),
        undefined=StrictUndefined,
        keep_trailing_newline=True,
    )
    figures = {
        f"fig-{bundle.id}-{name}": figure for bundle in bundles for name, figure in sorted(bundle.figures.items())
    }
    template = environment.get_template("report.html.j2")
    return template.render(
        title=title,
        model_rows=model_rows,
        bundles=bundles,
        versions=versions,
        pipeline=pipeline,
        figures=figures,
        triage=triage,
        score_labels=SCORE_LABELS,
        plotly_js=get_plotlyjs(),
    )
