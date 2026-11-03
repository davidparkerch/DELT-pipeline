"""Command-line interface used by the Nextflow processes.

Each subcommand is a thin wrapper: parse arguments, call the library, write files.
"""

import argparse
import json
import sys
from collections.abc import Callable, Sequence
from pathlib import Path

import polars as pl

from . import __version__
from .analyse import SCORES, absent_scores, read_heldout, score_table, validate_scores
from .fair import DatasetInfo, Provenance, describe_dataset, parse_sha256sums, read_dictionary
from .io import EXPORT_FORMATS, FeatureTable, export_table, normalise_name, parse_formats, read_table
from .library import read_building_blocks, select_sublibrary, simulate_library, write_counts
from .library.schema import ID
from .report import discover_ids, load_bundle, parse_versions, render_report
from .visualise import build_figures, figures_to_json


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _read_json(path: Path) -> object:
    return json.loads(path.read_text())


def _load_table(path: Path, id_column: str, target_column: str | None) -> FeatureTable:
    """Load a table whose column names are already export-safe."""
    target = normalise_name(target_column) if target_column else None
    return FeatureTable(read_table(path), normalise_name(id_column), target)


def _keywords(text: str) -> tuple[str, ...]:
    return tuple(word.strip() for word in text.replace(";", ",").split(",") if word.strip())


def cmd_simulate(args: argparse.Namespace) -> None:
    """Write a synthetic library in the KinDEL layout (screens, held-out KD and synthon tables)."""
    library = simulate_library(args.seed)
    for name in ("heldout", "synthons"):
        (args.outdir / name).mkdir(parents=True, exist_ok=True)
    for screen, frame in library.screens.items():
        frame.write_parquet(args.outdir / f"{screen}.parquet")
    for name, frame in library.heldout.items():
        frame.write_csv(args.outdir / "heldout" / f"{name}.csv")
    for position, frame in library.building_blocks.items():
        frame.write_csv(args.outdir / "synthons" / f"reacted_{position}.csv")


def cmd_select_sublibrary(args: argparse.Namespace) -> None:
    """Choose the BB-A building blocks of the sub-library (all held-out ones plus random ones)."""
    heldout = [pl.read_csv(path, infer_schema_length=None) for path in args.heldout]
    selection = select_sublibrary(read_building_blocks(args.synthons), heldout, args.n_bb_a, args.seed)
    _write_json(args.output, selection)


def cmd_subset(args: argparse.Namespace) -> None:
    """Stream a full KinDEL screen into the compact sub-library count table."""
    qc = write_counts(args.input, read_building_blocks(args.synthons), _read_json(args.sublibrary), args.output)
    _write_json(args.qc, qc)


def cmd_score(args: argparse.Namespace) -> None:
    """Score every compound; also write the top compounds as the published dataset."""
    if args.top_n < 1:
        msg = f"--top-n must be at least 1, got {args.top_n}"
        raise ValueError(msg)
    scores, qc = score_table(
        pl.read_parquet(args.input),
        _read_json(args.counts_qc),
        _read_json(args.sublibrary),
        alpha=args.alpha,
        level=args.level,
    )
    scores.write_parquet(args.output, statistics=True)
    scores.head(args.top_n).write_parquet(args.top)
    _write_json(args.qc, qc)


def cmd_validate(args: argparse.Namespace) -> None:
    """Correlate every score with measured KD (on- and off-DNA)."""
    assays = [("on", args.ondna), ("off", args.offdna)]
    results = []
    if any(path for _, path in assays):
        scores = pl.read_parquet(args.scores, columns=[ID, *SCORES])
        absent = absent_scores(_read_json(args.score_qc))
        for condition, path in assays:
            if path:
                heldout = read_heldout(pl.read_csv(path, infer_schema_length=None))
                results.append(validate_scores(scores, heldout, absent, args.target, condition))
    _write_json(args.output, results)


def cmd_visualise(args: argparse.Namespace) -> None:
    """Build plotly figures as JSON."""
    figures = build_figures(pl.read_parquet(args.top), _read_json(args.validation), args.screen)
    args.output.write_text(figures_to_json(figures))


def cmd_export(args: argparse.Namespace) -> None:
    """Write the cleaned table in downstream-ready formats."""
    table = _load_table(args.input, args.id_column, args.target_column)
    args.outdir.mkdir(parents=True, exist_ok=True)
    export_table(table, args.prefix, parse_formats(args.formats), args.outdir)


def cmd_describe(args: argparse.Namespace) -> None:
    """Write FAIR metadata: schema.org JSON-LD and a Frictionless data package."""
    table = _load_table(args.input, args.id_column, args.target_column)
    description = describe_dataset(
        table,
        args.input,
        list(args.exports),
        DatasetInfo(
            id=args.id,
            description=args.description,
            source=args.source,
            license=args.license,
            keywords=_keywords(args.keywords),
            creator_name=args.creator_name,
            creator_orcid=args.creator_orcid,
        ),
        Provenance(args.pipeline, args.pipeline_version, args.homepage),
        checksums=parse_sha256sums(args.sha256sums) if args.sha256sums else None,
        dictionary=read_dictionary(args.dictionary),
        exports_base=args.exports_base,
    )
    _write_json(args.jsonld, description.jsonld)
    _write_json(args.datapackage, description.datapackage)
    if args.identifier:
        args.identifier.write_text(description.identifier + "\n")


def cmd_report(args: argparse.Namespace) -> None:
    """Render the standalone HTML report from every dataset in a directory."""
    ids = discover_ids(args.input_dir)
    if not ids:
        msg = f"no '*.metadata.jsonld' files found in {args.input_dir}"
        raise ValueError(msg)
    bundles = [load_bundle(args.input_dir, dataset_id) for dataset_id in ids]
    pipeline = {"name": args.pipeline, "version": args.pipeline_version, "homepage": args.homepage}
    html = render_report(args.title, bundles, parse_versions(args.versions), pipeline)
    args.output.write_text(html)


def _table_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", type=Path, required=True, help="input table")
    parser.add_argument("--id-column", required=True, help="column with unique row identifiers")
    parser.add_argument("--target-column", default="", help="optional label/response column")


def _provenance_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--pipeline", required=True, help="generating pipeline name")
    parser.add_argument("--pipeline-version", required=True, help="generating pipeline version")
    parser.add_argument("--homepage", default="", help="pipeline homepage URL")


def build_parser() -> argparse.ArgumentParser:  # noqa: PLR0915 - one block per subcommand
    """Build the argument parser with one subcommand per pipeline step."""
    parser = argparse.ArgumentParser(prog=__package__, description=__doc__.splitlines()[0])
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)

    def command(name: str, handler: Callable[[argparse.Namespace], None]) -> argparse.ArgumentParser:
        sub = commands.add_parser(name, help=handler.__doc__, description=handler.__doc__)
        sub.set_defaults(handler=handler)
        return sub

    sub = command("simulate", cmd_simulate)
    sub.add_argument("--seed", type=int, default=7)
    sub.add_argument("--outdir", type=Path, required=True)

    sub = command("select-sublibrary", cmd_select_sublibrary)
    sub.add_argument("--synthons", type=Path, required=True, help="folder with reacted_{a,b,c}.csv")
    sub.add_argument("--heldout", type=Path, nargs="*", default=[], help="held-out KD tables")
    sub.add_argument("--n-bb-a", type=int, default=60, help="number of BB-A building blocks to keep")
    sub.add_argument("--seed", type=int, default=42)
    sub.add_argument("--output", type=Path, required=True, help="sub-library JSON")

    sub = command("subset", cmd_subset)
    sub.add_argument("--input", type=Path, required=True, help="full KinDEL screen (Parquet)")
    sub.add_argument("--synthons", type=Path, required=True, help="folder with reacted_{a,b,c}.csv")
    sub.add_argument("--sublibrary", type=Path, required=True, help="sub-library JSON")
    sub.add_argument("--output", type=Path, required=True, help="count table (Parquet)")
    sub.add_argument("--qc", type=Path, required=True, help="QC JSON (includes whole-library depth)")

    sub = command("score", cmd_score)
    sub.add_argument("--input", type=Path, required=True, help="count table (Parquet)")
    sub.add_argument("--counts-qc", type=Path, required=True)
    sub.add_argument("--sublibrary", type=Path, required=True)
    sub.add_argument("--alpha", type=float, default=0.5, help="pseudocount of the enrichment ratio")
    sub.add_argument("--level", type=float, default=0.95, help="confidence level of the rate-ratio interval")
    sub.add_argument("--top-n", type=int, default=2000, help="compounds in the top table")
    sub.add_argument("--output", type=Path, required=True, help="scores (Parquet)")
    sub.add_argument("--top", type=Path, required=True, help="top compounds (Parquet)")
    sub.add_argument("--qc", type=Path, required=True, help="score QC JSON")

    sub = command("validate", cmd_validate)
    sub.add_argument("--scores", type=Path, required=True)
    sub.add_argument("--score-qc", type=Path, required=True)
    sub.add_argument("--target", required=True, help="target name, e.g. mapk14")
    sub.add_argument("--ondna", type=Path, help="held-out on-DNA KD table")
    sub.add_argument("--offdna", type=Path, help="held-out off-DNA KD table")
    sub.add_argument("--output", type=Path, required=True, help="validation JSON")

    sub = command("visualise", cmd_visualise)
    sub.add_argument("--top", type=Path, required=True, help="top compounds (Parquet)")
    sub.add_argument("--validation", type=Path, required=True)
    sub.add_argument("--screen", required=True)
    sub.add_argument("--output", type=Path, required=True, help="figures JSON")

    sub = command("export", cmd_export)
    _table_options(sub)
    sub.add_argument("--formats", default="parquet", help=f"comma-separated: {', '.join(EXPORT_FORMATS)}")
    sub.add_argument("--prefix", required=True, help="output file name stem")
    sub.add_argument("--outdir", type=Path, default=Path())

    sub = command("describe", cmd_describe)
    _table_options(sub)
    _provenance_options(sub)
    sub.add_argument("--id", required=True, help="dataset id")
    sub.add_argument("--exports", type=Path, nargs="*", default=[], help="exported files to describe")
    sub.add_argument("--exports-base", default="", help="URL prefix from the JSON-LD file to the exports")
    sub.add_argument("--sha256sums", type=Path, help="sha256sum output covering the exports")
    sub.add_argument("--dictionary", type=Path, help="data dictionary CSV")
    sub.add_argument("--description", default="")
    sub.add_argument("--source", default="", help="DOI or URL of the original data")
    sub.add_argument("--license", default="CC-BY-4.0", help="SPDX license id")
    sub.add_argument("--keywords", default="", help="comma-separated keywords")
    sub.add_argument("--creator-name", default="")
    sub.add_argument("--creator-orcid", default="")
    sub.add_argument("--jsonld", type=Path, required=True)
    sub.add_argument("--datapackage", type=Path, required=True)
    sub.add_argument("--identifier", type=Path, help="also write the dataset identifier to this file")

    sub = command("report", cmd_report)
    _provenance_options(sub)
    sub.add_argument("--input-dir", type=Path, default=Path())
    sub.add_argument("--versions", type=Path, required=True, help="versions.yml")
    sub.add_argument("--title", default="Pipeline report")
    sub.add_argument("--output", type=Path, required=True, help="HTML report")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI; returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.handler(args)
    except (ValueError, FileNotFoundError) as error:
        print(f"{parser.prog} {args.command}: error: {error}", file=sys.stderr)
        return 1
    return 0
