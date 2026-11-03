"""Command-line interface used by the Nextflow processes.

Each subcommand is a thin wrapper: parse arguments, call the library, write files.
"""

import argparse
import json
import sys
from collections.abc import Callable, Sequence
from pathlib import Path

from . import __version__
from .analyse import analyse
from .fair import DatasetInfo, Provenance, describe_dataset, parse_sha256sums, read_dictionary
from .io import EXPORT_FORMATS, FeatureTable, export_table, parse_formats, read_table
from .report import discover_ids, load_bundle, parse_versions, render_report
from .visualise import build_figures, figures_to_json
from .wrangle import IMPUTE_STRATEGIES, normalise_name, wrangle


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _load_table(path: Path, id_column: str, target_column: str | None) -> FeatureTable:
    """Load a cleaned table (column names are already normalised)."""
    target = normalise_name(target_column) if target_column else None
    return FeatureTable(read_table(path), normalise_name(id_column), target)


def _keywords(text: str) -> tuple[str, ...]:
    return tuple(word.strip() for word in text.replace(";", ",").split(",") if word.strip())


def cmd_wrangle(args: argparse.Namespace) -> None:
    """Clean a raw table into Parquet plus a QC report."""
    result = wrangle(
        read_table(args.input),
        args.id_column,
        args.target_column or None,
        impute=args.impute,
        max_missing=args.max_missing,
    )
    result.table.frame.write_parquet(args.output)
    _write_json(args.qc, result.qc)


def cmd_analyse(args: argparse.Namespace) -> None:
    """Statistics, embedding, clustering and the supervised baseline."""
    table = _load_table(args.input, args.id_column, args.target_column)
    result = analyse(
        table,
        n_components=args.n_components,
        n_clusters=args.n_clusters,
        cv_folds=args.cv_folds,
        seed=args.seed,
        threads=args.threads,
    )
    result.embedding.write_parquet(args.embedding)
    _write_json(args.stats, result.stats)


def cmd_visualise(args: argparse.Namespace) -> None:
    """Build plotly figures as JSON."""
    table = _load_table(args.input, args.id_column, args.target_column)
    embedding = read_table(args.embedding)
    stats = json.loads(args.stats.read_text())
    args.output.write_text(figures_to_json(build_figures(table, embedding, stats)))


def cmd_export(args: argparse.Namespace) -> None:
    """Write the cleaned table in downstream-ready formats."""
    table = _load_table(args.input, args.id_column, args.target_column)
    args.outdir.mkdir(parents=True, exist_ok=True)
    export_table(table, args.prefix, parse_formats(args.formats), args.outdir)


def cmd_describe(args: argparse.Namespace) -> None:
    """Write FAIR metadata: schema.org JSON-LD and a Frictionless data package."""
    table = _load_table(args.input, args.id_column, args.target_column)
    stats = json.loads(args.stats.read_text()) if args.stats else {}
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
        feature_stats=stats.get("features"),
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

    sub = command("wrangle", cmd_wrangle)
    _table_options(sub)
    sub.add_argument("--impute", choices=IMPUTE_STRATEGIES, default="median")
    sub.add_argument("--max-missing", type=float, default=0.5, help="drop features missing more than this fraction")
    sub.add_argument("--output", type=Path, required=True, help="cleaned Parquet file")
    sub.add_argument("--qc", type=Path, required=True, help="QC summary JSON")

    sub = command("analyse", cmd_analyse)
    _table_options(sub)
    sub.add_argument("--n-components", type=int, default=2)
    sub.add_argument("--n-clusters", type=int, default=3)
    sub.add_argument("--cv-folds", type=int, default=5)
    sub.add_argument("--seed", type=int, default=42)
    sub.add_argument("--threads", type=int, default=1)
    sub.add_argument("--embedding", type=Path, required=True, help="embedding Parquet file")
    sub.add_argument("--stats", type=Path, required=True, help="statistics JSON")

    sub = command("visualise", cmd_visualise)
    _table_options(sub)
    sub.add_argument("--embedding", type=Path, required=True)
    sub.add_argument("--stats", type=Path, required=True)
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
    sub.add_argument("--stats", type=Path, help="statistics JSON (adds min/max per variable)")
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
