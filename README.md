# DELT-pipeline

Process DNA-Encoded Library Technology data and exports FAIR, analysis-ready datasets..

A [Nextflow](https://www.nextflow.io/) 26.10 pipeline that uses static types, records,
workflow outputs, topic channels and data lineage, with FAIR metadata for every output.
It profiles tabular feature matrices (compound descriptors, expression or abundance matrices, ML feature sets):

```
samplesheet ─► PREPARE_TABLE (gzip, tr) ─► WRANGLE (polars) ─┬─► ANALYSE (scipy, scikit-learn) ─► VISUALISE (plotly)
                                                             └─► EXPORT (parquet/arrow/csv/tsv/jsonl/npz/svmlight) ─► CHECKSUM (sha256sum)
                                                     ─► DESCRIBE (JSON-LD, data package) ─► REPORT (one standalone HTML)
```

## Quick start

Requirements: Nextflow ≥ 26.10, Podman, [just](https://just.systems) and [uv](https://docs.astral.sh/uv/).

```bash
just run                                  # example data, processes in the runtime image (Podman)
nextflow run davidparkerch/delt-pipeline -profile podman --input samplesheet.csv --outdir results
```

## Input

A CSV, JSON or YAML samplesheet with one row per feature table (see `assets/example/samplesheet.csv`):

| Column | Required | Description |
|---|---|---|
| `id` | yes | Unique dataset id |
| `table` | yes | CSV, TSV, Parquet, Arrow or JSON Lines table (optionally `.gz`), or a `lid://` lineage URI |
| `id_column` | yes | Column with unique row identifiers (sample, compound, …) |
| `target_column` | no | Label or response column, which enables the supervised baseline |
| `description` | no | Free-text description (FAIR metadata) |
| `source` | no | DOI or URL of the original data |
| `license` | no | SPDX license id (default `--data_license`, CC-BY-4.0) |
| `keywords` | no | `;`-separated keywords |
| `dictionary` | no | Data dictionary CSV: `column,description,unit,unit_iri,term_iri` |

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | (required) | Samplesheet |
| `--outdir` | `results` | Output directory |
| `--export_formats` | `parquet,tsv,npz` | Any of `parquet, arrow, csv, tsv, jsonl, npz, svmlight` |
| `--impute` | `median` | `median`, `mean` or `drop` |
| `--max_missing` | `0.5` | Drop features with more missing values than this |
| `--n_components` / `--n_clusters` | `2` / `3` | PCA components and k-means clusters |
| `--cv_folds` | `5` | Cross-validation folds for the random-forest baseline |
| `--seed` | `42` | Random seed |
| `--report_title` | `DELT-pipeline report` | HTML report title |
| `--data_license` | `CC-BY-4.0` | Default SPDX license of the data |
| `--creator_name` / `--creator_orcid` | `David Parker` | Data creator recorded in the metadata |

## Outputs

```
results/
├── report.html                    # standalone, offline HTML report (plotly + jinja)
├── datasets.json                  # index of every dataset and its identifier
├── datasets/<id>/                 # clean.parquet, qc.json, stats.json, embedding.parquet, figures.json, metadata.jsonld
├── exports/<id>/                  # analysis-ready exports + SHA256SUMS + datapackage.json
├── ro-crate-metadata.json         # Workflow Run RO-Crate (run provenance)
└── pipeline_info/                 # versions.yml, bco.json (BioCompute Object), dag.html, trace, report, timeline
```

Export formats for downstream work: **parquet** or **arrow** (polars, pandas, Spark, DuckDB), **csv** or **tsv**
(R/Bioconductor, spreadsheets), **npz** (`X`, `y`, `ids`, `feature_names` for NumPy, PyTorch or JAX),
**jsonl** (LLM and data-loader pipelines) and **svmlight** (XGBoost, LightGBM, scikit-learn).

## FAIR

- **Findable**: content-derived identifiers (`urn:sha256:…`), `lid://` lineage ids, rich JSON-LD metadata and a searchable `datasets.json` index.
- **Accessible**: plain files over open protocols; metadata kept separately from the data.
- **Interoperable**: schema.org JSON-LD, Frictionless Data Packages, EDAM formats, ontology IRIs from the data dictionary.
- **Reusable**: SPDX licenses on every dataset, plus provenance from Nextflow lineage, Workflow Run RO-Crate, BioCompute Object and `versions.yml`.

Check a run with `just validate-fair`. The repository itself carries `CITATION.cff`, `codemeta.json` and a
Workflow RO-Crate (ready for [WorkflowHub](https://workflowhub.eu)). Enable the
[Zenodo GitHub integration](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)
to mint a DOI for every release.

## Development

Open the folder in the Podman dev container (`.devcontainer/`; set the Dev Containers extension's container CLI to
`podman`), or use `just dev-shell`. Then:

```bash
just            # list commands
just fmt        # format Python and Nextflow
just ci         # lint and all tests (pytest 100% coverage, nf-test snapshots)
```

See [AGENTS.md](AGENTS.md) for conventions.

## Citation

See [CITATION.cff](CITATION.cff).

## License

Apache-2.0. See [LICENSE](LICENSE). Example data in `assets/example/` is synthetic and released under CC0-1.0.
