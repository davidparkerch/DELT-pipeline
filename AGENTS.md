# AGENTS.md: DELT-pipeline

Guidance for AI coding agents (and humans) working in this repository.

## What this is

A Nextflow 26.10 pipeline that scores DNA-encoded library (DEL) selections from
[KinDEL](https://github.com/insitro/kindel), validates the scores against measured KD and publishes FAIR
hit tables. Nextflow **only orchestrates**. The real work is done by:

- standard Linux tools (`sha256sum`, …), and
- the Python package `src/delt_pipeline/` (CLI: `delt_pipeline`), with these subpackages:
  `library/` (simulation, building blocks and their RDKit descriptors, streamed sub-library, depth and saturation),
  `analyse/` (enrichment scores, validation, building-block set enrichment, resynthesis shortlist, fingerprint model;
  scipy, scikit-learn, RDKit), `visualise/` (plotly), `report/` (jinja2 HTML), `explore/` (the static DELT Explorer page)
  and `fair/` (FAIR metadata).

The science is explained in `docs/methods.md`. Steps: `SELECT_SUBLIBRARY` → `SUBSET_LIBRARY` (also whole-library depth
and saturation) → `SCORE` → `AGGREGATE_SYNTHONS`, `VALIDATE`, `MODEL` and `BLOCKSET` per screen → one `BUILDING_BLOCKS`
and one `TRIAGE` over all screens (counter-screens flag the targets' artefacts) → `SHORTLIST` → `VISUALISE`, then
`EXPORT` / `CHECKSUM` / `DESCRIBE` of each screen's top compounds, the triaged hits, the shortlist and the building
blocks, and one `REPORT` and one `EXPLORER`.

Every process runs in the **runtime image** (Podman). The **dev image** is the same base plus
Java, Nextflow, nf-test, just and dev dependencies, so inside the dev container processes run natively.

## Commands (always use `just`)

| Task | Command |
|---|---|
| Format everything | `just fmt` |
| Lint (ruff and `nextflow lint`, including static types) | `just lint` |
| Python tests (100% coverage gate) | `just test-py` |
| nf-test (CI mode, snapshots must match) | `just test-nf` (e.g. `just test-nf --tag score`) |
| Re-record snapshots after an intended change | `just snapshot-update` |
| Every Nextflow file has a test | `just nf-coverage` |
| Everything CI runs (lint, tests, coverage gates) | `just ci` |
| Run on the synthetic example library | `just run` |
| Download KinDEL (~21 GB) / run on it | `just fetch-kindel` / `just run-kindel` |
| Validate FAIR outputs and repository metadata | `just validate-fair` |
| Rebuild README figures / publish report and explorer to `docs/site/` | `just readme-figures` / `just publish-report` |
| Test the explorer page in a headless browser (needs network, not in CI) | `just smoke-explorer results/explorer.html [--chromium PATH]` |
| Where do the sulfonamide blocks rank in each screen (after a KinDEL run) | `just known-chemistry` |

**Before finishing any change, run `just fmt && just ci` and make sure it passes.**
Outside the dev container, `just test-nf` and `just run` build the runtime image and use `-profile podman`.

## Nextflow rules (Nextflow ≥ 26.10, strict syntax)

- Every `.nf` file starts with `nextflow.enable.types = true`. Do not add legacy (untyped) processes or workflows.
- Use **records**, not tuples. Shared record types and helpers live in `common.nf`.
  Re-assemble per-sample data with `join(other, by: 'id')`.
- Process I/O uses typed sections: `input:` (`record(...)`, `name: Type`), `stage:`, `output:`
  (`file()`, `files()`, `eval()`) and `topic:`.
- Pass user-supplied text to commands with `shellQuote()` (from `common.nf`). `stage: env` values are
  **not** quoted inside containers.
- **Versions**: each process emits `tuple('<tool>', eval('<tool> --version')) >> 'versions'` in its
  `topic:` section. A process that consumes the versions file (`REPORT`) must **never** emit to the
  `versions` topic, or the run hangs.
- Use typed operators only: `map`, `filter`, `flatMap`, `join(by:)`, `mix`, `collect`, `groupBy`, `unique`.
  Avoid `branch`, `set`, `tap`, `collectFile`, `splitCsv`, `groupTuple`, `first`.
- Parameters are declared in the typed `params {}` block in `main.nf` and **referenced only in the entry
  workflow**. Pass them on as records (`Options`, `RunInfo`). `params.outdir` is config-only.
- **Never use `publishDir`.** Publish through the entry workflow's `publish:` section and the `output {}` block.
- Every process has a `stub:` block.
- A process call takes **at most one channel argument** (the linter enforces it): collect the others into a `Bag<Path>`
  and address the staged files by name in the script (see `TRIAGE`, `EXPLORER`).
- **A published table's file name must differ from its export stem.** `EXPORT` writes `<stem>.<format>` into a work
  directory that already holds the staged input; the same name clobbers the input (the table is `blocks.parquet`,
  its exports are `building_blocks.*`).
- In `*.nf.test` files, `${...}` inside a `process { """ ... """ }` block is interpolated by nf-test, not Nextflow:
  use constants, or `\${...}`.
- Containers: **Podman only** (`conf/podman.config`, `Containerfile`). Never add files, profiles or images for any other container engine.

## Data rules

- KinDEL data is licensed **CC BY-NC-SA 4.0 (insitro terms), non-commercial**. Never commit it: it lives in
  the gitignored `data/kindel/`. Outputs derived from it carry that license (samplesheet `license` column).
- Test fixtures in `assets/example/` are **synthetic** (`delt_pipeline simulate`, seeded). If you change the
  simulator, regenerate them and re-record snapshots.
- `application/` holds the private CV and motivation letter and is gitignored. Never commit it.
- Library-wide statistics (the depth used for normalisation) come from the **whole** screen, never from
  the sub-library. Priors are fitted on random-BB-A compounds only.

## Python rules

- Keep the CLI thin: `cli.py` parses arguments and calls pure functions in the subpackages.
- Outputs must be **deterministic** so nf-test md5 snapshots are stable: fixed seeds, rounded floats
  (`_numbers.significant`), sorted JSON keys, no timestamps, fixed zip dates.
- Inside the package, import with relative imports (`from ..io import ...`), so the code never needs templating.
- Ruff config is in `pyproject.toml` (line length 120, Google docstrings, type annotations).
- **RDKit and scikit-learn outputs must stay byte-stable**: canonical SMILES only, descriptors rounded
  (`round_sig_figs`), one fingerprint definition in `_chem.py` (Morgan radius 2, 2048 bits, the same as RDKit.js in the
  explorer), seeded models. Model results are floats from a multi-threaded fit: they are listed in `tests/.nftignore`
  and their fields are asserted explicitly instead of by md5.
- **The explorer** (`src/delt_pipeline/explore/`) is one static HTML file. Python builds the payload and the page and is
  covered 100 %. `static/app.js` is glue with no unit tests: after changing it, build a page and run
  `just smoke-explorer`, which checks the page against RDKit and numpy in Python. Keep the page under 12 MB, and keep
  `RDKIT_JS_VERSION` equal to the Python `rdkit` version so both compute the same fingerprints.
- Dependencies are managed with uv: `uv add <pkg>` (runtime) or `uv add --dev <pkg>`, and commit `uv.lock`.

## Tests: 100% is the rule

- **Python**: `pytest` with `--cov-fail-under=100` (line and branch). Every new function needs tests.
- **Nextflow**: every module, workflow and `common.nf` has a `tests/main.nf.test` with snapshots
  (`just nf-coverage` must report 100%). Each module test has a real-data test, a stub test and the
  versions topic (`topics "versions"`, `Fixtures.versions(...)`).
- `tests/lib/Fixtures.groovy` holds shared helpers: linking `assets/` for samplesheets, masking the
  pipeline's own version, canonical ordering of `files()` outputs, and a report summary.
- Files with run-specific content are listed in `tests/.nftignore`; assert their contents explicitly.
- A multi-threaded or timing-dependent output is not md5-snapshotted: add it to `tests/.nftignore` and assert its content.
- Snapshots are made in the runtime image. If you change outputs on purpose, run `just snapshot-update`
  and review every `.snap` diff.

## Adding a process

1. `modules/local/<name>/main.nf`: typed process with `tag`, `label`, `input:`, `output:`, `topic:`, `script:` and `stub:`.
2. If it runs Python, add a subcommand in `cli.py` plus library code and pytest tests.
3. `modules/local/<name>/tests/main.nf.test`: real-data test, stub test and versions.
4. Wire it into `workflows/delt_selection/main.nf` (or a new typed workflow), and publish its outputs in `main.nf`.
5. Describe every new output in the FAIR metadata (see below), then run `just snapshot-update` and `just ci`.

## FAIR rules

- Every published dataset file is listed in `<id>.metadata.jsonld` (schema.org) and, for exports, in
  `datapackage.json` (Frictionless), with its SHA-256 checksum, media type and EDAM format.
- A new export format needs an entry in `src/delt_pipeline/fair/vocab.py` (a test enforces this).
- Never drop licenses (`license` column, `params.data_license`), source references (`source`) or
  creator ORCIDs. Run provenance comes from nf-prov (`ro-crate-metadata.json`, `pipeline_info/bco.json`)
  and Nextflow lineage.
- Keep `CITATION.cff`, `codemeta.json` and `ro-crate-metadata.json` in sync: use `just bump-version X.Y.Z`.

## Data lineage

`lineage.enabled = true`. List runs with `just lineage-list` and inspect outputs with
`just lineage-view lid://<hash>`. Outputs can be reused as inputs by their `lid://` URI
(for example, in a samplesheet `table` column).

## Releasing

`just bump-version X.Y.Z`, commit, then `git tag vX.Y.Z && git push --tags`. The release workflow
checks the versions, pushes `ghcr.io/davidparkerch/delt-pipeline:X.Y.Z`
and creates the GitHub release.
