# DELT-pipeline: development commands. Run `just` to list them.
#
# Inside the dev container, pipeline processes run natively (the dev image is a superset of
# the runtime image). Anywhere else they run in the runtime image with Podman.

set shell := ["bash", "-euo", "pipefail", "-c"]
set dotenv-load

export NXF_VER := "26.10.5"

runtime_image := "localhost/delt-pipeline:dev"
dev_image := "localhost/delt-pipeline-dev:latest"
roc_validator := "roc-validator==0.13.1"
cffconvert := "cffconvert==2.0.0"
nf_core := "nf-core==4.2.1"

# Nextflow sources checked by `nextflow lint` (results/ and work/ are never linted)
nf_sources := "main.nf common.nf nextflow.config conf modules workflows"
in_container := path_exists("/run/.containerenv")

# List recipes
default:
    @just --list --unsorted

# ---------------------------------------------------------------- setup

# Create or refresh the Python environment (uv)
setup:
    uv sync

# ---------------------------------------------------------------- quality

# Format Python and Nextflow code in place
fmt:
    uv run ruff format .
    uv run ruff check --fix .
    nextflow lint -format -spaces 4 {{ nf_sources }}

# Run every linter without changing files
lint: lint-py lint-nf

# Ruff lint and format check
lint-py:
    uv run ruff check .
    uv run ruff format --check .

# Nextflow syntax, static types and formatting
lint-nf:
    #!/usr/bin/env bash
    set -euo pipefail
    nextflow lint -o concise {{ nf_sources }}
    # `nextflow lint -format` has no check mode: format a copy and compare.
    scratch="$(mktemp -d)"
    trap 'rm -rf "$scratch"' EXIT
    cp -r {{ nf_sources }} "$scratch"/
    (cd "$scratch" && nextflow lint -format -spaces 4 {{ nf_sources }} > /dev/null)
    for path in {{ nf_sources }}; do
        diff -ru "$path" "$scratch/$path" || { echo "Nextflow code is not formatted: run 'just fmt'" >&2; exit 1; }
    done

# ---------------------------------------------------------------- tests

# Lint, every test and the test-coverage gates (what CI runs)
ci: lint test nf-coverage

# Python and Nextflow tests
test: test-py test-nf

# pytest with the 100% coverage gate
test-py *args:
    uv run pytest {{ args }}

# nf-test in CI mode (fails on missing or changed snapshots), e.g. `just test-nf --tag wrangle`
test-nf *args: _runtime
    {{ _nf_env }} nf-test test --ci {{ _nf_profile }} {{ args }}

# Re-record nf-test snapshots after an intended change; review the diff before committing
snapshot-update *args: _runtime
    {{ _nf_env }} nf-test test --update-snapshot {{ _nf_profile }} {{ args }}

# Which Nextflow files have tests (must stay at 100%)
nf-coverage:
    #!/usr/bin/env bash
    set -euo pipefail
    report="$(nf-test coverage --silent)"
    echo "$report"
    grep -q '(100%)' <<< "$report"

# ---------------------------------------------------------------- run

# Run the pipeline on the bundled example data, e.g. `just run --export_formats parquet,npz`
run *args: _runtime
    {{ _nf_env }} nextflow run . -profile test{{ if in_container == "true" { "" } else { ",podman" } }} {{ args }}

# Open the HTML report from the last run
report dir="results":
    uv run python -m webbrowser "file://$(realpath {{ dir }}/report.html)"

# Validate the FAIR metadata of a run (RO-Crate, data packages, identifiers) and of the repository
validate-fair dir="results":
    #!/usr/bin/env bash
    set -euo pipefail
    # Inputs are referenced, not copied, so the payload-availability check (ro-crate-1.1_12.1) is skipped.
    uvx --from {{ roc_validator }} rocrate-validator validate -p provenance-run-crate-0.5 \
        --skip-checks ro-crate-1.1_12.1 --no-paging "{{ dir }}"
    for package in "{{ dir }}"/exports/*/datapackage.json; do
        uv run frictionless validate "$package"
    done
    uv run python - "{{ dir }}" <<'PY'
    import hashlib, json, pathlib, sys
    for path in sorted(pathlib.Path(sys.argv[1]).glob("datasets/*/*.metadata.jsonld")):
        metadata = json.loads(path.read_text())
        clean = path.with_name(metadata["distribution"][0]["name"])
        expected = "urn:sha256:" + hashlib.sha256(clean.read_bytes()).hexdigest()
        assert metadata["@id"] == expected, f"{path}: identifier does not match {clean.name}"
        print(f"{path}: identifier ok")
    PY
    just validate-repo

# Validate the repository's own metadata: citation, workflow RO-Crate and versions
validate-repo:
    uvx --from {{ cffconvert }} cffconvert --validate -i CITATION.cff
    uvx --from {{ roc_validator }} rocrate-validator validate -p workflow-ro-crate-1.0 --no-paging .
    just release-check

# ---------------------------------------------------------------- lineage

# List runs recorded in the lineage store
lineage-list:
    nextflow lineage list

# Show the lineage record of a lid:// identifier
lineage-view lid:
    nextflow lineage view {{ lid }}

# Render the lineage graph of an output as lineage.html
lineage-render lid:
    nextflow lineage render {{ lid }}

# ---------------------------------------------------------------- containers

# Build the runtime image used by pipeline processes
build-runtime:
    podman build --target runtime --build-arg REVISION="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)" -t {{ runtime_image }} .

# Build the development image
build-dev:
    podman build --target dev -t {{ dev_image }} .

# Interactive shell in the dev container (workspace mounted at /workspace)
dev-shell: (dev-run "bash")

# Run a command in the dev container, e.g. `just dev-run just ci`
dev-run +command:
    podman run --rm $([[ -t 0 ]] && echo --interactive --tty) --userns=keep-id:uid=1000,gid=1000 \
        -v "$PWD":/workspace:Z -w /workspace {{ dev_image }} \
        bash -c 'uv sync --quiet && {{ command }}'

# ---------------------------------------------------------------- release

# Set a new version everywhere it is recorded, e.g. `just bump-version 0.2.0`
bump-version new:
    #!/usr/bin/env bash
    set -euo pipefail
    old="$(uv version --short --color never)"
    uv version "{{ new }}"
    sed -i "s/version = '$old'/version = '{{ new }}'/" nextflow.config
    sed -i "s/^version: $old$/version: {{ new }}/" CITATION.cff
    sed -i "s/\"version\": \"$old\"/\"version\": \"{{ new }}\"/" codemeta.json
    sed -i "s/\"version\": \"$old\"/\"version\": \"{{ new }}\"/" ro-crate-metadata.json
    sed -i "s/:$old'/:{{ new }}'/" conf/podman.config
    sed -i "s/^ARG VERSION=$old$/ARG VERSION={{ new }}/" Containerfile
    just release-check
    echo "Bumped $old -> {{ new }}. Commit, then tag v{{ new }} to release."

# Check that every recorded version agrees (and matches a tag, if given)
release-check tag="":
    #!/usr/bin/env bash
    set -euo pipefail
    version="$(uv version --short --color never)"
    status=0
    check() {
        if ! grep -qF -- "$2" "$1"; then echo "✗ $1 does not contain: $2" >&2; status=1; else echo "✓ $1"; fi
    }
    check nextflow.config "version = '$version'"
    check CITATION.cff "version: $version"
    check codemeta.json "\"version\": \"$version\""
    check ro-crate-metadata.json "\"version\": \"$version\""
    check conf/podman.config ":$version'"
    check Containerfile "ARG VERSION=$version"
    if [[ -n "{{ tag }}" && "{{ tag }}" != "v$version" ]]; then
        echo "✗ tag {{ tag }} does not match version $version" >&2
        status=1
    fi
    exit $status

# ---------------------------------------------------------------- utilities

# nf-core tools, e.g. `just nf-core test-datasets list`
nf-core *args:
    uvx --from {{ nf_core }} nf-core {{ args }}

# Remove run, test and cache artefacts
clean:
    rm -rf work results .nextflow .nextflow.log* .nf-test .nf-test.log .lineage .pytest_cache .ruff_cache coverage.xml lineage.html

# ---------------------------------------------------------------- internals

# Podman needs the runtime image unless we are inside the dev container
_runtime:
    {{ if in_container == "true" { "true" } else { "just build-runtime" } }}

_nf_env := if in_container == "true" { "" } else { "PIPELINE_IMAGE=" + runtime_image }
_nf_profile := if in_container == "true" { "" } else { "--profile +podman" }

