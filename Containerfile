# DELT-pipeline images (build with Podman):
#   podman build --target runtime -t localhost/delt-pipeline:dev .      # runs pipeline processes
#   podman build --target dev -t localhost/delt-pipeline-dev:latest .   # development environment
# Both share the `base` stage, so the dev image can run every process natively.

FROM ghcr.io/astral-sh/uv:0.13.6-python3.14-trixie-slim AS base

ARG VERSION=0.1.0
ARG REVISION=unknown
LABEL org.opencontainers.image.title="delt-pipeline" \
      org.opencontainers.image.description="Process DNA-Encoded Library Technology data and exports FAIR, analysis-ready datasets." \
      org.opencontainers.image.source="https://github.com/davidparkerch/delt-pipeline" \
      org.opencontainers.image.url="https://github.com/davidparkerch/delt-pipeline" \
      org.opencontainers.image.licenses="Apache-2.0" \
      org.opencontainers.image.authors="David Parker" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${REVISION}"

# procps: Nextflow task metrics; gzip and coreutils: the Linux tools used by the pipeline
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps gzip coreutils \
    && rm -rf /var/lib/apt/lists/*

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    PATH=/opt/venv/bin:$PATH

WORKDIR /opt/delt-pipeline

# ------------------------------------------------------------------ runtime
FROM base AS runtime

# Dependencies first (cached layer), then the project itself (non-editable).
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-dev
COPY pyproject.toml uv.lock README.md ./
COPY src ./src
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-editable \
    && delt_pipeline --version

# ------------------------------------------------------------------ dev
FROM base AS dev

ARG NXF_VER=26.10.5
ARG NF_TEST_VERSION=0.9.5
ARG USERNAME=dev
ARG UID=1000
ARG GID=1000

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash-completion ca-certificates curl git just less openjdk-21-jre-headless \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid ${GID} ${USERNAME} \
    && useradd --uid ${UID} --gid ${GID} --create-home --shell /bin/bash ${USERNAME} \
    && echo "export PATH=/home/${USERNAME}/.local/bin:/opt/venv/bin:\$PATH" > /etc/profile.d/delt-pipeline.sh

# Project dependencies including the dev group; the project itself is installed
# editable from the mounted workspace by `uv sync` (devcontainer postCreateCommand).
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project \
    && chown -R ${UID}:${GID} /opt/venv

USER ${USERNAME}
ENV NXF_VER=${NXF_VER} \
    PATH=/home/${USERNAME}/.local/bin:/opt/venv/bin:$PATH

# Nextflow and nf-test (pinned), with the Nextflow runtime and plugins pre-fetched
RUN mkdir -p ~/.local/bin && cd ~/.local/bin \
    && curl -fsSL https://get.nextflow.io | bash \
    && nextflow -version \
    && nextflow plugin install nf-prov@1.8.0 \
    && curl -fsSL https://get.nf-test.com | bash -s ${NF_TEST_VERSION} \
    && nf-test version

WORKDIR /workspace
