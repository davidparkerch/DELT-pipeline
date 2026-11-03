"""Shared vocabularies for formats and licenses (FAIR I1, I2, R1.1).

Every export format in :data:`..io.EXPORT_FORMATS` must be described here; a test enforces it.
"""

from dataclasses import dataclass

from ..io import EXPORT_FORMATS

EDAM = "http://edamontology.org/"
EDAM_TEXT = f"{EDAM}format_2330"
EDAM_BINARY = f"{EDAM}format_2333"
EDAM_CSV = f"{EDAM}format_3752"
EDAM_TSV = f"{EDAM}format_3475"
EDAM_JSON = f"{EDAM}format_3464"
EDAM_JSONLD = f"{EDAM}format_3749"
EDAM_NUMPY = f"{EDAM}format_4003"


@dataclass(frozen=True)
class FormatInfo:
    """How a file format is described in metadata."""

    media_type: str
    """IANA (or widely used) media type."""
    edam: str
    """EDAM format IRI (the closest generic term when EDAM has no specific one)."""
    documentation: str
    """URL of the format specification."""
    tabular: bool
    """Whether Frictionless can read it as a table (adds a Table Schema to the data package)."""


FORMATS: dict[str, FormatInfo] = {
    "parquet": FormatInfo(
        "application/vnd.apache.parquet", EDAM_BINARY, "https://parquet.apache.org/docs/file-format/", tabular=False
    ),
    "arrow": FormatInfo(
        "application/vnd.apache.arrow.file", EDAM_BINARY, "https://arrow.apache.org/docs/format/Columnar.html", False
    ),
    "csv": FormatInfo("text/csv", EDAM_CSV, "https://www.rfc-editor.org/rfc/rfc4180", tabular=True),
    "tsv": FormatInfo(
        "text/tab-separated-values",
        EDAM_TSV,
        "https://www.iana.org/assignments/media-types/text/tab-separated-values",
        True,
    ),
    "jsonl": FormatInfo("application/jsonl", EDAM_TEXT, "https://jsonlines.org/", tabular=True),
    "npz": FormatInfo(
        "application/zip", EDAM_NUMPY, "https://numpy.org/doc/stable/reference/generated/numpy.savez.html", False
    ),
    "svmlight": FormatInfo("text/plain", EDAM_TEXT, "https://www.cs.cornell.edu/people/tj/svm_light/", tabular=False),
}


def format_for(filename: str) -> str | None:
    """Export format name for a file name, or ``None`` if it is not an export."""
    for name, export in EXPORT_FORMATS.items():
        if filename.endswith(export.suffix):
            return name
    return None


def spdx_url(license_id: str) -> str:
    """Resolvable URL for an SPDX license identifier (already-URLs pass through)."""
    if license_id.startswith(("http://", "https://")):
        return license_id
    return f"https://spdx.org/licenses/{license_id}"


def resolvable(reference: str) -> str:
    """Turn a DOI, ``doi:`` CURIE or URL into a resolvable URL (FAIR A1)."""
    reference = reference.strip()
    if reference.lower().startswith("doi:"):
        return f"https://doi.org/{reference[4:]}"
    if reference.startswith("10."):
        return f"https://doi.org/{reference}"
    return reference


def orcid_url(orcid: str) -> str:
    """Resolvable ORCID URL from a bare iD or URL."""
    orcid = orcid.strip()
    return orcid if orcid.startswith("https://") else f"https://orcid.org/{orcid}"
