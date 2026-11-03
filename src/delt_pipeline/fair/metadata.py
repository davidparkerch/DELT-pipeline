"""Dataset metadata: schema.org JSON-LD and Frictionless Data Packages (FAIR F2-F3, I1-I3, R1)."""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import polars as pl

from ..io import FeatureTable
from ..wrangle import normalise_name
from .identifiers import content_urn, sha256_file
from .vocab import FORMATS, format_for, orcid_url, resolvable, spdx_url

JSONLD_CONTEXT = {"@vocab": "https://schema.org/", "prov": "http://www.w3.org/ns/prov#"}
DATAPACKAGE_PROFILE = "data-package"
DICTIONARY_FIELDS = ("description", "unit", "unit_iri", "term_iri")
ROLE_DESCRIPTIONS = {
    "identifier": "Unique row identifier",
    "annotation": "Non-numeric annotation carried along with the data",
    "feature": "Numeric feature",
    "target": "Prediction target (label or response)",
}


@dataclass(frozen=True)
class DatasetInfo:
    """Descriptive metadata supplied with a dataset (usually from the samplesheet)."""

    id: str
    description: str = ""
    source: str = ""
    license: str = "CC-BY-4.0"
    keywords: tuple[str, ...] = ()
    creator_name: str = ""
    creator_orcid: str = ""


@dataclass(frozen=True)
class Provenance:
    """The software that generated the dataset."""

    pipeline: str
    version: str
    homepage: str = ""


@dataclass(frozen=True)
class DatasetDescription:
    """Identifier plus the two metadata documents for one dataset."""

    identifier: str
    jsonld: dict[str, Any]
    datapackage: dict[str, Any]


def read_dictionary(path: Path | None) -> dict[str, dict[str, str]]:
    """Read a data dictionary CSV (``column, description, unit, unit_iri, term_iri``).

    Column names are normalised like the data so the two always match.
    """
    if path is None:
        return {}
    frame = pl.read_csv(path, infer_schema_length=0).fill_null("")
    if "column" not in frame.columns:
        msg = f"data dictionary {path} needs a 'column' column"
        raise ValueError(msg)
    return {
        normalise_name(row["column"]): {key: row[key] for key in DICTIONARY_FIELDS if row.get(key)}
        for row in frame.iter_rows(named=True)
    }


def slug(text: str) -> str:
    """Frictionless-compatible name: lowercase letters, digits, ``-``, ``.`` and ``_``."""
    return re.sub(r"[^a-z0-9._-]+", "-", text.lower()).strip("-") or "resource"


def _roles(table: FeatureTable) -> dict[str, str]:
    roles = {table.id_column: "identifier"}
    roles |= dict.fromkeys(table.annotations, "annotation")
    roles |= dict.fromkeys(table.features, "feature")
    if table.target_column is not None:
        roles[table.target_column] = "target"
    return roles


def frictionless_type(dtype: pl.DataType) -> str:
    """Frictionless Table Schema field type for a polars dtype."""
    if dtype.is_integer():
        return "integer"
    if dtype.is_float():
        return "number"
    if dtype == pl.Boolean:
        return "boolean"
    if dtype == pl.String:
        return "string"
    return "any"


def _variable(name: str, role: str, entry: dict[str, str], summary: dict[str, Any] | None) -> dict[str, Any]:
    variable: dict[str, Any] = {
        "@type": "PropertyValue",
        "name": name,
        "description": entry.get("description") or ROLE_DESCRIPTIONS[role],
    }
    if "unit" in entry:
        variable["unitText"] = entry["unit"]
    if "unit_iri" in entry:
        variable["unitCode"] = entry["unit_iri"]
    if "term_iri" in entry:
        variable["propertyID"] = entry["term_iri"]
    if summary:
        variable["minValue"] = summary["min"]
        variable["maxValue"] = summary["max"]
    return variable


def _schema_field(name: str, dtype: pl.DataType, role: str, entry: dict[str, str]) -> dict[str, Any]:
    schema_field: dict[str, Any] = {
        "name": name,
        "type": frictionless_type(dtype),
        "description": entry.get("description") or ROLE_DESCRIPTIONS[role],
    }
    if "term_iri" in entry:
        schema_field["rdfType"] = entry["term_iri"]
    if "unit" in entry:
        schema_field["unit"] = entry["unit"]
    if role == "identifier":
        schema_field["constraints"] = {"required": True, "unique": True}
    return schema_field


def _file_entry(path: Path, checksums: dict[str, str]) -> dict[str, Any]:
    name = format_for(path.name)
    if name is None:
        msg = f"cannot describe {path.name}: not a known export format"
        raise ValueError(msg)
    info = FORMATS[name]
    return {
        "name": path.name,
        "format": name,
        "media_type": info.media_type,
        "edam": info.edam,
        "documentation": info.documentation,
        "tabular": info.tabular,
        "sha256": checksums.get(path.name) or sha256_file(path),
        "bytes": path.stat().st_size,
    }


def describe_dataset(
    table: FeatureTable,
    clean_path: Path,
    exports: list[Path],
    info: DatasetInfo,
    provenance: Provenance,
    *,
    checksums: dict[str, str] | None = None,
    dictionary: dict[str, dict[str, str]] | None = None,
    feature_stats: list[dict[str, Any]] | None = None,
    exports_base: str = "",
) -> DatasetDescription:
    """Describe a cleaned dataset and its exports.

    Args:
        table: The cleaned table (from ``clean_path``).
        clean_path: The canonical cleaned Parquet file; its hash is the dataset identifier.
        exports: Exported files (same content in other formats).
        info: Descriptive metadata.
        provenance: Generating pipeline.
        checksums: Known SHA-256 digests by file name (e.g. from ``sha256sum``).
        dictionary: Data dictionary from :func:`read_dictionary`.
        feature_stats: Per-feature summaries from the analysis (adds min/max values).
        exports_base: URL prefix from the metadata file to the exports, e.g. ``../../exports/x/``.

    Returns:
        The identifier and both metadata documents.
    """
    checksums = checksums or {}
    dictionary = dictionary or {}
    stats_by_name = {item["feature"]: item for item in feature_stats or []}
    identifier = content_urn(clean_path)
    roles = _roles(table)
    schema = table.frame.schema
    license_url = spdx_url(info.license)

    files = [_file_entry(clean_path, checksums)] + [_file_entry(path, checksums) for path in sorted(exports)]
    distribution = []
    for entry in files:
        base = "" if entry["name"] == clean_path.name else exports_base
        distribution.append(
            {
                "@type": "DataDownload",
                "name": entry["name"],
                "contentUrl": f"{base}{entry['name']}",
                "encodingFormat": [entry["media_type"], entry["edam"]],
                "usageInfo": entry["documentation"],
                "sha256": entry["sha256"],
                "contentSize": f"{entry['bytes']} B",
            }
        )

    jsonld: dict[str, Any] = {
        "@context": JSONLD_CONTEXT,
        "@type": "Dataset",
        "@id": identifier,
        "identifier": [identifier, {"@type": "PropertyValue", "propertyID": "dataset_id", "value": info.id}],
        "name": info.id,
        "description": info.description or f"Cleaned feature table '{info.id}'",
        "keywords": list(info.keywords),
        "license": license_url,
        "variableMeasured": [
            _variable(name, roles[name], dictionary.get(name, {}), stats_by_name.get(name))
            for name in table.frame.columns
        ],
        "distribution": distribution,
        "prov:wasGeneratedBy": {
            "@type": "SoftwareApplication",
            "name": provenance.pipeline,
            "softwareVersion": provenance.version,
            **({"url": provenance.homepage} if provenance.homepage else {}),
        },
    }
    if info.source:
        jsonld["isBasedOn"] = resolvable(info.source)
    if info.creator_name:
        creator: dict[str, Any] = {"@type": "Person", "name": info.creator_name}
        if info.creator_orcid:
            creator["@id"] = orcid_url(info.creator_orcid)
        jsonld["creator"] = creator

    resources = []
    for entry in files[1:]:
        resource: dict[str, Any] = {
            "name": slug(entry["name"]),
            "type": "table" if entry["tabular"] else "file",
            "path": entry["name"],
            "format": entry["format"],
            "mediatype": entry["media_type"],
            "hash": f"sha256:{entry['sha256']}",
            "bytes": entry["bytes"],
        }
        if entry["tabular"]:
            resource["profile"] = "tabular-data-resource"
            resource["encoding"] = "utf-8"
            resource["schema"] = {
                "fields": [
                    _schema_field(name, schema[name], roles[name], dictionary.get(name, {}))
                    for name in table.frame.columns
                ],
                "primaryKey": [table.id_column],
                "missingValues": [""],
            }
        resources.append(resource)

    datapackage: dict[str, Any] = {
        "profile": DATAPACKAGE_PROFILE,
        "id": identifier,
        "name": slug(info.id),
        "title": info.id,
        "description": jsonld["description"],
        "keywords": list(info.keywords),
        "licenses": [{"name": info.license, "path": license_url, "title": info.license}],
        "resources": resources,
    }
    if info.source:
        datapackage["sources"] = [{"title": info.source, "path": resolvable(info.source)}]
    if info.creator_name:
        contributor = {"title": info.creator_name, "role": "author"}
        if info.creator_orcid:
            contributor["path"] = orcid_url(info.creator_orcid)
        datapackage["contributors"] = [contributor]
    return DatasetDescription(identifier, jsonld, datapackage)
