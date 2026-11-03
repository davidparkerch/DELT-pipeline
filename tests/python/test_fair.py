import hashlib
import json
from pathlib import Path

import polars as pl
import pytest
from frictionless import Package

from delt_pipeline.fair import (
    FORMATS,
    DatasetInfo,
    Provenance,
    content_urn,
    describe_dataset,
    format_for,
    frictionless_type,
    orcid_url,
    parse_sha256sums,
    read_dictionary,
    resolvable,
    sha256_file,
    slug,
    spdx_url,
)
from delt_pipeline.io import EXPORT_FORMATS, FeatureTable, export_table

PROVENANCE = Provenance("demo-pipeline", "1.2.3", "https://example.org/demo")


def test_identifiers(tmp_path: Path) -> None:
    path = tmp_path / "data.bin"
    path.write_bytes(b"fair")
    digest = hashlib.sha256(b"fair").hexdigest()
    assert sha256_file(path) == digest
    assert content_urn(path) == f"urn:sha256:{digest}"


def test_parse_sha256sums(tmp_path: Path) -> None:
    digest = "a" * 64
    path = tmp_path / "SHA256SUMS"
    path.write_text(f"{digest}  one.csv\n\n{digest} *two.npz\n")
    assert parse_sha256sums(path) == {"one.csv": digest, "two.npz": digest}
    path.write_text("not-a-checksum\n")
    with pytest.raises(ValueError, match="not a sha256sum line"):
        parse_sha256sums(path)


def test_every_export_format_is_described() -> None:
    assert set(FORMATS) == set(EXPORT_FORMATS)
    assert all(info.edam.startswith("http://edamontology.org/format_") for info in FORMATS.values())


def test_vocab_helpers() -> None:
    assert format_for("x.clean.parquet") == "parquet"
    assert format_for("SHA256SUMS") is None
    assert spdx_url("CC0-1.0") == "https://spdx.org/licenses/CC0-1.0"
    assert spdx_url("https://example.org/license") == "https://example.org/license"
    assert resolvable("doi:10.1000/xyz") == "https://doi.org/10.1000/xyz"
    assert resolvable(" 10.1000/xyz ") == "https://doi.org/10.1000/xyz"
    assert resolvable("https://example.org/data") == "https://example.org/data"
    assert orcid_url("0000-0002-1825-0097") == "https://orcid.org/0000-0002-1825-0097"
    assert orcid_url("https://orcid.org/0000-0002-1825-0097") == "https://orcid.org/0000-0002-1825-0097"
    assert slug("My Data_set!") == "my-data_set"
    assert slug("!!!") == "resource"


@pytest.mark.parametrize(
    ("dtype", "expected"),
    [(pl.Int64, "integer"), (pl.Float32, "number"), (pl.Boolean, "boolean"), (pl.String, "string"), (pl.Date, "any")],
)
def test_frictionless_type(dtype: pl.DataType, expected: str) -> None:
    assert frictionless_type(dtype) == expected


def test_read_dictionary(tmp_path: Path) -> None:
    assert read_dictionary(None) == {}
    path = tmp_path / "dictionary.csv"
    path.write_text("column,description,unit,unit_iri,term_iri\nfeat alpha,Alpha,Da,http://u/1,\n")
    assert read_dictionary(path) == {"feat_alpha": {"description": "Alpha", "unit": "Da", "unit_iri": "http://u/1"}}
    path.write_text("name,description\nx,y\n")
    with pytest.raises(ValueError, match="needs a 'column' column"):
        read_dictionary(path)


@pytest.fixture
def described(tmp_path: Path, table: FeatureTable):
    clean = tmp_path / "ds.clean.parquet"
    table.frame.write_parquet(clean)
    exports = export_table(table, "ds", list(EXPORT_FORMATS), tmp_path)
    info = DatasetInfo(
        id="ds",
        description="Demo set",
        source="doi:10.1000/demo",
        license="CC0-1.0",
        keywords=("chemistry", "ml"),
        creator_name="Ada Lovelace",
        creator_orcid="0000-0002-1825-0097",
    )
    dictionary = {
        "feat_alpha": {"description": "Alpha", "unit": "Da", "unit_iri": "http://u/1", "term_iri": "http://t/1"}
    }
    stats = [{"feature": "feat_alpha", "min": -3.0, "max": 3.0}]
    checksums = {exports[0].name: sha256_file(exports[0])}
    description = describe_dataset(
        table, clean, exports, info, PROVENANCE, checksums=checksums, dictionary=dictionary, feature_stats=stats
    )
    return description, clean, exports


def test_jsonld(described) -> None:
    description, clean, exports = described
    jsonld = description.jsonld
    assert description.identifier == jsonld["@id"] == content_urn(clean)
    assert jsonld["isBasedOn"] == "https://doi.org/10.1000/demo"
    assert jsonld["license"] == "https://spdx.org/licenses/CC0-1.0"
    assert jsonld["creator"]["@id"] == "https://orcid.org/0000-0002-1825-0097"
    assert jsonld["prov:wasGeneratedBy"]["url"] == PROVENANCE.homepage
    alpha = next(item for item in jsonld["variableMeasured"] if item["name"] == "feat_alpha")
    assert alpha == {
        "@type": "PropertyValue",
        "name": "feat_alpha",
        "description": "Alpha",
        "unitText": "Da",
        "unitCode": "http://u/1",
        "propertyID": "http://t/1",
        "minValue": -3.0,
        "maxValue": 3.0,
    }
    assert {item["name"] for item in jsonld["distribution"]} == {clean.name} | {path.name for path in exports}
    assert all(item["sha256"] == sha256_file(clean.parent / item["name"]) for item in jsonld["distribution"])


def test_datapackage_validates(described, tmp_path: Path) -> None:
    description, _, exports = described
    descriptor = tmp_path / "datapackage.json"
    descriptor.write_text(json.dumps(description.datapackage))
    report = Package(str(descriptor)).validate()
    assert report.valid, report.flatten(["type", "note"])
    tabular = [item for item in description.datapackage["resources"] if "schema" in item]
    assert {item["format"] for item in tabular} == {"csv", "tsv", "jsonl"}
    assert len(description.datapackage["resources"]) == len(exports)
    assert description.datapackage["contributors"][0]["path"].startswith("https://orcid.org/")


def test_minimal_description(tmp_path: Path, unlabelled: FeatureTable) -> None:
    clean = tmp_path / "u.clean.parquet"
    unlabelled.frame.write_parquet(clean)
    description = describe_dataset(
        unlabelled,
        clean,
        [],
        DatasetInfo(id="u", creator_name="Anon"),
        Provenance("p", "0.0.1"),
        exports_base="../x/",
    )
    assert "isBasedOn" not in description.jsonld
    assert "@id" not in description.jsonld["creator"]
    assert "url" not in description.jsonld["prov:wasGeneratedBy"]
    assert "sources" not in description.datapackage
    assert "path" not in description.datapackage["contributors"][0]
    anonymous = describe_dataset(unlabelled, clean, [], DatasetInfo(id="u"), Provenance("p", "0.0.1"))
    assert "creator" not in anonymous.jsonld
    assert "contributors" not in anonymous.datapackage


def test_describe_rejects_unknown_files(tmp_path: Path, unlabelled: FeatureTable) -> None:
    clean = tmp_path / "u.clean.parquet"
    unlabelled.frame.write_parquet(clean)
    other = tmp_path / "notes.txt"
    other.write_text("x")
    with pytest.raises(ValueError, match="not a known export format"):
        describe_dataset(unlabelled, clean, [other], DatasetInfo(id="u"), PROVENANCE)
