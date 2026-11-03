"""Content-derived identifiers and checksums (FAIR F1, F3)."""

import hashlib
from pathlib import Path

_CHUNK = 1 << 20


def sha256_file(path: Path) -> str:
    """Hex SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def content_urn(path: Path) -> str:
    """Globally unique, persistent identifier for a file's content: ``urn:sha256:<hex>``."""
    return f"urn:sha256:{sha256_file(path)}"


def parse_sha256sums(path: Path) -> dict[str, str]:
    """Parse ``sha256sum`` output (``<hex>  <name>`` per line) into ``{name: hex}``.

    Raises:
        ValueError: On a malformed line.
    """
    checksums = {}
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        digest, separator, name = line.partition(" ")
        if not separator or len(digest) != 64:  # noqa: PLR2004 - SHA-256 hex length
            msg = f"{path}:{number}: not a sha256sum line: {line!r}"
            raise ValueError(msg)
        checksums[name.lstrip(" *")] = digest
    return checksums
