"""The data behind the explorer page, as compact JSON.

The page has no server, so everything it shows travels in the HTML: the triaged hits of the target
screens (one row array per hit, columns named once), Morgan fingerprints of the hits for similarity
search, the building blocks with their descriptors, and each block's set enrichment in every screen.
Numbers are rounded and rows are arrays, which keeps two screens of 2,000 hits at about 3 MB.
"""

import base64
from typing import Any

import numpy as np
import polars as pl

from .._chem import FP_BITS, morgan_generator, parse_smiles
from .._numbers import significant

DIGITS = 4
"""Significant digits of the numbers in the page."""
HIT_COLUMNS = (
    "hit_id",
    "screen",
    "rank",
    "tier",
    "molecule_hash",
    "smiles",
    "bb_a",
    "bb_b",
    "bb_c",
    "target_1",
    "target_2",
    "target_3",
    "matrix_1",
    "matrix_2",
    "matrix_3",
    "load",
    "target_sum",
    "control_sum",
    "target_support",
    "count_diff",
    "er",
    "poisson_ratio",
    "poisson_lower",
    "poisson_upper",
    "z_faver",
    "eb_log2_ratio",
    "matrix_lower",
    "counter_lower",
    "disynthon",
    "disynthon_lower",
    "replicate_support",
    "disynthon_support",
    "counter_binder",
    "matrix_binder",
    "shared_target",
    "shared_with",
    "reasons",
    "pick_rank",
    "cluster_id",
)
"""Columns of a hit row, in order; the page reads them by name from the payload."""
BLOCK_COLUMNS = (
    "block_id",
    "position",
    "bb_id",
    "fragment_smiles",
    "mw",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "heavy_atoms",
    "rings",
    "murcko_scaffold",
    "n_same_fragment",
    "mean_load",
    "low_load",
)
"""Columns of a building-block row, in order."""
AUC_DECIMALS = 3
BLOCKSET_COLUMNS = ("auc", "q_value", "n_top")
"""What the page knows about a block in a screen."""


def fingerprint_bytes(smiles: list[str]) -> list[bytes]:
    """Morgan fingerprints (radius 2, 2048 bits) packed to 256 bytes each, bit 0 first in each byte's top bit.

    The same fingerprint RDKit.js computes in the browser, so a query can be compared with every hit there.

    Raises:
        ValueError: If a SMILES cannot be parsed.
    """
    generator = morgan_generator()
    packed = []
    for text in smiles:
        bits = np.zeros(FP_BITS, dtype=np.uint8)
        bits[list(generator.GetFingerprint(parse_smiles(text)).GetOnBits())] = 1
        packed.append(np.packbits(bits).tobytes())
    return packed


def _clean(value: Any) -> Any:  # noqa: ANN401 - a cell of any type
    """Round floats; turn booleans into 0/1 and non-finite numbers into null."""
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, float):
        return significant(value, DIGITS)
    return value


def hit_rows(hits: pl.DataFrame, picks: pl.DataFrame) -> list[list[Any]]:
    """One row array per hit, in :data:`HIT_COLUMNS` order.

    ``pick_rank`` and ``cluster_id`` come from the shortlist (0 and null when the hit was not picked).
    """
    joined = hits.join(
        picks.select("hit_id", "pick_rank", "cluster_id"), on="hit_id", how="left", maintain_order="left"
    ).with_columns(pl.col("pick_rank").fill_null(0))
    return [
        [_clean(row[name]) for name in HIT_COLUMNS]
        for row in joined.sort("screen", "rank").select(HIT_COLUMNS).iter_rows(named=True)
    ]


def block_rows(blocks: pl.DataFrame) -> list[list[Any]]:
    """One row array per building block, in :data:`BLOCK_COLUMNS` order."""
    return [
        [_clean(row[name]) for name in BLOCK_COLUMNS]
        for row in blocks.sort("position", "bb_id").select(BLOCK_COLUMNS).iter_rows(named=True)
    ]


def blockset_map(blocksets: dict[str, pl.DataFrame]) -> dict[str, dict[str, list[Any]]]:
    """Set-enrichment results per screen: ``{screen: {block_id: [auc, q_value, n_top]}}``."""
    result: dict[str, dict[str, list[Any]]] = {}
    for screen, table in sorted(blocksets.items()):
        keyed = table.with_columns(
            block_id=pl.col("position").str.to_uppercase() + pl.col("bb_id").cast(pl.String).str.zfill(4)
        )
        result[screen] = {
            # The AUC is rounded once, to the 3 decimals the page shows, so Python and the page agree.
            row["block_id"]: [round(row["auc"], AUC_DECIMALS), _clean(row["q_value"]), _clean(row["n_top"])]
            for row in keyed.sort("block_id").iter_rows(named=True)
        }
    return result


def build_payload(
    hits: pl.DataFrame,
    picks: pl.DataFrame,
    blocks: pl.DataFrame,
    blocksets: dict[str, pl.DataFrame],
    meta: dict[str, Any],
) -> dict[str, Any]:
    """Everything the page needs.

    Args:
        hits: The triage table (``hits.parquet``).
        picks: The resynthesis shortlist.
        blocks: The building-block resource.
        blocksets: Set-enrichment table of every screen, by screen id (counter-screens included).
        meta: Free-form facts for the page (title, license, links, parameters).

    Returns:
        A JSON-serialisable dict. ``fingerprints`` holds one base64 string per hit row.
    """
    rows = hit_rows(hits, picks)
    smiles = [row[HIT_COLUMNS.index("smiles")] for row in rows]
    return {
        "meta": meta,
        "screens": sorted(hits.get_column("screen").unique().to_list()),
        "hit_columns": list(HIT_COLUMNS),
        "hits": rows,
        "fingerprints": [base64.b64encode(packed).decode("ascii") for packed in fingerprint_bytes(smiles)],
        "block_columns": list(BLOCK_COLUMNS),
        "blocks": block_rows(blocks),
        "blocksets": blockset_map(blocksets),
        "blockset_columns": list(BLOCKSET_COLUMNS),
    }
