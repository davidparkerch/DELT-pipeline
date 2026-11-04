"""Sub-library selection and the compact count table, streamed from a full KinDEL screen.

The full KinDEL tables are ~7 GB each (68 M compounds). They are never loaded into memory: polars scans
the Parquet file lazily, once to total every count column over the **whole** library (the depth that
normalises rates), and once to keep the chosen sub-library.
"""

from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

from .depth import depth_qc
from .schema import BB_COLUMNS, COUNTS, ID, KINDEL_CONTROL, KINDEL_COUNTS, KINDEL_TARGET, POSITIONS, RENAME


def read_building_blocks(directory: Path) -> dict[str, pl.DataFrame]:
    """Read ``reacted_{a,b,c}.csv`` (KinDEL synthon tables); the row number is the building-block id.

    Returns:
        One frame per position with columns ``bb_<p>`` (UInt16 id) and ``smiles_<p>``.

    Raises:
        ValueError: If a table lacks its ``smiles_<p>`` column or lists a SMILES twice.
    """
    tables = {}
    for position in POSITIONS:
        path = directory / f"reacted_{position}.csv"
        column = f"smiles_{position}"
        frame = pl.read_csv(path, infer_schema_length=0)
        if column not in frame.columns:
            msg = f"{path} has no '{column}' column"
            raise ValueError(msg)
        if frame.get_column(column).is_duplicated().any():
            msg = f"{path} lists a building block SMILES more than once"
            raise ValueError(msg)
        tables[position] = frame.select(pl.int_range(pl.len(), dtype=pl.UInt16).alias(f"bb_{position}"), pl.col(column))
    return tables


def library_diversity(building_blocks: dict[str, pl.DataFrame]) -> int:
    """Number of compounds the combinatorial library encodes (product of building-block counts)."""
    return int(np.prod([building_blocks[position].height for position in POSITIONS]))


def heldout_bb_a(building_blocks: dict[str, pl.DataFrame], heldout: list[pl.DataFrame]) -> tuple[list[int], int]:
    """BB-A ids of every in-library held-out compound.

    Returns:
        ``(ids, unmatched)``: sorted unique ids, and how many in-library compounds had an unknown BB-A.
    """
    lookup = building_blocks["a"]
    frames = [frame.select(pl.col("smiles_a").cast(pl.String)) for frame in heldout]
    smiles = pl.concat(frames).drop_nulls() if frames else pl.DataFrame(schema={"smiles_a": pl.String})
    matched = smiles.join(lookup, on="smiles_a", how="left")
    ids = sorted(set(matched.get_column("bb_a").drop_nulls().to_list()))
    return ids, int(matched.get_column("bb_a").null_count())


def select_sublibrary(
    building_blocks: dict[str, pl.DataFrame], heldout: list[pl.DataFrame], n_bb_a: int, seed: int
) -> dict[str, Any]:
    """Choose the BB-A building blocks of the sub-library.

    All BB-As of held-out compounds are kept (so they can be validated), then random BB-As are added
    until there are ``n_bb_a``. The two groups are reported separately: the random group is an
    unbiased sample of the library and is what library-wide priors are fitted on.

    Args:
        building_blocks: From :func:`read_building_blocks`.
        heldout: Held-out KD tables (KinDEL ``heldout/*.csv``); may be empty.
        n_bb_a: Target number of BB-As (at least the held-out ones are kept).
        seed: Seed for the random BB-As.

    Returns:
        JSON-serialisable description with ``heldout`` and ``random`` id lists.
    """
    if n_bb_a < 1:
        msg = f"n_bb_a must be at least 1, got {n_bb_a}"
        raise ValueError(msg)
    total = building_blocks["a"].height
    from_heldout, unmatched = heldout_bb_a(building_blocks, heldout)
    remaining = np.setdiff1d(np.arange(total), from_heldout)
    n_random = min(max(n_bb_a - len(from_heldout), 0), remaining.size)
    chosen = np.random.default_rng(seed).choice(remaining, size=n_random, replace=False)
    sizes = {position: building_blocks[position].height for position in POSITIONS}
    return {
        "building_blocks": sizes,
        "diversity": library_diversity(building_blocks),
        "heldout": from_heldout,
        "random": sorted(int(index) for index in chosen),
        "unmatched_heldout": unmatched,
        "seed": seed,
    }


def library_depth(table: Path) -> dict[str, Any]:
    """Total of every count column, and the number of compounds, over the whole screen (streamed)."""
    totals = (
        pl.scan_parquet(table)
        .select(pl.len().alias("n_compounds"), *(pl.col(name).sum() for name in KINDEL_COUNTS))
        .collect(engine="streaming")
        .row(0, named=True)
    )
    return {
        "n_compounds": int(totals["n_compounds"]),
        "depth": {RENAME[name]: int(totals[name]) for name in KINDEL_COUNTS},
    }


def count_histograms(table: Path) -> dict[str, dict[int, int]]:
    """Number of compounds with exactly ``k >= 1`` counts, for every count column (streamed).

    The pooled target and matrix counts (``target``, ``matrix``) are included because the scores use
    them. Each column is scanned on its own, so only that column is read from the Parquet file and the
    group-by holds one entry per distinct count.
    """
    columns = {RENAME[name]: pl.col(name) for name in KINDEL_COUNTS}
    columns["target"] = pl.sum_horizontal(KINDEL_TARGET)
    columns["matrix"] = pl.sum_horizontal(KINDEL_CONTROL)
    histograms = {}
    for name, expression in columns.items():
        frame = (
            pl.scan_parquet(table)
            .select(expression.alias("count"))
            .filter(pl.col("count") > 0)
            .group_by("count")
            .len()
            .collect(engine="streaming")
        )
        histograms[name] = dict(
            zip(frame.get_column("count").to_list(), frame.get_column("len").to_list(), strict=True)
        )
    return histograms


def _bb_id(building_blocks: dict[str, pl.DataFrame], position: str) -> pl.Expr:
    """Map a synthon SMILES column to its building-block id (null when unknown), element by element."""
    table = building_blocks[position]
    return pl.col(f"smiles_{position}").replace_strict(
        table.get_column(f"smiles_{position}"), table.get_column(f"bb_{position}"), default=None, return_dtype=pl.UInt16
    )


def subset_library(table: Path, building_blocks: dict[str, pl.DataFrame], bb_a: list[int]) -> pl.LazyFrame:
    """Lazily select the sub-library and map synthon SMILES to building-block ids.

    The BB-A filter comes first, so polars decodes the wide string columns only for kept rows; the
    mapping is element-wise (no joins), so the query streams in constant memory and keeps the
    screen's row order. Compounds with an unknown (blank) building block get a null id.
    """
    keep = building_blocks["a"].filter(pl.col("bb_a").is_in(bb_a)).get_column("smiles_a")
    return (
        pl.scan_parquet(table)
        .filter(pl.col("smiles_a").is_in(keep.implode()))
        .select(
            ID,
            "smiles",
            *(_bb_id(building_blocks, position).alias(f"bb_{position}") for position in POSITIONS),
            *(pl.col(name).cast(pl.Int32).alias(RENAME[name]) for name in KINDEL_COUNTS),
        )
    )


def write_counts(
    table: Path, building_blocks: dict[str, pl.DataFrame], sublibrary: dict[str, Any], output: Path
) -> dict[str, Any]:
    """Stream the compact sub-library count table to Parquet and return its QC summary.

    Nothing larger than a few row groups is held in memory: the table is written by a streaming
    sink, and the QC is computed by scanning the written file.

    Args:
        table: Full KinDEL screen (Parquet).
        building_blocks: From :func:`read_building_blocks`.
        sublibrary: From :func:`select_sublibrary`.
        output: Parquet file to write.

    Returns:
        QC: whole-library depth (used for normalisation) and saturation (counts per library member,
        fraction of the library observed, rarefaction), sub-library size, depth and sparsity.
    """
    bb_a = sorted(sublibrary["heldout"] + sublibrary["random"])
    selected = subset_library(table, building_blocks, bb_a)
    unmapped = selected.select(pl.any_horizontal(pl.col(BB_COLUMNS).is_null()).sum()).collect(engine="streaming").item()
    selected.drop_nulls(BB_COLUMNS).sink_parquet(output, statistics=True, maintain_order=True)
    subset = (
        pl.scan_parquet(output)
        .select(
            pl.len().alias("n_compounds"),
            *(pl.col(name).sum().alias(f"depth_{name}") for name in COUNTS),
            *((pl.col(name) == 0).mean().alias(f"zero_{name}") for name in COUNTS),
        )
        .collect(engine="streaming")
        .row(0, named=True)
    )
    diversity = int(sublibrary["diversity"])
    library = library_depth(table)
    library["diversity"] = diversity
    library["saturation"] = {
        name: depth_qc(histogram, diversity) for name, histogram in count_histograms(table).items()
    }
    return {
        "library": library,
        "sublibrary": {
            "n_bb_a": len(bb_a),
            "n_compounds": int(subset["n_compounds"]),
            "n_unmapped": int(unmapped),
            "depth": {name: int(subset[f"depth_{name}"] or 0) for name in COUNTS},
            "fraction_zero": {name: round(float(subset[f"zero_{name}"] or 0.0), 6) for name in COUNTS},
        },
    }
