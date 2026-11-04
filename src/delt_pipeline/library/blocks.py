"""The building-block resource: one row per building block, with RDKit descriptors and yield proxies.

KinDEL lists every building block twice: as the free reagent (``smiles_<p>``; position C lists a
reagent pair) and as the fragment that ends up in the compound (``reacted_<p>``), where pseudo-atoms
(francium, astatine) mark the attachment points. The resource describes the **fragment**, capped with
hydrogens, because that is what the compound contains. The pre-selection counts add a synthesis-yield
proxy per block: a block whose compounds are rarely sequenced before selection is a yield or truncation
suspect, and its silence in a selection is weak evidence of non-binding.
"""

from pathlib import Path
from typing import Any

import polars as pl
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from .._numbers import significant
from .schema import BB_COLUMNS, LOAD, POSITIONS

RDLogger.DisableLog("rdApp.*")

ATTACHMENT_ATOMIC_NUMBERS = frozenset({85, 87})
"""Astatine and francium: KinDEL's markers for the attachment points of a reacted building block."""
DIGITS = 5
"""Significant digits of the descriptors (keeps outputs byte-stable)."""
LOW_LOAD_RATIO = 0.5
"""A block is flagged ``low_load`` when its mean pre-selection load is below this share of its position's median."""

DESCRIPTORS = (
    "mw",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "heavy_atoms",
    "fraction_csp3",
    "rings",
)
"""Descriptor columns, in output order."""


def describe_block(reacted: str) -> dict[str, Any]:
    """Describe one reacted building block.

    Args:
        reacted: SMILES with attachment pseudo-atoms (``[Fr]``, ``[At]``).

    Returns:
        ``fragment_smiles`` (canonical, attachment points as ``*``), ``canonical_smiles`` and ``inchikey``
        of the hydrogen-capped fragment, the number of ``attachments``, RDKit descriptors of the capped
        fragment and its Murcko scaffold (``None`` when acyclic).

    Raises:
        ValueError: If RDKit cannot parse the SMILES.
    """
    molecule = Chem.MolFromSmiles(reacted)
    if molecule is None:
        msg = f"cannot parse building block {reacted!r}"
        raise ValueError(msg)
    marked = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() in ATTACHMENT_ATOMIC_NUMBERS]
    fragment, capped = Chem.RWMol(molecule), Chem.RWMol(molecule)
    for index in marked:
        fragment.GetAtomWithIdx(index).SetAtomicNum(0)
        capped.GetAtomWithIdx(index).SetAtomicNum(1)  # replaced by a hydrogen, then removed (keeps stereo)
    capped = Chem.RemoveHs(capped.GetMol())
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=capped)
    values = (
        Descriptors.MolWt(capped),
        Crippen.MolLogP(capped),
        rdMolDescriptors.CalcTPSA(capped),
        rdMolDescriptors.CalcNumHBD(capped),
        rdMolDescriptors.CalcNumHBA(capped),
        rdMolDescriptors.CalcNumRotatableBonds(capped),
        capped.GetNumHeavyAtoms(),
        rdMolDescriptors.CalcFractionCSP3(capped),
        rdMolDescriptors.CalcNumRings(capped),
    )
    return {
        "fragment_smiles": Chem.MolToSmiles(fragment),
        "canonical_smiles": Chem.MolToSmiles(capped),
        "inchikey": Chem.MolToInchiKey(capped),
        "attachments": len(marked),
        **{name: significant(float(value), DIGITS) for name, value in zip(DESCRIPTORS, values, strict=True)},
        "murcko_scaffold": scaffold or None,
    }


def read_reacted(directory: Path) -> pl.DataFrame:
    """Read ``reacted_{a,b,c}.csv`` into one table: ``position``, ``bb_id`` (row number), ``reagent``, ``reacted``.

    Raises:
        ValueError: If a table lacks its ``smiles_<p>`` or ``reacted_<p>`` column.
    """
    frames = []
    for position in POSITIONS:
        path = directory / f"reacted_{position}.csv"
        frame = pl.read_csv(path, infer_schema_length=0)
        for column in (f"smiles_{position}", f"reacted_{position}"):
            if column not in frame.columns:
                msg = f"{path} has no '{column}' column"
                raise ValueError(msg)
        frames.append(
            frame.select(
                pl.lit(position).alias("position"),
                pl.int_range(pl.len(), dtype=pl.UInt16).alias("bb_id"),
                pl.col(f"smiles_{position}").alias("reagent"),
                pl.col(f"reacted_{position}").alias("reacted"),
            )
        )
    return pl.concat(frames)


def describe_blocks(reacted: pl.DataFrame) -> pl.DataFrame:
    """Descriptors for every building block, and which blocks share a fragment.

    Blocks with the same capped fragment (different reagents, for example another protecting group) are
    chemically the same compound part under different DNA codes. ``same_fragment_as`` is the lowest id
    among them and ``n_same_fragment`` how many there are.
    """
    rows = [describe_block(smiles) for smiles in reacted.get_column("reacted")]
    frame = pl.concat([reacted.drop("reacted"), pl.DataFrame(rows, infer_schema_length=None)], how="horizontal")
    return frame.with_columns(
        block_id=pl.col("position").str.to_uppercase() + pl.col("bb_id").cast(pl.String).str.zfill(4),
        same_fragment_as=pl.col("bb_id").min().over("position", "inchikey"),
        n_same_fragment=pl.len().over("position", "inchikey").cast(pl.UInt16),
    )


def load_statistics(counts: pl.DataFrame, sublibrary: dict[str, Any], sizes: dict[str, int]) -> pl.DataFrame:
    """Pre-selection load per building block of the sub-library, a proxy for synthesis yield.

    Args:
        counts: Compact count table (``bb_a``, ``bb_b``, ``bb_c``, ``load``; one row per observed compound).
        sublibrary: The BB-A selection.
        sizes: Number of building blocks per position.

    Returns:
        One row per building block that occurs in the sub-library: ``n_possible`` (compounds it
        could be in), ``n_with_load`` (compounds sequenced before selection), ``mean_load`` (load per possible
        compound) and ``fraction_with_load``. Only the pre-selection load is used: it is one experiment shared by
        all screens, whereas which compounds appear in a screen's table depends on that screen's counts.
    """
    n_a = len(sublibrary["heldout"]) + len(sublibrary["random"])
    sub_size = {"a": sizes["b"] * sizes["c"], "b": n_a * sizes["c"], "c": n_a * sizes["b"]}
    frames = []
    for position in POSITIONS:
        frames.append(
            counts.group_by(f"bb_{position}")
            .agg(
                (pl.col(LOAD) > 0).sum().alias("n_with_load"),
                pl.col(LOAD).sum().alias("load_sum"),
            )
            .select(
                pl.lit(position).alias("position"),
                pl.col(f"bb_{position}").cast(pl.UInt16).alias("bb_id"),
                pl.lit(sub_size[position]).alias("n_possible"),
                "n_with_load",
                pl.col("load_sum"),
            )
        )
    return (
        pl.concat(frames)
        .with_columns(
            mean_load=pl.col("load_sum") / pl.col("n_possible"),
            fraction_with_load=pl.col("n_with_load") / pl.col("n_possible"),
        )
        .drop("load_sum")
    )


def building_block_table(
    reacted: pl.DataFrame, counts: list[Path], sublibrary: dict[str, Any]
) -> tuple[pl.DataFrame, dict[str, Any]]:
    """The searchable building-block resource and a QC summary.

    Args:
        reacted: From :func:`read_reacted`.
        counts: Compact count tables (Parquet), one per screen. The pre-selection load is one experiment
            shared by all screens, so the first table gives the statistics; the others are only checked.
        sublibrary: The BB-A selection.

    Returns:
        ``(table, qc)``: one row per building block, sorted by position and id, and a summary.
    """
    if not counts:
        msg = "need at least one count table"
        raise ValueError(msg)
    sizes = {position: int(reacted.filter(pl.col("position") == position).height) for position in POSITIONS}
    statistics = [
        load_statistics(pl.read_parquet(path, columns=[*BB_COLUMNS, LOAD]), sublibrary, sizes).sort("position", "bb_id")
        for path in sorted(counts)
    ]
    blocks = describe_blocks(reacted).join(statistics[0], on=["position", "bb_id"], how="left")
    source = (
        pl.when(pl.col("position") != "a")
        .then(pl.lit("all"))
        .when(pl.col("bb_id").cast(pl.Int64).is_in(sublibrary["heldout"]))
        .then(pl.lit("heldout"))
        .when(pl.col("bb_id").cast(pl.Int64).is_in(sublibrary["random"]))
        .then(pl.lit("random"))
        .otherwise(pl.lit("excluded"))
    )
    median = pl.col("mean_load").median().over("position")
    blocks = blocks.with_columns(
        sublibrary_source=source,
        load_ratio=pl.col("mean_load") / median,
    ).with_columns(low_load=pl.col("load_ratio") < LOW_LOAD_RATIO)
    columns = [
        "block_id",
        "position",
        "bb_id",
        "reagent",
        "fragment_smiles",
        "canonical_smiles",
        "inchikey",
        "attachments",
        "same_fragment_as",
        "n_same_fragment",
        *DESCRIPTORS,
        "murcko_scaffold",
        "sublibrary_source",
        "n_possible",
        "n_with_load",
        "mean_load",
        "fraction_with_load",
        "load_ratio",
        "low_load",
    ]
    table = blocks.select(columns).sort("position", "bb_id")
    table = table.with_columns(pl.col("mean_load", "fraction_with_load", "load_ratio").round_sig_figs(DIGITS))
    qc = {
        "n_blocks": {position: int(table.filter(pl.col("position") == position).height) for position in POSITIONS},
        "n_in_sublibrary": int(table.filter(pl.col("sublibrary_source") != "excluded").height),
        "n_shared_fragments": int(table.filter(pl.col("n_same_fragment") > 1).height),
        "n_low_load": int(table.filter(pl.col("low_load")).height),
        "low_load_ratio": LOW_LOAD_RATIO,
        "load_shared_across_screens": all(frame.equals(statistics[0]) for frame in statistics[1:]),
    }
    return table, qc
