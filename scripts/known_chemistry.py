"""Check building-block set enrichment against known chemistry: carbonic anhydrase binds primary sulfonamides.

KinDEL's BCA screen is a carbonic anhydrase selection. Primary sulfonamides are the classic zinc-binding group of
carbonic anhydrase inhibitors, so a method that finds real binders should rank the sulfonamide building blocks
high in BCA without being told. This script finds the primary-sulfonamide blocks with RDKit, looks up their
set-enrichment rank in every screen and computes the chance of a ranking this good (all sulfonamide blocks within
the best ``k`` of the position's blocks).

Usage: ``python scripts/known_chemistry.py RESULTS_DIR`` (or ``just known-chemistry``).
"""

import sys
from math import comb
from pathlib import Path

import polars as pl
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
PRIMARY_SULFONAMIDE = Chem.MolFromSmarts("[NX3;H2]S(=O)(=O)")


def sulfonamide_blocks(blocks: pl.DataFrame) -> pl.DataFrame:
    """The building blocks whose fragment contains a primary sulfonamide."""
    has = [Chem.MolFromSmiles(smiles).HasSubstructMatch(PRIMARY_SULFONAMIDE) for smiles in blocks["canonical_smiles"]]
    return blocks.filter(pl.Series(has))


def chance_all_in_top(n_blocks: int, n_hits: int, k: int) -> float:
    """Probability that ``n_hits`` blocks picked at random from ``n_blocks`` all fall within the best ``k``."""
    return comb(k, n_hits) / comb(n_blocks, n_hits)


def main(results: Path) -> None:
    """Print the rank of every primary-sulfonamide block in every screen."""
    blocks = pl.read_parquet(results / "blocks/blocks.parquet")
    found = sulfonamide_blocks(blocks)
    print(f"{found.height} primary-sulfonamide building blocks: {', '.join(found['block_id'])}")
    for path in sorted(results.glob("screens/*/*.blockset.parquet")):
        screen = path.name.removesuffix(".blockset.parquet")
        table = pl.read_parquet(path).sort("auc", descending=True)
        for position in sorted(found["position"].unique()):
            ranked = table.filter(pl.col("position") == position).with_row_index("rank", offset=1)
            joined = ranked.join(found.select("position", "bb_id", "block_id"), on=["position", "bb_id"])
            ranks = sorted(joined["rank"].to_list())
            print(f"{screen:8s} position {position}: ranks {ranks} of {ranked.height} blocks")
            if screen == "bca":
                k = max(ranks)
                chance = chance_all_in_top(ranked.height, len(ranks), k)
                print(f"         chance that all {len(ranks)} fall within the top {k} at random: {chance:.1e}")
                print(joined.select("block_id", "rank", "auc", "n_top").sort("rank"))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
