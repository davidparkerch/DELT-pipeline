"""A resynthesis shortlist: the best hits, one family at a time.

Hit validation (resynthesis, then a binding assay) is the expensive step of a DEL campaign, and 2,000 ranked
hits are not a decision. Ranking alone also buys redundancy: the best compounds of a selection are usually
close analogues, because an enriched disynthon carries all its members up with it. The shortlist therefore
clusters the trusted (tier 1) hits by structure and picks across the clusters, best first. A cap per
disynthon keeps one family from filling the plate.

Clustering is Butina on Morgan fingerprints (radius 2, 2048 bits), the same family of methods KinDEL used
to choose its own validation compounds.
"""

from typing import Any

import numpy as np
import polars as pl
from rdkit import DataStructs, RDLogger
from rdkit.Chem import MolFromSmiles, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

from .._numbers import significant

RDLogger.DisableLog("rdApp.*")

N_PICKS = 96
"""Compounds on the shortlist: one 96-well plate."""
SIMILARITY = 0.65
"""Tanimoto similarity at or above which two hits fall in one cluster."""
MAX_PER_CLUSTER = 2
MAX_PER_DISYNTHON = 3
ELIGIBLE_TIERS = ("1",)
"""Only trusted hits: replicate and disynthon support, and no counter-screen or bead binder."""
COLUMNS = (
    "screen",
    "pick_rank",
    "pick_round",
    "cluster_id",
    "cluster_size",
    "hit_id",
    "rank",
    "tier",
    "molecule_hash",
    "smiles",
    "bb_a",
    "bb_b",
    "bb_c",
    "disynthon_key",
    "target_sum",
    "control_sum",
    "poisson_lower",
    "eb_log2_ratio",
    "reasons",
)
"""Columns of the shortlist, in order."""


def fingerprints(smiles: list[str]) -> list[DataStructs.ExplicitBitVect]:
    """Morgan fingerprints (radius 2, 2048 bits) of SMILES strings.

    Raises:
        ValueError: If a SMILES cannot be parsed.
    """
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    result = []
    for text in smiles:
        molecule = MolFromSmiles(text)
        if molecule is None:
            msg = f"cannot parse SMILES {text!r}"
            raise ValueError(msg)
        result.append(generator.GetFingerprint(molecule))
    return result


def cluster(prints: list[DataStructs.ExplicitBitVect], similarity: float = SIMILARITY) -> list[int]:
    """Butina cluster label of every fingerprint, numbered by the first (best-ranked) member.

    Args:
        prints: Fingerprints in rank order (best first); ties in the Butina algorithm resolve by position.
        similarity: Tanimoto similarity at or above which a compound joins a cluster centre.

    Returns:
        One label per fingerprint; label 0 is the cluster of the best-ranked compound.
    """
    if not 0 < similarity < 1:
        msg = f"similarity must be in (0, 1), got {similarity}"
        raise ValueError(msg)
    n = len(prints)
    distances: list[float] = []
    for index in range(1, n):
        distances.extend(1.0 - np.array(DataStructs.BulkTanimotoSimilarity(prints[index], prints[:index])))
    groups = Butina.ClusterData(distances, n, 1.0 - similarity, isDistData=True) if n > 1 else [(0,)] * n
    labels = np.zeros(n, dtype=int)
    for members in groups:
        labels[list(members)] = min(members)
    order = {first: label for label, first in enumerate(sorted(set(labels.tolist())))}
    return [order[first] for first in labels.tolist()]


def mean_similarity(prints: list[DataStructs.ExplicitBitVect]) -> float | None:
    """Mean pairwise Tanimoto similarity (lower is more diverse); ``None`` for fewer than two compounds."""
    if len(prints) < 2:  # noqa: PLR2004 - a pair needs two
        return None
    total = sum(float(np.sum(DataStructs.BulkTanimotoSimilarity(prints[i], prints[:i]))) for i in range(1, len(prints)))
    return total / (len(prints) * (len(prints) - 1) / 2)


def _disynthon_key() -> pl.Expr:
    """``AB:12-7``-style identity of the compound's best disynthon (null when it has none)."""
    pairs = {"AB": ("bb_a", "bb_b"), "AC": ("bb_a", "bb_c"), "BC": ("bb_b", "bb_c")}
    expression = pl.lit(None, dtype=pl.String)
    for name, (first, second) in pairs.items():
        expression = (
            pl.when(pl.col("disynthon") == name)
            .then(pl.format("{}:{}-{}", pl.lit(name), pl.col(first), pl.col(second)))
            .otherwise(expression)
        )
    return expression.alias("disynthon_key")


def select_picks(
    clusters: list[int],
    keys: list[str | None],
    n_picks: int = N_PICKS,
    max_per_cluster: int = MAX_PER_CLUSTER,
    max_per_disynthon: int = MAX_PER_DISYNTHON,
) -> list[tuple[int, int]]:
    """Choose compounds across clusters, best first.

    Round 1 takes the best compound of every cluster (clusters ordered by their best rank), round 2 the
    next best of every cluster, and so on up to ``max_per_cluster``. A compound is skipped when its
    disynthon already has ``max_per_disynthon`` picks.

    Args:
        clusters: Cluster label per eligible compound, in rank order.
        keys: Disynthon key per compound (``None`` has no cap).
        n_picks: Compounds to pick at most.
        max_per_cluster: Picks per cluster at most.
        max_per_disynthon: Picks per disynthon at most.

    Returns:
        ``(index, round)`` pairs in pick order.
    """
    if min(n_picks, max_per_cluster, max_per_disynthon) < 1:
        msg = "n_picks, max_per_cluster and max_per_disynthon must be at least 1"
        raise ValueError(msg)
    members: dict[int, list[int]] = {}
    for index, label in enumerate(clusters):
        members.setdefault(label, []).append(index)
    taken: set[int] = set()
    per_disynthon: dict[str, int] = {}
    picks: list[tuple[int, int]] = []
    for pick_round in range(1, max_per_cluster + 1):
        before = len(picks)
        for label in sorted(members, key=lambda key: members[key][0]):
            if len(picks) >= n_picks:
                return picks
            for index in members[label]:
                key = keys[index]
                if index in taken or (key is not None and per_disynthon.get(key, 0) >= max_per_disynthon):
                    continue
                taken.add(index)
                if key is not None:
                    per_disynthon[key] = per_disynthon.get(key, 0) + 1
                picks.append((index, pick_round))
                break
        if len(picks) == before:
            break
    return picks


def shortlist(
    hits: pl.DataFrame,
    *,
    n_picks: int = N_PICKS,
    similarity: float = SIMILARITY,
    max_per_cluster: int = MAX_PER_CLUSTER,
    max_per_disynthon: int = MAX_PER_DISYNTHON,
    tiers: tuple[str, ...] = ELIGIBLE_TIERS,
) -> tuple[pl.DataFrame, dict[str, Any]]:
    """Pick a diverse resynthesis shortlist for every screen of the triaged hits.

    Args:
        hits: The triage table (``triage.triage``).
        n_picks: Compounds per screen at most.
        similarity: Tanimoto similarity at or above which two hits share a cluster.
        max_per_cluster: Picks per cluster at most.
        max_per_disynthon: Picks per disynthon at most.
        tiers: Triage tiers that are eligible.

    Returns:
        ``(picks, qc)``: the shortlist (``COLUMNS``, sorted by screen and ``pick_rank``) and, per screen, how
        many clusters and disynthons it spans next to the naive top ``n_picks`` by rank, and the mean
        pairwise similarity of both.
    """
    frames = []
    qc: dict[str, Any] = {
        "parameters": {
            "n_picks": n_picks,
            "similarity": similarity,
            "max_per_cluster": max_per_cluster,
            "max_per_disynthon": max_per_disynthon,
            "tiers": list(tiers),
        },
        "screens": {},
    }
    for screen in sorted(hits.get_column("screen").unique().to_list()):
        eligible = (
            hits.filter((pl.col("screen") == screen) & pl.col("tier").is_in(tiers))
            .sort("rank")
            .with_columns(_disynthon_key())
        )
        prints = fingerprints(eligible.get_column("smiles").to_list())
        ranks = eligible.get_column("rank").to_list()
        labels = cluster(prints, similarity) if prints else []
        keys = eligible.get_column("disynthon_key").to_list()
        chosen = select_picks(labels, keys, n_picks, max_per_cluster, max_per_disynthon)
        sizes = np.bincount(labels) if labels else np.array([], dtype=int)
        frame = eligible[[index for index, _ in chosen], :].with_columns(
            pick_rank=pl.Series(range(1, len(chosen) + 1), dtype=pl.Int32),
            pick_round=pl.Series([pick_round for _, pick_round in chosen], dtype=pl.Int8),
            cluster_id=pl.Series([labels[index] for index, _ in chosen], dtype=pl.Int32),
            cluster_size=pl.Series([int(sizes[labels[index]]) for index, _ in chosen], dtype=pl.Int32),
        )
        frames.append(frame.select(COLUMNS))
        naive = list(range(min(n_picks, len(labels))))
        picked = [index for index, _ in chosen]
        qc["screens"][screen] = {
            "n_eligible": eligible.height,
            "n_clusters": int(sizes.size),
            "n_picks": len(chosen),
            "picks": _coverage(labels, keys, picked, prints, ranks),
            "naive_top": _coverage(labels, keys, naive, prints, ranks),
        }
    picks = pl.concat(frames) if frames else pl.DataFrame(schema=dict.fromkeys(COLUMNS, pl.String))
    return picks, qc


def _coverage(
    labels: list[int],
    keys: list[str | None],
    indices: list[int],
    prints: list[DataStructs.ExplicitBitVect],
    ranks: list[int],
) -> dict[str, Any]:
    """What a set of compounds covers: clusters, disynthons, how alike they are, and how deep in the ranking."""
    similarity = mean_similarity([prints[index] for index in indices])
    return {
        "clusters": len({labels[index] for index in indices}),
        "disynthons": len({keys[index] for index in indices if keys[index] is not None}),
        "mean_similarity": significant(similarity) if similarity is not None else None,
        "max_rank": max((ranks[index] for index in indices), default=None),
    }
