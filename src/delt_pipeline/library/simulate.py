"""Synthetic DEL screens in the KinDEL file layout, with planted signals for tests and the quick start.

The library is small (15 x 25 x 16 = 6,000 compounds) but has the features that make DEL data hard:

* almost every count is zero, and compounds with no counts at all are absent from the table;
* uneven pre-selection abundance (``seq_load``);
* a **binder family**: every compound sharing one disynthon (BB-A x BB-B) binds the target. It is
  well synthesised, so it also has many control counts (non-specific background);
* a **low-count outlier**: 6 target counts in one replicate, 0 control counts, which naive enrichment ratios love;
* a **bead binder**: binds the matrix, so it is enriched in every target *and* in the control;
* a **promiscuous compound**: enriched in both kinases and in the BCA counter-screen;
* a **kinase-shared family**: binds both kinases (MAPK14, DDR1) but not BCA, like a hinge binder.

SMILES are simple valid fragments joined with ``.``: chemically meaningless, but parseable.
"""

import hashlib
from dataclasses import dataclass, field

import numpy as np
import polars as pl

from .schema import KINDEL_CONTROL, KINDEL_LOAD, KINDEL_TARGET, SMILES_COLUMNS

N_A, N_B, N_C = 15, 25, 16

BACKGROUND_RATE = 0.16
"""Expected sequence count per replicate for an average compound with no binding."""

NOISE = 0.7
"""Standard deviation of the log-normal compound-level noise in each arm."""

HELDOUT_BACKGROUND_A = (0, 2, 7)
"""BB-As that the held-out non-binders are drawn from."""


@dataclass(frozen=True)
class Planted:
    """Where the planted signals are, as building-block indices (0-based)."""

    binder: tuple[int, int] = (2, 5)
    """Disynthon (A, B) that binds MAPK14."""
    ddr1_binder: tuple[int, int] = (4, 17)
    """Disynthon (A, B) that binds DDR1."""
    bca_binder: tuple[int, int] = (9, 20)
    """Disynthon (A, B) that binds BCA."""
    shared: tuple[int, int] = (7, 11)
    """Disynthon (A, B) that binds both kinases."""
    outlier: tuple[int, int, int] = (1, 1, 1)
    """Compound with 6 target counts (one replicate) and no control counts in MAPK14."""
    bead: tuple[int, int, int] = (11, 3, 6)
    """Compound that binds the bead matrix."""
    promiscuous: tuple[int, int, int] = (12, 8, 2)
    """Compound enriched in all three proteins."""


@dataclass(frozen=True)
class SimulatedLibrary:
    """Building blocks, one KinDEL-layout table per screen and held-out KD tables."""

    building_blocks: dict[str, pl.DataFrame]
    screens: dict[str, pl.DataFrame]
    heldout: dict[str, pl.DataFrame]
    planted: Planted = field(default_factory=Planted)


def synthon_smiles() -> dict[str, list[str]]:
    """Distinct, valid SMILES for each position: amino acids (A), benzoic acids (B), pyridyl amines (C)."""
    return {
        "a": [f"NC{'C' * index}C(=O)O" for index in range(N_A)],
        "b": [f"OC(=O)c1ccc({'C' * (index + 1)})cc1" for index in range(N_B)],
        "c": [f"NC{'C' * index}c1ccncc1" for index in range(N_C)],
    }


def molecule_hash(smiles: str) -> str:
    """SHA-256 of the SMILES, the same kind of identifier as KinDEL's ``molecule_hash``."""
    return hashlib.sha256(smiles.encode()).hexdigest()


def _grid() -> pl.DataFrame:
    a, b, c = np.meshgrid(np.arange(N_A), np.arange(N_B), np.arange(N_C), indexing="ij")
    return pl.DataFrame({"a": a.ravel(), "b": b.ravel(), "c": c.ravel()}, schema=dict.fromkeys("abc", pl.Int64))


def _disynthon(grid: pl.DataFrame, ab: tuple[int, int]) -> np.ndarray:
    return ((grid["a"] == ab[0]) & (grid["b"] == ab[1])).to_numpy()


def _compound(grid: pl.DataFrame, abc: tuple[int, int, int]) -> np.ndarray:
    return ((grid["a"] == abc[0]) & (grid["b"] == abc[1]) & (grid["c"] == abc[2])).to_numpy()


def _enrichment(grid: pl.DataFrame, planted: Planted) -> dict[str, np.ndarray]:
    """True fold enrichment of each compound in each protein selection (1 = no binding)."""
    shared = _disynthon(grid, planted.shared)
    promiscuous = _compound(grid, planted.promiscuous)
    folds = {}
    for screen, binder, fold in (
        ("mapk14", planted.binder, 4.0),
        ("ddr1", planted.ddr1_binder, 4.0),
        ("bca", planted.bca_binder, 4.0),
    ):
        enrichment = np.ones(grid.height)
        enrichment[_disynthon(grid, binder)] = fold
        if screen != "bca":
            enrichment[shared] = 4.0
        enrichment[promiscuous] = 20.0
        folds[screen] = enrichment
    return folds


def _heldout(grid: pl.DataFrame, products: list[str], rates: np.ndarray, rng: np.random.Generator) -> pl.DataFrame:
    """KD (nM) for some in-library compounds, decreasing with their true enrichment, plus external compounds."""
    families = np.flatnonzero(rates > 1.0)
    # Non-binders are taken from a few BB-As only, so most of the library stays unbiased (random BB-As).
    candidates = (rates == 1.0) & np.isin(grid["a"].to_numpy(), HELDOUT_BACKGROUND_A)
    background = rng.choice(np.flatnonzero(candidates), size=6, replace=False)
    chosen = np.concatenate([families[:: max(1, families.size // 8)], background])
    kd = 20_000.0 / rates[chosen] * rng.lognormal(0.0, 0.35, chosen.size)
    smiles = synthon_smiles()
    rows = {
        "smiles": [products[index] for index in chosen] + ["CCOc1ccccc1", "CCN(CC)C(=O)c1ccccc1"],
        "molecule_hash": [molecule_hash(products[index]) for index in chosen] + [None, None],
        "kd": [round(float(value), 1) for value in kd] + [150.0, 9000.0],
    }
    for position in "abc":
        indices = grid[position].to_numpy()[chosen]
        rows[f"smiles_{position}"] = [smiles[position][index] for index in indices] + [None, None]
    frame = pl.DataFrame(
        rows, schema_overrides={"molecule_hash": pl.String, **dict.fromkeys(SMILES_COLUMNS, pl.String)}
    )
    return frame.with_row_index("").select("", "smiles", "molecule_hash", "kd", *SMILES_COLUMNS)


def simulate_library(seed: int = 7) -> SimulatedLibrary:
    """Simulate three screens (MAPK14, DDR1, BCA) of one library, plus held-out KD for the kinases.

    Args:
        seed: Random seed; the result is fully determined by it.

    Returns:
        The simulated building blocks, screens and held-out tables.
    """
    rng = np.random.default_rng(seed)
    planted = Planted()
    grid = _grid()
    smiles = synthon_smiles()
    products = [
        f"{smiles['a'][a]}.{smiles['b'][b]}.{smiles['c'][c]}"
        for a, b, c in zip(grid["a"].to_list(), grid["b"].to_list(), grid["c"].to_list(), strict=True)
    ]
    # Pre-selection abundance: uneven synthesis yields; the planted families are well made.
    abundance = rng.gamma(shape=2.0, scale=0.5, size=grid.height)
    abundance[_disynthon(grid, planted.binder)] *= 12.0
    for family in (planted.ddr1_binder, planted.bca_binder, planted.shared):
        abundance[_disynthon(grid, family)] *= 6.0
    abundance[_compound(grid, planted.promiscuous)] *= 5.0
    load = rng.poisson(abundance * 6.0)

    folds = _enrichment(grid, planted)
    bead = _compound(grid, planted.bead)
    outlier = _compound(grid, planted.outlier)
    screens = {}
    for screen, enrichment in folds.items():
        # Compound-level noise in each arm (synthesis, PCR): the overdispersion a DEL prior has to learn.
        control_rate = abundance * BACKGROUND_RATE * rng.lognormal(0.0, NOISE, grid.height)
        control_rate[bead] *= 60.0
        target_rate = control_rate * enrichment * rng.lognormal(0.0, NOISE, grid.height)
        target = rng.poisson(target_rate[:, None], size=(grid.height, 3))
        control = rng.poisson(control_rate[:, None], size=(grid.height, 3))
        if screen == "mapk14":
            target[outlier] = (0, 6, 0)
            control[outlier] = (0, 0, 0)
        columns: dict[str, object] = {
            "smiles": products,
            "molecule_hash": [molecule_hash(product) for product in products],
            **{f"smiles_{position}": [smiles[position][index] for index in grid[position]] for position in "abc"},
            **{name: target[:, rep] for rep, name in enumerate(KINDEL_TARGET)},
            **{name: control[:, rep] for rep, name in enumerate(KINDEL_CONTROL)},
            KINDEL_LOAD: load,
        }
        frame = pl.DataFrame(columns).with_columns(target_sum=pl.sum_horizontal(KINDEL_TARGET))
        # Like KinDEL, compounds that were never sequenced in any sample are absent.
        observed = frame.filter(pl.sum_horizontal(*KINDEL_TARGET, *KINDEL_CONTROL, KINDEL_LOAD) > 0)
        screens[screen] = observed

    heldout = {}
    for screen in ("mapk14", "ddr1"):
        for condition, noise_seed in (("ondna", 1), ("offdna", 2)):
            sub_rng = np.random.default_rng(seed * 100 + noise_seed + (0 if screen == "mapk14" else 10))
            heldout[f"{screen}_{condition}"] = _heldout(grid, products, folds[screen], sub_rng)

    building_blocks = {
        position: pl.DataFrame({f"smiles_{position}": names, f"reacted_{position}": names})
        for position, names in smiles.items()
    }
    return SimulatedLibrary(building_blocks, screens, heldout, planted)
