"""A fingerprint model on DEL counts, judged by the same measured KD as the statistical scores.

The question this answers is practical: does a structure-based model beat a well-specified count statistic?
KinDEL tested that with machine-learning models and reported that Counts and Poisson enrichment rank the
on-DNA validation compounds better than every model, and that the models lose most of that when the test
compounds share no disynthon with the training data. This module repeats the experiment with the simplest
reasonable model, a gradient-boosted Poisson regression of the target counts on Morgan fingerprints, and
reports it next to the paper's table.

Two targets are compared. **Counts** regresses the pooled target counts alone, like a model that ignores the
control. **Control-weighted** regresses the target count per unit of the compound's own control count (a
Poisson regression with the control as exposure), so that synthesis abundance and bead binding are not
mistaken for binding.

Two splits are compared. **Random**: only the validated compounds themselves are held out. **Disynthon**: no
training compound shares a disynthon (AB, AC or BC) with any validated compound, so the model must
generalise to building-block pairs it has not seen.

The training set is a sample, because the sub-library has about 11 M compounds: the compounds with the most
target counts (the part with signal, like KinDEL's top-1M training set) plus a deterministic random sample
(the background). It is biased on purpose, and the prediction is a ranking score, not a calibrated count.
"""

from collections.abc import Iterator
from importlib import resources
from io import StringIO
from typing import Any

import numpy as np
import polars as pl
from sklearn.ensemble import HistGradientBoostingRegressor

from .._chem import FP_BITS, morgan_generator, parse_smiles
from ..library.schema import ID, POSITIONS
from .validate import PUBLISHED, negative_spearman

N_TOP = 60_000
"""Compounds with the most target counts in the training sample."""
N_RANDOM = 60_000
"""Compounds drawn at random from the rest (about; a deterministic hash filter)."""
SPLITS = ("random", "disynthon")
PARAMETERS = {
    "max_iter": 150,
    "learning_rate": 0.1,
    "max_leaf_nodes": 31,
    "min_samples_leaf": 20,
    "l2_regularization": 1.0,
}
"""Settings of the gradient-boosted model (the same for every fit)."""
DISYNTHONS = (("bb_a", "bb_b"), ("bb_a", "bb_c"), ("bb_b", "bb_c"))
COLUMNS = (ID, "smiles", "bb_a", "bb_b", "bb_c", "target_sum", "control_sum")
VARIANTS = ("counts", "control_weighted")
CONTROL_PSEUDOCOUNT = 1.0
"""Added to the control count to make an exposure (a compound with no control counts still has some background)."""


def published_ml(target: str) -> list[dict[str, Any]]:
    """KinDEL's machine-learning results for a target (in-library, -Spearman vs KD, mean and s.d. of 5 runs).

    Source: Chen et al., ICML 2025, Tables 1 and 2, random and disynthon splits. Empty for other targets.
    """
    text = resources.files(__package__).joinpath("data/kindel_ml.csv").read_text(encoding="utf-8")
    return pl.read_csv(StringIO(text)).filter(pl.col("target") == target.lower()).drop("target").to_dicts()


def training_sample(scores: pl.LazyFrame, n_top: int = N_TOP, n_random: int = N_RANDOM, seed: int = 42) -> pl.DataFrame:
    """The compounds the model learns from: the top by target counts plus a random background.

    Args:
        scores: Scored sub-library (``COLUMNS``).
        n_top: Compounds with the most target counts (ties broken by the hash, so the choice is deterministic).
        n_random: About this many more, drawn with a hash filter from the rest of the observed compounds.
        seed: Seed of the hash filter.

    Returns:
        The sample, with ``COLUMNS``.
    """
    base = scores.select(COLUMNS)
    top = base.sort(["target_sum", ID], descending=[True, False]).head(n_top).collect()
    total = int(base.select(pl.len()).collect().item())
    keep = max(1, round(total / max(n_random, 1)))
    background = (
        base.filter(pl.col(ID).hash(seed=seed) % keep == 0)
        .filter(~pl.col(ID).is_in(top.get_column(ID).implode()))
        .collect()
    )
    return pl.concat([top, background]).sort(ID)


def fingerprint_matrix(smiles: list[str]) -> np.ndarray:
    """Morgan fingerprints (radius 2, 2048 bits) as a float32 matrix, one row per SMILES.

    Raises:
        ValueError: If a SMILES cannot be parsed.
    """
    generator = morgan_generator()
    matrix = np.zeros((len(smiles), FP_BITS), dtype=np.float32)
    for row, text in enumerate(smiles):
        matrix[row] = generator.GetFingerprintAsNumPy(parse_smiles(text))
    return matrix


def heldout_blocks(heldout: pl.DataFrame, building_blocks: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """Building-block ids (``bb_a``, ``bb_b``, ``bb_c``) of the in-library held-out compounds.

    Args:
        heldout: A KinDEL ``heldout/*.csv`` table (``molecule_hash``, ``smiles_a``, ``smiles_b``, ``smiles_c``).
        building_blocks: From :func:`~delt_pipeline.library.subset.read_building_blocks`.

    Returns:
        One row per held-out compound that has all three building blocks in the library tables.
    """
    frame = heldout.filter(pl.col(ID).is_not_null()).select(ID, *(f"smiles_{p}" for p in POSITIONS))
    for position in POSITIONS:
        frame = frame.join(building_blocks[position], on=f"smiles_{position}", how="left")
    return frame.drop_nulls([f"bb_{p}" for p in POSITIONS]).select(ID, *(f"bb_{p}" for p in POSITIONS))


def sharing_a_disynthon(sample: pl.DataFrame, held: pl.DataFrame) -> pl.Series:
    """Hashes of the sample compounds that share a disynthon (AB, AC or BC) with a held-out compound."""
    shared = [
        sample.join(held.select(list(pair)).unique().cast(sample.select(list(pair)).schema), on=list(pair), how="semi")
        for pair in DISYNTHONS
    ]
    return pl.concat(shared).get_column(ID).unique() if shared else pl.Series(ID, [], dtype=pl.String)


def _fit(matrix: np.ndarray, train: pl.DataFrame, variant: str, seed: int) -> HistGradientBoostingRegressor:
    """Gradient-boosted Poisson regression of target counts on fingerprint bits.

    ``counts`` fits the counts as they are. ``control_weighted`` fits the rate ``target / (control + 1)`` with
    weight ``control + 1``, which is the Poisson likelihood of the count with the control as exposure.
    """
    counts = train.get_column("target_sum").to_numpy().astype(float)
    model = HistGradientBoostingRegressor(loss="poisson", early_stopping=False, random_state=seed, **PARAMETERS)
    if variant == "counts":
        return model.fit(matrix, counts)
    exposure = train.get_column("control_sum").to_numpy().astype(float) + CONTROL_PSEUDOCOUNT
    return model.fit(matrix, counts / exposure, sample_weight=exposure)


def _splits(sample: pl.DataFrame, held: pl.DataFrame, held_hashes: set[str]) -> Iterator[tuple[str, pl.DataFrame]]:
    """The training set of each split."""
    free = sample.filter(~pl.col(ID).is_in(sorted(held_hashes)))
    yield "random", free
    shared = sharing_a_disynthon(free, held)
    yield "disynthon", free.filter(~pl.col(ID).is_in(shared.implode()))


def model_report(
    scores: pl.LazyFrame,
    building_blocks: dict[str, pl.DataFrame],
    heldout: dict[str, pl.DataFrame],
    target: str,
    *,
    n_top: int = N_TOP,
    n_random: int = N_RANDOM,
    seed: int = 42,
) -> dict[str, Any]:
    """Train the fingerprint model under both splits and correlate its predictions with measured KD.

    Args:
        scores: Scored sub-library.
        building_blocks: From :func:`~delt_pipeline.library.subset.read_building_blocks`.
        heldout: Held-out KD tables by condition (``on``, ``off``); may be empty.
        target: Screen target (``mapk14``, ``ddr1``) to look up the published results.
        n_top: Compounds with the most counts in the training sample.
        n_random: Random background compounds in the training sample (about).
        seed: Seed of the sample and the model.

    Returns:
        ``target``, the ``model`` settings, ``training`` sizes per split, ``results`` (variant, split, condition,
        number of compounds and -Spearman vs KD), KinDEL's ``published`` Counts and Poisson baselines and its
        ``published_ml`` rows. ``results`` is empty without held-out data.
    """
    report: dict[str, Any] = {
        "target": target,
        "model": {
            "estimator": "HistGradientBoostingRegressor",
            "loss": "poisson",
            "fingerprint": "Morgan r2 2048",
            **PARAMETERS,
        },
        "results": [],
        "published": {f"{condition}": PUBLISHED.get((target.lower(), condition), {}) for condition in ("on", "off")},
        "published_ml": published_ml(target),
    }
    held_in_library = {
        condition: table.filter(pl.col(ID).is_not_null() & pl.col("smiles").is_not_null())
        for condition, table in heldout.items()
    }
    if not held_in_library:
        report["training"] = {}
        return report
    sample = training_sample(scores, n_top, n_random, seed)
    held = pl.concat([heldout_blocks(table, building_blocks) for table in held_in_library.values()]).unique()
    held_hashes = {hash_ for table in held_in_library.values() for hash_ in table.get_column(ID).to_list()}
    report["training"] = {}
    for split, train in _splits(sample, held, held_hashes):
        matrix = fingerprint_matrix(train.get_column("smiles").to_list())
        report["training"][split] = {
            "n_train": train.height,
            "n_with_counts": int((train.get_column("target_sum") > 0).sum()),
        }
        for variant in VARIANTS:
            model = _fit(matrix, train, variant, seed)
            for condition, table in held_in_library.items():
                predicted = model.predict(fingerprint_matrix(table.get_column("smiles").to_list()))
                kd = table.get_column("kd").cast(pl.Float64).to_numpy()
                report["results"].append(
                    {
                        "variant": variant,
                        "split": split,
                        "condition": condition,
                        "n_compounds": table.height,
                        "spearman": negative_spearman(predicted, kd),
                    }
                )
    report["training"]["n_excluded_disynthon"] = (
        report["training"]["random"]["n_train"] - report["training"]["disynthon"]["n_train"]
    )
    report["training"]["sample"] = {"n_top": n_top, "n_random": n_random, "n_compounds": sample.height}
    return report
