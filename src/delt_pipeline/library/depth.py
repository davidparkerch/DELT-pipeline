"""Sequencing depth and library saturation: how much of the library the counts actually sampled.

A DEL selection is only informative for a compound that was sampled. The protocol literature asks for
sequence counts that exceed the library size, ideally tenfold (Decurtins et al., Nat. Protoc. 2016).
This module turns a count histogram (compounds per exact count) into the numbers a reviewer needs:
counts per compound, the fraction of the library seen at all, and a rarefaction curve, the same
library-complexity curve RNA-seq pipelines draw.
"""

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from .._numbers import significant

K_MAX = 10
MIN_K_MAX = 2
"""Counts at or above this value share the last histogram bucket in the QC file."""
FRACTIONS = (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
"""Fractions of the observed counts at which the rarefaction curve is evaluated."""
RECOMMENDED_OVERSAMPLING = 10.0
"""Counts per library member recommended for statistically sound selections (Decurtins 2016)."""


def rarefaction(histogram: Mapping[int, int], diversity: int, fractions: tuple[float, ...] = FRACTIONS) -> list[float]:
    """Expected fraction of the library observed after keeping only a fraction of the counts.

    Keeping each count independently with probability ``f`` leaves a compound with ``k`` counts
    unobserved with probability ``(1 - f)**k``, so the expected number of observed compounds is
    ``sum_k n_k * (1 - (1 - f)**k)``. This is exact for binomial thinning, needs no resampling and is
    deterministic.

    Args:
        histogram: Number of compounds ``n_k`` with exactly ``k >= 1`` counts.
        diversity: Library size (the number of possible compounds).
        fractions: Fractions ``f`` in ``(0, 1]``.

    Returns:
        The expected observed fraction of the library for every ``f``.
    """
    if diversity < 1:
        msg = f"diversity must be at least 1, got {diversity}"
        raise ValueError(msg)
    if not all(0 < f <= 1 for f in fractions):
        msg = f"fractions must lie in (0, 1], got {fractions}"
        raise ValueError(msg)
    counts = np.array(sorted(histogram), dtype=float)
    members = np.array([histogram[int(k)] for k in counts], dtype=float)
    curve = []
    for f in fractions:
        seen = np.ones_like(counts) if f == 1 else -np.expm1(counts * np.log1p(-f))
        curve.append(float(significant(float(np.sum(members * seen)) / diversity)))
    return curve


def depth_qc(histogram: Mapping[int, int], diversity: int, k_max: int = K_MAX) -> dict[str, Any]:
    """Depth and saturation summary of one count column.

    Args:
        histogram: Number of compounds with exactly ``k >= 1`` counts (compounds without counts are
            implied by ``diversity``).
        diversity: Library size.
        k_max: Counts of at least this value share the last bucket of the stored histogram.

    Returns:
        ``total`` counts, ``oversampling`` (counts per library member), ``observed`` compounds and
        their fraction, the fraction a perfectly uniform library would show at this depth
        (``1 - exp(-oversampling)``), the share of observed compounds seen exactly once, the largest
        count, the bucketed ``histogram`` (index ``k`` is the number of compounds with exactly ``k``
        counts, the last bucket is ``>= k_max``) and the ``rarefaction`` curve.
    """
    if k_max < MIN_K_MAX:
        msg = f"k_max must be at least {MIN_K_MAX}, got {k_max}"
        raise ValueError(msg)
    observed = sum(histogram.values())
    if observed > diversity:
        msg = f"{observed} observed compounds exceed the library size {diversity}"
        raise ValueError(msg)
    total = sum(k * n for k, n in histogram.items())
    oversampling = total / diversity
    buckets = [diversity - observed] + [histogram.get(k, 0) for k in range(1, k_max)]
    buckets.append(sum(n for k, n in histogram.items() if k >= k_max))
    return {
        "total": int(total),
        "oversampling": significant(oversampling),
        "observed": int(observed),
        "fraction_observed": significant(observed / diversity),
        "fraction_observed_uniform": significant(-math.expm1(-oversampling)),
        "singleton_fraction": significant(histogram.get(1, 0) / observed) if observed else None,
        "max_count": max(histogram, default=0),
        "histogram": buckets,
        "rarefaction": {"fraction": list(FRACTIONS), "observed": rarefaction(histogram, diversity)},
    }
