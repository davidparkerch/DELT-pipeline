"""Numeric helpers for reproducible, JSON-friendly results."""

import math

SIGNIFICANT_DIGITS = 6


def significant(value: float, digits: int = SIGNIFICANT_DIGITS) -> float | None:
    """Round to ``digits`` significant digits; NaN and infinities become ``None``.

    Rounding keeps outputs byte-stable across CPUs and BLAS builds, which is what
    lets nf-test snapshot them.
    """
    number = float(value)
    if not math.isfinite(number):
        return None
    return float(f"{number:.{digits}g}")
