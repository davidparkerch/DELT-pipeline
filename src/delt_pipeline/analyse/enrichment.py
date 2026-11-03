r"""Enrichment scores for DEL selections: target counts versus no-protein (matrix) control counts.

Every scorer takes pooled counts per compound (``c_t`` in the target, ``c_c`` in the control) and the
*exposure* of each arm: its total sequence count over the **whole** library, so rates stay comparable
with the published KinDEL analysis even when only a sub-library is scored.

=================  =========================================================================
Score              What it is
=================  =========================================================================
``count_diff``     mean target - mean control count (KinDEL's "Counts" baseline)
``er``             enrichment ratio with a pseudocount :math:`\alpha`
``poisson_lower``  lower confidence bound of the Poisson rate ratio (exact, conditional)
``z_faver``        normalised z-score of Faver et al. (2019), target arm only
``eb_log2_ratio``  empirical-Bayes shrunken log2 rate ratio (Beta-Binomial prior)
=================  =========================================================================

The rate ratio is :math:`R = (\lambda_t / E_t) / (\lambda_c / E_c)`. Conditioning on the compound's total
:math:`n = c_t + c_c` turns two Poisson counts into one binomial,
:math:`c_t \mid n \sim \mathrm{Bin}(n, \pi)` with :math:`\pi = R E_t / (R E_t + E_c)`. So an exact
interval for :math:`\pi` (Clopper-Pearson) maps to an exact interval for :math:`R` through
:math:`R = \frac{\pi}{1 - \pi} \cdot \frac{E_c}{E_t}`. "No enrichment" is :math:`\pi_0 = E_t/(E_t+E_c)`,
not 0.5.
"""

from dataclasses import dataclass

import numpy as np
from scipy import optimize, special

MAX_LOG_PRIOR = 12.0
"""Bound on log(a), log(b) of the Beta prior: a nearly-binomial data set otherwise drifts to infinity."""


@dataclass(frozen=True)
class Exposure:
    """Total sequence counts of the target and control arms over the whole library."""

    target: float
    control: float

    def __post_init__(self) -> None:
        if self.target <= 0 or self.control <= 0:
            msg = f"exposures must be positive, got target={self.target}, control={self.control}"
            raise ValueError(msg)

    @property
    def ratio(self) -> float:
        """Control over target exposure; multiplies odds into a rate ratio."""
        return self.control / self.target

    @property
    def null_fraction(self) -> float:
        """Expected fraction of a compound's counts in the target when it is not enriched."""
        return self.target / (self.target + self.control)


def count_difference(target_mean: np.ndarray, control_mean: np.ndarray) -> np.ndarray:
    """KinDEL's "Counts" baseline: mean target count minus mean control count."""
    return np.asarray(target_mean, dtype=float) - np.asarray(control_mean, dtype=float)


def enrichment_ratio(c_t: np.ndarray, c_c: np.ndarray, exposure: Exposure, alpha: float = 0.5) -> np.ndarray:
    r"""Depth-normalised enrichment ratio with a pseudocount.

    .. math:: \mathrm{ER} = \frac{(c_t + \alpha) / E_t}{(c_c + \alpha) / E_c}

    Cheap and intuitive, but a compound with 3 target and 0 control counts can outrank a real binder
    with hundreds of counts in both: the score ignores how much evidence there is.
    """
    if alpha <= 0:
        msg = f"the pseudocount must be positive, got {alpha}"
        raise ValueError(msg)
    c_t = np.asarray(c_t, dtype=float)
    c_c = np.asarray(c_c, dtype=float)
    return (c_t + alpha) / (c_c + alpha) * exposure.ratio


def _odds_to_ratio(pi: np.ndarray, exposure: Exposure) -> np.ndarray:
    with np.errstate(divide="ignore"):
        return pi / (1.0 - pi) * exposure.ratio


def poisson_ratio(
    c_t: np.ndarray, c_c: np.ndarray, exposure: Exposure, level: float = 0.95
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rate ratio estimate with an exact (Clopper-Pearson, conditional binomial) confidence interval.

    Ranking by the **lower bound** rewards compounds with both a high ratio *and* enough counts to be
    sure of it, which is what removes low-count noise.

    Args:
        c_t: Pooled target counts.
        c_c: Pooled control counts.
        exposure: Whole-library totals of both arms.
        level: Two-sided confidence level.

    Returns:
        ``(estimate, lower, upper)``. The estimate is ``inf`` when ``c_c == 0 < c_t`` and ``nan`` when
        both are 0; the lower bound is 0 when ``c_t == 0`` and the upper bound is ``inf`` when ``c_c == 0``.
    """
    if not 0 < level < 1:
        msg = f"the confidence level must be in (0, 1), got {level}"
        raise ValueError(msg)
    c_t = np.asarray(c_t, dtype=float)
    c_c = np.asarray(c_c, dtype=float)
    tail = (1.0 - level) / 2.0
    # Beta quantiles need positive shape parameters; the guarded cases are filled in by np.where.
    pi_lower = np.where(c_t > 0, special.betaincinv(np.maximum(c_t, 1.0), c_c + 1.0, tail), 0.0)
    pi_upper = np.where(c_c > 0, special.betaincinv(c_t + 1.0, np.maximum(c_c, 1.0), 1.0 - tail), 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        estimate = c_t / c_c * exposure.ratio
    return estimate, _odds_to_ratio(pi_lower, exposure), _odds_to_ratio(pi_upper, exposure)


def faver_z(c_t: np.ndarray, target_total: float, diversity: int) -> np.ndarray:
    r"""Normalised z-score of Faver et al. (2019), ACS Comb. Sci. 21, 75-82.

    .. math:: z_n = \frac{p_o - p_e}{\sqrt{p_e (1 - p_e)}}, \quad p_o = c_t / N_t, \quad p_e = 1 / D

    It compares the observed fraction of target reads with a uniform library of diversity :math:`D`, and
    is independent of sequencing depth. It does not use the control, so bead binders score high too.
    """
    if diversity < 2 or target_total <= 0:  # noqa: PLR2004 - a library needs two members
        msg = f"need diversity >= 2 and a positive total, got {diversity} and {target_total}"
        raise ValueError(msg)
    expected = 1.0 / diversity
    observed = np.asarray(c_t, dtype=float) / target_total
    return (observed - expected) / np.sqrt(expected * (1.0 - expected))


@dataclass(frozen=True)
class BetaPrior:
    """Beta(a, b) prior on the target fraction pi, shared by all compounds."""

    a: float
    b: float
    n_compounds: int
    """How many compounds (with at least one count) the prior was fitted on."""

    @property
    def mean(self) -> float:
        """Prior mean of pi."""
        return self.a / (self.a + self.b)


def fit_beta_prior(c_t: np.ndarray, c_c: np.ndarray) -> BetaPrior:
    """Fit a Beta prior to the target fractions by maximising the Beta-Binomial likelihood.

    Only compounds with at least one count carry information. Identical ``(c_t, c_c)`` pairs are
    grouped, so millions of compounds reduce to a few thousand terms.

    Raises:
        ValueError: If no compound has a count.
    """
    c_t = np.asarray(c_t, dtype=np.int64)
    c_c = np.asarray(c_c, dtype=np.int64)
    informative = (c_t + c_c) > 0
    if not informative.any():
        msg = "cannot fit the prior: no compound has any count"
        raise ValueError(msg)
    pairs, weights = np.unique(np.column_stack([c_t[informative], c_c[informative]]), axis=0, return_counts=True)
    k, m = pairs[:, 0].astype(float), pairs[:, 1].astype(float)

    def negative_log_likelihood(log_ab: np.ndarray) -> float:
        a, b = np.exp(log_ab)
        return -float(np.sum(weights * (special.betaln(k + a, m + b) - special.betaln(a, b))))

    # Start from the method of moments on the pooled fraction (a weak prior centred on it).
    start_mean = float(np.clip(k.dot(weights) / (k + m).dot(weights), 0.01, 0.99))
    start = np.log([2.0 * start_mean, 2.0 * (1.0 - start_mean)])
    bounds = [(-MAX_LOG_PRIOR, MAX_LOG_PRIOR)] * 2
    result = optimize.minimize(negative_log_likelihood, start, method="L-BFGS-B", bounds=bounds)
    a, b = np.exp(result.x)
    return BetaPrior(float(a), float(b), int(weights.sum()))


def shrunk_log2_ratio(c_t: np.ndarray, c_c: np.ndarray, prior: BetaPrior, exposure: Exposure) -> np.ndarray:
    r"""Posterior mean of the log2 rate ratio under the Beta prior (empirical Bayes shrinkage).

    The posterior of :math:`\pi` is :math:`\mathrm{Beta}(a + c_t, b + c_c)`, and the posterior mean of the
    log-odds has a closed form, :math:`\psi(a + c_t) - \psi(b + c_c)`. Compounds with few counts are pulled
    toward the library-wide prior; compounds with many counts keep their own ratio. This is the DEL
    analogue of DESeq2's ``lfcShrink``.
    """
    c_t = np.asarray(c_t, dtype=float)
    c_c = np.asarray(c_c, dtype=float)
    log_odds = special.digamma(prior.a + c_t) - special.digamma(prior.b + c_c)
    return (log_odds + np.log(exposure.ratio)) / np.log(2.0)
