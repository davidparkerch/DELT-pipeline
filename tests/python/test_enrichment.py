import numpy as np
import pytest
from scipy import stats

from delt_pipeline.analyse.enrichment import (
    BetaPrior,
    Exposure,
    count_difference,
    enrichment_ratio,
    faver_z,
    fit_beta_prior,
    poisson_ratio,
    shrunk_log2_ratio,
)

EVEN = Exposure(target=1000.0, control=1000.0)


def odds(pi: float) -> float:
    return pi / (1 - pi)


def test_exposure() -> None:
    exposure = Exposure(target=200.0, control=100.0)
    assert exposure.ratio == 0.5
    assert exposure.null_fraction == pytest.approx(2 / 3)
    with pytest.raises(ValueError, match="positive"):
        Exposure(target=0.0, control=1.0)


def test_count_difference() -> None:
    np.testing.assert_allclose(count_difference(np.array([3.0, 1.0]), np.array([1.0, 1.0])), [2.0, 0.0])


def test_enrichment_ratio_by_hand() -> None:
    # (3 + 0.5) / (0 + 0.5) * (200 / 100)
    ratio = enrichment_ratio(np.array([3]), np.array([0]), Exposure(target=100.0, control=200.0), alpha=0.5)
    assert ratio[0] == pytest.approx(14.0)
    with pytest.raises(ValueError, match="pseudocount"):
        enrichment_ratio(np.array([1]), np.array([1]), EVEN, alpha=0)


def test_poisson_interval_matches_clopper_pearson() -> None:
    estimate, lower, upper = poisson_ratio(np.array([10]), np.array([10]), EVEN)
    reference = stats.binomtest(10, 20).proportion_ci(confidence_level=0.95, method="exact")
    assert estimate[0] == pytest.approx(1.0)
    assert lower[0] == pytest.approx(odds(reference.low))
    assert upper[0] == pytest.approx(odds(reference.high))
    assert lower[0] == pytest.approx(1 / upper[0])  # symmetric when the counts are equal


def test_poisson_interval_uses_the_exposure_not_one_half() -> None:
    """Equal true rates, but twice the target depth: the interval must contain 1.

    Treating pi = 0.5 as "no enrichment" would call this compound enriched (odds lower bound ~ 1.6).
    """
    exposure = Exposure(target=2000.0, control=1000.0)
    estimate, lower, upper = poisson_ratio(np.array([200]), np.array([100]), exposure)
    naive = stats.binomtest(200, 300).proportion_ci(method="exact").low
    assert odds(naive) > 1.0
    assert estimate[0] == pytest.approx(1.0)
    assert lower[0] < 1.0 < upper[0]


def test_poisson_interval_edges() -> None:
    estimate, lower, upper = poisson_ratio(np.array([0, 5, 0]), np.array([4, 0, 0]), EVEN)
    assert estimate[0] == 0.0
    assert lower[0] == 0.0
    assert np.isinf(estimate[1])
    assert np.isinf(upper[1])
    assert lower[1] > 0
    assert np.isnan(estimate[2])
    assert (lower[2], upper[2]) == (0.0, np.inf)
    for level in (0.0, 1.0):
        with pytest.raises(ValueError, match="confidence level"):
            poisson_ratio(np.array([1]), np.array([1]), EVEN, level=level)


def test_lower_bound_prefers_evidence_over_a_lucky_ratio() -> None:
    """The core point: 3 vs 0 has a higher ratio than 200 vs 40, but much less evidence."""
    c_t, c_c = np.array([3, 200]), np.array([0, 40])
    ratio = enrichment_ratio(c_t, c_c, EVEN)
    _, lower, _ = poisson_ratio(c_t, c_c, EVEN)
    assert ratio[0] > ratio[1]
    assert lower[0] < 1.0 < lower[1]


def test_faver_z_by_hand() -> None:
    z = faver_z(np.array([10, 20]), target_total=1000, diversity=100)
    assert z[0] == pytest.approx(0.0)
    assert z[1] == pytest.approx(0.01 / np.sqrt(0.01 * 0.99))
    with pytest.raises(ValueError, match="diversity"):
        faver_z(np.array([1]), target_total=10, diversity=1)


def test_fit_beta_prior_recovers_parameters() -> None:
    rng = np.random.default_rng(3)
    pi = rng.beta(2.0, 3.0, 20_000)
    n = rng.poisson(20, pi.size) + 1
    c_t = rng.binomial(n, pi)
    prior = fit_beta_prior(c_t, n - c_t)
    assert prior.mean == pytest.approx(0.4, abs=0.01)
    assert prior.a == pytest.approx(2.0, rel=0.1)
    assert prior.b == pytest.approx(3.0, rel=0.1)
    assert prior.n_compounds == pi.size


def test_fit_beta_prior_ignores_empty_compounds_and_needs_some() -> None:
    prior = fit_beta_prior(np.array([0, 3, 1]), np.array([0, 1, 2]))
    assert prior.n_compounds == 2
    with pytest.raises(ValueError, match="no compound"):
        fit_beta_prior(np.array([0, 0]), np.array([0, 0]))


def test_shrinkage() -> None:
    prior = BetaPrior(a=4.0, b=4.0, n_compounds=100)
    c_t, c_c = np.array([0, 3, 300]), np.array([0, 0, 0])
    shrunk = shrunk_log2_ratio(c_t, c_c, prior, EVEN)
    assert shrunk[0] == pytest.approx(0.0)  # no data: the prior (centred on no enrichment)
    assert 0 < shrunk[1] < np.log2(enrichment_ratio(c_t, c_c, EVEN)[1])  # 3 vs 0 is pulled in hard
    assert shrunk[2] > 5  # 300 vs 0 keeps a large ratio
    assert prior.mean == 0.5
