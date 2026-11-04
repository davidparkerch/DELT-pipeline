import pytest

from delt_pipeline.library import FRACTIONS, K_MAX, depth_qc, rarefaction


def test_rarefaction_is_hand_computable() -> None:
    # Two compounds with one count, one with two, in a library of four.
    histogram = {1: 2, 2: 1}
    assert rarefaction(histogram, 4, (0.5, 1.0)) == [0.4375, 0.75]  # (2 * 0.5 + 0.75) / 4, and 3 / 4
    curve = rarefaction(histogram, 4)
    assert len(curve) == len(FRACTIONS)
    assert curve == sorted(curve)  # more counts never observe fewer compounds
    assert rarefaction({}, 4, (0.5,)) == [0.0]


def test_rarefaction_rejects_bad_input() -> None:
    with pytest.raises(ValueError, match="at least 1"):
        rarefaction({1: 1}, 0)
    for bad in (0.0, 1.5):
        with pytest.raises(ValueError, match="fractions"):
            rarefaction({1: 1}, 4, (bad,))


def test_depth_qc() -> None:
    qc = depth_qc({1: 2, 2: 1}, 4, k_max=3)
    assert qc["total"] == 4
    assert qc["oversampling"] == 1.0
    assert qc["observed"] == 3
    assert qc["fraction_observed"] == 0.75
    assert qc["fraction_observed_uniform"] == 0.632121  # 1 - exp(-1)
    assert qc["singleton_fraction"] == 0.666667
    assert qc["max_count"] == 2
    assert qc["histogram"] == [1, 2, 1, 0]  # none, 1, 2, >= 3 counts
    assert qc["rarefaction"]["observed"][-1] == 0.75


def test_depth_qc_buckets_the_tail() -> None:
    qc = depth_qc({1: 5, 9: 2, 10: 1, 50: 1}, 100)
    assert len(qc["histogram"]) == K_MAX + 1
    assert qc["histogram"][0] == 100 - 9
    assert qc["histogram"][9] == 2
    assert qc["histogram"][K_MAX] == 2  # 10 and 50 share the last bucket
    assert sum(qc["histogram"]) == 100
    assert qc["max_count"] == 50


def test_depth_qc_of_an_empty_column() -> None:
    qc = depth_qc({}, 10)
    assert qc["total"] == qc["observed"] == qc["max_count"] == 0
    assert qc["singleton_fraction"] is None
    assert qc["rarefaction"]["observed"] == [0.0] * len(FRACTIONS)


def test_depth_qc_rejects_bad_input() -> None:
    with pytest.raises(ValueError, match="k_max"):
        depth_qc({1: 1}, 4, k_max=1)
    with pytest.raises(ValueError, match="exceed the library size"):
        depth_qc({1: 5}, 4)
