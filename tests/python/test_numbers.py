import math

from delt_pipeline._numbers import significant


def test_significant() -> None:
    assert significant(1234.5678) == 1234.57
    assert significant(0.000123456789, 3) == 0.000123
    assert significant(math.nan) is None
    assert significant(math.inf) is None
