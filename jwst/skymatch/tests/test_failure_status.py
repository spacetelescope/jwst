import pytest

from jwst.skymatch.skymatch_step import _sky_match_status


@pytest.mark.parametrize(
    "is_valid, expected",
    [(True, "COMPLETE"), (False, "FAILED")],
)
def test_sky_match_status(is_valid, expected):
    assert _sky_match_status(is_valid) == expected
