"""Tests for set_telescope_pointing that does not need DB connection."""

import pytest

pytest.importorskip("pysiaf")

from astropy.time import Time  # noqa: E402

from jwst.lib import set_telescope_pointing as stp  # noqa: E402


@pytest.mark.parametrize("method", [method for method in stp.Methods])
def test_method_string(method):
    """Ensure that the value of the method is the string representation"""
    assert f"{method}" == method.value


def test_change_engdb_url_fail():
    """Test changing the engineering database by call"""
    with pytest.raises(Exception):
        stp.get_pointing(
            Time("2019-06-03T17:25:40", format="isot").mjd,
            Time("2019-06-03T17:25:56", format="isot").mjd,
            engdb_url="http://localhost",
        )
