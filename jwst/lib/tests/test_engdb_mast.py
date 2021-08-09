"""Test the MAST Engineering interface"""
import pytest
import requests

from jwst.lib import engdb_mast


def is_alive(url):
    """Check if a url is alive

    Parameters
    ----------
    url: str
        The URL to check.

    Returns
    -------
    is_alive: bool
        True if alive
    """
    is_alive = False
    try:
        r = requests.get(url)
        is_alive = (r.status_code == requests.codes.ok)
    except Exception:
        pass
    return is_alive


@pytest.mark.skipif(not is_alive(engdb_mast.MAST_BASE_URL),
                    reason=f'MAST url {engdb_mast.MAST_BASE_URL} not available. Skipping.')
def test_aliveness():
    """Check connection creation

    Failure is any failure from instantiation.
    """
    engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL, token='dummytoken')


def test_negative_aliveness():
    """Ensure failure occurs with a bad url"""
    with pytest.raises(RuntimeError):
        engdb_mast.EngdbMast(base_url='https://127.0.0.1/_engdb_mast_test', token='dummytoken')


def test_notoken():
    """Check that failure occurs without a token"""
    with pytest.raises(RuntimeError):
        engdb_mast.EngdbMast()
