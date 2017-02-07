"""
Test mkpool
"""
from glob import glob
import os
import pytest

from astropy.io import fits

from . import helpers
from .. import (AssociationRegistry, AssociationPool)
from ..mkpool import mkpool

REQUIRED_PARAMS = set(('PROGRAM', 'FILENAME'))


@pytest.fixture(scope='module')
def env():
    rules = AssociationRegistry()
    exposure_path = helpers.t_path(
        'data/exposures'
    )
    exposures = glob(os.path.join(exposure_path, '*.fits'))
    return rules, exposures


def test_mkpool(env):
    rules, exposures = env
    pool = mkpool(exposures)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)
    filenames = [
        filename
        for filename in pool['FILENAME']
    ]
    assert set(exposures) == set(filenames)


def test_hdulist(env):
    rules, exposures = env
    hduls = [
        fits.open(exposure)
        for exposure in exposures
    ]
    pool = mkpool(hduls)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)


def test_hdu(env):
    rules, exposures = env
    hdus = [
        fits.open(exposure)[0]
        for exposure in exposures
    ]
    pool = mkpool(hdus)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)
