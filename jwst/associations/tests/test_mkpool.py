"""
Test mkpool
"""
from glob import glob
import os
import pytest

from astropy.io import fits

from jwst.associations.tests import helpers
from jwst.associations import AssociationPool
from jwst.associations.mkpool import mkpool

REQUIRED_PARAMS = set(('program', 'filename'))


@pytest.fixture(scope='module')
def exposures():
    exposure_path = helpers.t_path(
        'data/exposures'
    )
    exposures = glob(os.path.join(exposure_path, '*.fits'))
    return exposures


def test_mkpool(exposures):
    pool = mkpool(exposures)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)
    filenames = [
        filename
        for filename in pool['filename']
    ]
    assert set(exposures) == set(filenames)


def test_hdulist(exposures):
    hduls = [
        fits.open(exposure)
        for exposure in exposures
    ]
    pool = mkpool(hduls)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)


def test_hdu(exposures):
    hdus = [
        fits.open(exposure)[0]
        for exposure in exposures
    ]
    pool = mkpool(hdus)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)
