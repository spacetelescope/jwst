import os

from asdf import schema as mschema
from numpy.testing import assert_array_equal
import numpy as np
import pytest

from jwst.datamodels import JwstDataModel, RampModel


@pytest.fixture(scope="module")
def fits_file():
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    path = os.path.join(data_dir, 'test.fits')
    return path


def test_from_fits(fits_file):
    with RampModel(fits_file) as dm:
        assert dm.meta.instrument.name == 'MIRI'
        assert dm.shape == (5, 35, 40, 32)


def test_delete():
    with JwstDataModel() as dm:
        dm.meta.instrument.name = 'NIRCAM'
        assert dm.meta.instrument.name == 'NIRCAM'
        del dm.meta.instrument.name
        assert dm.meta.instrument.name is None


def test_fits_without_sci():
    from astropy.io import fits
    schema = {
        "allOf": [
            mschema.load_schema("http://stsci.edu/schemas/jwst_datamodel/core.schema",
                resolve_references=True),
            {
                "type": "object",
                "properties": {
                    "coeffs": {
                        'max_ndim': 1,
                        'fits_hdu': 'COEFFS',
                        'datatype': 'float32'
                    }
                }
            }
        ]
    }

    fits = fits.HDUList(
        [fits.PrimaryHDU(),
         fits.ImageHDU(name='COEFFS', data=np.array([0.0], np.float32))])

    with JwstDataModel(fits, schema=schema) as dm:
        assert_array_equal(dm.coeffs, [0.0])


def test_casting():
    with RampModel((5, 5, 10, 10)) as dm:
        total = dm.data.sum()
        dm.data = dm.data + 2
        assert dm.data.sum() > total
