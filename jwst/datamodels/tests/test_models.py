from __future__ import absolute_import, unicode_literals, division, print_function

from astropy.extern import six
import datetime
import os
import shutil
import tempfile

import pytest
try:
    import yaml
    has_yaml = True
except ImportError:
    has_yaml = False

from astropy.time import Time

import numpy as np
from numpy.testing.decorators import knownfailureif
from numpy.testing import assert_array_equal

from .. import DataModel, ImageModel, QuadModel, MultiSlitModel, open
from .. import schema


FITS_FILE = None
TMP_FITS = None
TMP_YAML = None
TMP_JSON = None
TMP_DIR = None
TMP_FITS2 = None


def setup():
    global FITS_FILE, TMP_DIR, TMP_FITS, TMP_YAML, TMP_JSON, TMP_FITS2
    ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')
    FITS_FILE = os.path.join(ROOT_DIR, 'test.fits')

    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = os.path.join(TMP_DIR, 'tmp.fits')
    TMP_YAML = os.path.join(TMP_DIR, 'tmp.yaml')
    TMP_JSON = os.path.join(TMP_DIR, 'tmp.json')
    TMP_FITS2 = os.path.join(TMP_DIR, 'tmp2.fits')


def teardown():
    shutil.rmtree(TMP_DIR)


def test_set_shape():
    with pytest.raises(AttributeError):
        with ImageModel((50, 50)) as dm:
            assert dm.shape == (50, 50)
            dm.shape = (42, 23)


def test_broadcast():
    with ImageModel((50, 50)) as dm:
        data = np.empty((50,))
        dm.dq = data


def test_broadcast2():
    with ImageModel() as dm:
        data = np.empty((52, 50))
        dm.data = data

        dq = np.empty((50,))
        dm.dq = dq


def test_from_hdulist():
    from astropy.io import fits
    with fits.open(FITS_FILE) as hdulist:
        with open(hdulist) as dm:
            pass
        assert hdulist.fileinfo(0)['file'].closed == False


def delete_array():
    with open() as dm:
        del dm.data


def test_subarray():
    with DataModel(FITS_FILE) as dm:
        x = dm.meta.subarray.xstart


def roundtrip(func):
    def _create_source():
        dm = DataModel(FITS_FILE)

        assert dm.meta.instrument.name == 'MIRI'

        dm.meta.instrument.name = 'NIRCAM'
        dm.meta.subarray.xstart = 42
        return dm

    def _check_output(dm):
        assert dm.meta.instrument.name == 'NIRCAM'
        assert dm.meta.subarray.xstart == 42

    def test():
        with _create_source() as dm:
            with func(dm) as dm2:
                _check_output(dm2)

    test.__name__ = func.__name__

    return test


@roundtrip
def test_from_fits_write(dm):
    dm.to_fits(TMP_FITS, clobber=True)
    return DataModel.from_fits(TMP_FITS)


# @knownfailureif(not has_yaml)
# @roundtrip
# def test_from_fits_to_yaml(dm):
#     dm.to_yaml(TMP_YAML)
#     return DataModel.from_yaml(TMP_YAML)


# @roundtrip
# def test_from_fits_to_json(dm):
#     dm.to_json(TMP_JSON)
#     return DataModel.from_json(TMP_JSON)


def test_delete():
    with DataModel() as dm:
        dm.meta.instrument.name = 'NIRCAM'
        assert dm.meta.instrument.name == 'NIRCAM'
        del dm.meta.instrument.name
        assert dm.meta.instrument.name is None


def test_open():
    with open() as dm:
        pass

    with open((50, 50)) as dm:
        pass

    with open(FITS_FILE) as dm:
        assert isinstance(dm, QuadModel)


def test_copy():
    with ImageModel((50, 50)) as dm:
        dm.meta.instrument.name = "NIRCAM"
        dm.meta.foo = "BAR"

        with dm.copy() as dm2:
            dm2.data[0, 0] = 42
            assert np.sum(dm.data.flatten()) == 0

            assert dm2.meta.instrument.name == "NIRCAM"
            assert dm2.meta.foo == "BAR"
            dm2.meta.foo = "BAZ"
            assert dm.meta.foo == "BAR"
            dm2.meta.observation.obs_id = "FOO"
            assert dm.meta.observation.obs_id is None


def test_section():
    with QuadModel((5, 35, 40, 32)) as dm:
        section = dm.get_section('data')[3:4, 1:3]
        assert section.shape == (1, 2, 40, 32)


def test_init_with_array():
    array = np.empty((50, 50))
    with open(array) as dm:
        assert dm.data.shape == (50, 50)
        assert isinstance(dm, ImageModel)


def test_init_with_array2():
    array = np.empty((50, 50))
    with ImageModel(array) as dm:
        assert dm.data.shape == (50, 50)


def test_init_with_array3():
    with pytest.raises(ValueError):
        array = np.empty((50,))
        with ImageModel(array) as dm:
            pass


def test_set_array():
    with pytest.raises(ValueError):
        with ImageModel() as dm:
            data = np.empty((50,))
            dm.data = data


def test_set_array2():
    with ImageModel() as dm:
        data = np.empty((50, 50))
        dm.data = data


def test_base_model_has_no_arrays():
    with pytest.raises(AttributeError):
        with DataModel() as dm:
            dm.data


def test_array_type():
    with ImageModel() as dm:
        assert dm.dq.dtype == np.uint32


def test_copy_model():
    with DataModel() as dm:
        with DataModel(dm) as dm2:
            assert hasattr(dm2, 'meta')


def test_dtype_match():
    with ImageModel() as dm:
        dm.data = np.array([[1, 2, 3]], np.float32)


def test_multislit():
    with MultiSlitModel() as dm:
        dm.slits.append(dm.slits.item())
        slit = dm.slits[-1]
        slit.data = np.random.rand(5, 5)
        slit.dm = np.random.rand(5, 5)
        slit.err = np.random.rand(5, 5)


def test_secondary_shapes():
    with ImageModel((256, 256)) as dm:
        assert dm.dq.shape == (256, 256)
        dm.dq = [1]


def test_multislit_from_image():
    with ImageModel((64, 64)) as im:
        with MultiSlitModel(im) as ms:
            assert len(ms.slits) == 1
            assert ms.slits[0].data.shape == (64, 64)


def test_multislit_from_fits_image():
    with ImageModel((64, 64)) as im:
        im.save(TMP_FITS, clobber=True)

    with MultiSlitModel(TMP_FITS) as ms:
        assert len(ms.slits) == 1
        assert ms.slits[0].data.shape == (64, 64)

        for i, slit in enumerate(ms.slits):
            assert slit.data is ms.slits[i].data

        ms2 = ms.copy()
        ms2.slits = ms.slits
        assert len(ms2.slits) == 1


def test_multislit_metadata():
    with MultiSlitModel() as ms:
        with ImageModel((64, 64)) as im:
            ms.slits.append(ms.slits.item())
            ms.slits[-1].data = im.data
        im = ms.slits[0]
        im.subarray.name = "FULL"
        assert ms.slits[0].subarray.name == "FULL"


def test_multislit_metadata():
    with MultiSlitModel() as ms:
        ms.slits.append(ms.slits.item())
        for key, val in ms.iteritems():
            if six.PY2:
                assert isinstance(val, (bytes, unicode, int, long, float, bool,
                                        Time))
            else:
                assert isinstance(val, (bytes, str, int, float, bool, Time))


def test_multislit_copy():
    with MultiSlitModel() as input:
        for i in range(4):
            input.slits.append(input.slits.item(
                data=np.empty((50, 50), dtype=np.float32)))

        i = 0
        for slit in input.slits:
            i += 1
        assert i == 4

        input.save(TMP_FITS)

        output = input.copy()
        assert len(output.slits) == 4

        i = 0
        for slit in output.slits:
            i += 1
        assert i == 4

    from astropy.io import fits
    with fits.open(TMP_FITS) as hdulist:
        assert len(hdulist) == 6

    with MultiSlitModel(TMP_FITS) as input:
        i = 0
        for slit in input.slits:
            i += 1
        assert i == 4

        output = input.copy()
        assert len(output.slits) == 4


def test_model_with_nonstandard_primary_array():
    ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')

    class NonstandardPrimaryArrayModel(DataModel):
        schema_url = os.path.join(
            ROOT_DIR, "nonstandard_primary_array.schema.yaml")

        def __init__(self, init=None, wavelength=None, alpha=None, **kwargs):
            super(NonstandardPrimaryArrayModel, self).__init__(init=init, **kwargs)

            if wavelength is not None:
                self.wavelength = wavelength

            if alpha is not None:
                self.alpha = alpha

        # The wavelength array is the primary array.
        # Try commenting this function out and the problem goes away.
        def get_primary_array_name(self):
            return 'wavelength'

    m = NonstandardPrimaryArrayModel()
    list(m.keys())


def test_relsens():
    with ImageModel() as im:
        assert len(im.relsens.dtype) == 2


def test_image_with_extra_keyword_to_multislit():
    with ImageModel(data=np.empty((32, 32))) as im:
        im.save(TMP_FITS, clobber=True)

    from astropy.io import fits
    with fits.open(TMP_FITS, mode="update") as hdulist:
        hdulist[1].header['BUNIT'] = 'x'

    with ImageModel(TMP_FITS) as im:
        with MultiSlitModel() as ms:
            ms.update(im)
            for i in range(3):
                ms.slits.append(ImageModel(data=np.empty((4, 4))))
            assert len(ms.slits) == 3

            ms.save(TMP_FITS2, clobber=True)

    with MultiSlitModel(TMP_FITS2) as ms:
        assert len(ms.slits) == 3
        for slit in ms.slits:
            assert slit.data.shape == (4, 4)
