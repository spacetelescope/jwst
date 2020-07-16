import os
import shutil
import tempfile

import pytest

import numpy as np
from numpy.testing import assert_array_equal

from asdf import schema as mschema

from .. import DataModel, ImageModel, RampModel
from ..util import open

ROOT_DIR = None
FITS_FILE = None
TMP_FITS = None
TMP_FITS2 = None
TMP_YAML = None
TMP_JSON = None
TMP_DIR = None


def setup():
    global ROOT_DIR, FITS_FILE, TMP_DIR, TMP_FITS, TMP_YAML, TMP_JSON, TMP_FITS2
    ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')
    FITS_FILE = os.path.join(ROOT_DIR, 'test.fits')

    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = os.path.join(TMP_DIR, 'tmp.fits')
    TMP_YAML = os.path.join(TMP_DIR, 'tmp.yaml')
    TMP_JSON = os.path.join(TMP_DIR, 'tmp.json')
    TMP_FITS2 = os.path.join(TMP_DIR, 'tmp2.fits')


def teardown():
    shutil.rmtree(TMP_DIR)

def records_equal(a, b):
    a = a.item()
    b = b.item()
    a_size = len(a)
    b_size = len(b)
    equal = a_size == b_size
    for i in range(a_size):
        if not equal: break
        equal = a[i] == b[i]
    return equal

def test_from_new_hdulist():
    with pytest.raises(AttributeError):
        from astropy.io import fits
        hdulist = fits.HDUList()
        with open(hdulist) as dm:
            dm.data


def test_from_new_hdulist2():
    from astropy.io import fits
    hdulist = fits.HDUList()
    data = np.empty((50, 50), dtype=np.float32)
    primary = fits.PrimaryHDU()
    hdulist.append(primary)
    science = fits.ImageHDU(data=data, name='SCI')
    hdulist.append(science)
    with open(hdulist) as dm:
        dq = dm.dq
        assert dq is not None


def test_setting_arrays_on_fits():
    from astropy.io import fits
    hdulist = fits.HDUList()
    data = np.empty((50, 50), dtype=np.float32)
    primary = fits.PrimaryHDU()
    hdulist.append(primary)
    science = fits.ImageHDU(data=data, name='SCI')
    hdulist.append(science)
    with open(hdulist) as dm:
        dm.data = np.empty((50, 50), dtype=np.float32)
        dm.dq = np.empty((10, 50, 50), dtype=np.uint32)


def delete_array():
    with pytest.raises(AttributeError):
        from astropy.io import fits
        hdulist = fits.HDUList()
        data = np.empty((50, 50))
        science = fits.ImageHDU(data=data, name='SCI')
        hdulist.append(science)
        hdulist.append(science)
        with open(hdulist) as dm:
            del dm.data
            assert len(hdulist) == 1


def test_from_fits():
    with RampModel(FITS_FILE) as dm:
        assert dm.meta.instrument.name == 'MIRI'
        assert dm.shape == (5, 35, 40, 32)


def test_from_scratch():
    with ImageModel((50, 50)) as dm:
        data = np.asarray(np.random.rand(50, 50), np.float32)
        dm.data[...] = data

        dm.meta.instrument.name = 'NIRCAM'

        dm.to_fits(TMP_FITS, overwrite=True)

        with ImageModel.from_fits(TMP_FITS) as dm2:
            assert dm2.shape == (50, 50)
            assert dm2.meta.instrument.name == 'NIRCAM'
            assert dm2.dq.dtype.name == 'uint32'
            assert np.all(dm2.data == data)


def test_delete():
    with DataModel(FITS_FILE) as dm:
        dm.meta.instrument.name = 'NIRCAM'
        assert dm.meta.instrument.name == 'NIRCAM'
        del dm.meta.instrument.name
        assert dm.meta.instrument.name is None


# def test_section():
#     with RampModel((5, 35, 40, 32)) as dm:
#         section = dm.get_section('data')[3:4, 1:3]
#         assert section.shape == (1, 2, 40, 32)


# def test_date_obs():
#     with DataModel(FITS_FILE) as dm:
#         assert dm.meta.observation.date.microsecond == 314592


def test_fits_without_sci():
    from astropy.io import fits
    schema = {
        "allOf": [
            mschema.load_schema(
                os.path.join(os.path.dirname(__file__),
                             "../schemas/core.schema.yaml"),
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

    with DataModel(fits, schema=schema) as dm:
        assert_array_equal(dm.coeffs, [0.0])


def _header_to_dict(x):
    return dict((a, b) for (a, b, c) in x)


def test_extra_fits():
    path = os.path.join(ROOT_DIR, "headers.fits")

    assert os.path.exists(path)

    with DataModel(path) as dm:
        assert 'BITPIX' not in _header_to_dict(dm.extra_fits.PRIMARY.header)
        assert _header_to_dict(dm.extra_fits.PRIMARY.header)['SCIYSTRT'] == 705
        dm2 = dm.copy()
        dm2.to_fits(TMP_FITS, overwrite=True)

    with DataModel(TMP_FITS) as dm:
        assert 'BITPIX' not in _header_to_dict(dm.extra_fits.PRIMARY.header)
        assert _header_to_dict(dm.extra_fits.PRIMARY.header)['SCIYSTRT'] == 705


def test_hdu_order():
    from astropy.io import fits

    with ImageModel(data=np.array([[0.0]]),
                    dq=np.array([[0.0]]),
                    err=np.array([[0.0]])) as dm:
        dm.save(TMP_FITS)

    with fits.open(TMP_FITS, memmap=False) as hdulist:
        assert hdulist[1].header['EXTNAME'] == 'SCI'
        assert hdulist[2].header['EXTNAME'] == 'DQ'
        assert hdulist[3].header['EXTNAME'] == 'ERR'


def test_casting():
    with RampModel(FITS_FILE) as dm:
        sum = np.sum(dm.data)
        dm.data[:] = dm.data + 2
        assert np.sum(dm.data) > sum


# def test_comments():
#     with RampModel(FITS_FILE) as dm:
#         assert 'COMMENT' in (x[0] for x in dm._extra_fits.PRIMARY)
#         dm._extra_fits.PRIMARY.COMMENT = ['foobar']
#         assert dm._extra_fits.PRIMARY.COMMENT == ['foobar']


def test_fits_comments():
    with ImageModel() as dm:
        dm.meta.subarray.xstart = 42
        dm.save(TMP_FITS, overwrite=True)

    from astropy.io import fits
    with fits.open(TMP_FITS, memmap=False) as hdulist:
        header = hdulist[0].header
        find = ['Subarray parameters']
        found = 0

        for card in header.cards:
            if card[1] in find:
                found += 1

        assert found == len(find)


def test_metadata_doesnt_override():
    with ImageModel() as dm:
        dm.save(TMP_FITS, overwrite=True)

    from astropy.io import fits
    with fits.open(TMP_FITS, mode='update', memmap=False) as hdulist:
        hdulist[0].header['FILTER'] = 'F150W2'

    with ImageModel(TMP_FITS) as dm:
        assert dm.meta.instrument.filter == 'F150W2'


def test_table_with_metadata():
    schema = {
        "allOf": [
            mschema.load_schema(
                os.path.join(os.path.dirname(__file__),
                             "../schemas/core.schema.yaml"),
                resolve_references=True),
            {"type": "object",
            "properties": {
                "flux_table": {
                    "title": "Photometric flux conversion table",
                    "fits_hdu": "FLUX",
                    "datatype":
                    [
                        {"name": "parameter", "datatype": ['ascii', 7]},
                        {"name": "factor", "datatype": "float64"},
                        {"name": "uncertainty", "datatype": "float64"}
                    ]
                },
                "meta": {
                    "type": "object",
                    "properties": {
                        "fluxinfo": {
                            "title": "Information about the flux conversion",
                            "type": "object",
                            "properties": {
                                "exposure": {
                                    "title": "Description of exposure analyzed",
                                    "type": "string",
                                    "fits_hdu": "FLUX",
                                    "fits_keyword": "FLUXEXP"
                                }
                            }
                        }
                    }
                }
            }
         }
        ]
    }

    class FluxModel(DataModel):
        def __init__(self, init=None, flux_table=None, **kwargs):
            super(FluxModel, self).__init__(init=init, schema=schema, **kwargs)

            if flux_table is not None:
                self.flux_table = flux_table

    flux_im = [
        ('F560W', 1.0e-5, 1.0e-7),
        ('F770W', 1.1e-5, 1.6e-7),
        ]
    with FluxModel(flux_table=flux_im) as datamodel:
        datamodel.meta.fluxinfo.exposure = 'Exposure info'
        datamodel.save(TMP_FITS, overwrite=True)
        del datamodel

    from astropy.io import fits
    with fits.open(TMP_FITS, memmap=False) as hdulist:
        assert len(hdulist) == 3
        assert isinstance(hdulist[1], fits.BinTableHDU)
        assert hdulist[1].name == 'FLUX'
        assert hdulist[2].name == 'ASDF'


def test_replace_table():
    from astropy.io import fits

    schema_narrow = {
        "allOf": [
            mschema.load_schema(
                os.path.join(os.path.dirname(__file__),
                             "../schemas/core.schema.yaml"),
                resolve_references=True),
            {
                "type": "object",
                "properties": {
                    "data": {
                        "title": "relative sensitivity table",
                        "fits_hdu": "RELSENS",
                        "datatype": [
                            {"name": "TYPE", "datatype": ["ascii", 16]},
                            {"name": "T_OFFSET", "datatype": "float32"},
                            {"name": "DECAY_PEAK", "datatype": "float32"},
                            {"name": "DECAY_FREQ", "datatype": "float32"},
                            {"name": "TAU", "datatype": "float32"}
                        ]
                    }
                }
            }
        ]
    }

    schema_wide = {
        "allOf": [
            mschema.load_schema(
                os.path.join(os.path.dirname(__file__),
                             "../schemas/core.schema.yaml"),
                resolve_references=True),
            {
                "type": "object",
                "properties": {
                    "data": {
                        "title": "relative sensitivity table",
                        "fits_hdu": "RELSENS",
                        "datatype": [
                            {"name": "TYPE", "datatype": ["ascii", 16]},
                            {"name": "T_OFFSET", "datatype": "float64"},
                            {"name": "DECAY_PEAK", "datatype": "float64"},
                            {"name": "DECAY_FREQ", "datatype": "float64"},
                            {"name": "TAU", "datatype": "float64"}
                        ]
                    }
                }
            }
        ]
    }

    x = np.array([("string", 1., 2., 3., 4.)],
                 dtype=[('TYPE', 'S16'),
                        ('T_OFFSET', np.float32),
                        ('DECAY_PEAK', np.float32),
                        ('DECAY_FREQ', np.float32),
                        ('TAU', np.float32)])

    m = DataModel(schema=schema_narrow)
    m.data = x
    m.to_fits(TMP_FITS, overwrite=True)

    with fits.open(TMP_FITS, memmap=False) as hdulist:
        assert records_equal(x, np.asarray(hdulist[1].data))
        assert hdulist[1].data.dtype[1].str == '>f4'
        assert hdulist[1].header['TFORM2'] == 'E'

    with DataModel(TMP_FITS, schema=schema_wide) as m:
        m.to_fits(TMP_FITS2, overwrite=True)

    with fits.open(TMP_FITS2, memmap=False) as hdulist:
        assert records_equal(x, np.asarray(hdulist[1].data))
        assert hdulist[1].data.dtype[1].str == '>f8'
        assert hdulist[1].header['TFORM2'] == 'D'


def test_table_with_unsigned_int():
    schema = {
        'title': 'Test data model',
        '$schema': 'http://stsci.edu/schemas/fits-schema/fits-schema',
        'type': 'object',
        'properties': {
            'meta': {
                'type': 'object',
                'properties': {}
            },
            'test_table': {
                'title': 'Test table',
                'fits_hdu': 'TESTTABL',
                'datatype': [
                    {'name': 'FLOAT64_COL', 'datatype': 'float64'},
                    {'name': 'UINT32_COL', 'datatype': 'uint32'}
                ]
            }
        }
    }

    with DataModel(schema=schema) as dm:

        float64_info = np.finfo(np.float64)
        float64_arr = np.random.uniform(size=(10,))
        float64_arr[0] = float64_info.min
        float64_arr[-1] = float64_info.max

        uint32_info = np.iinfo(np.uint32)
        uint32_arr = np.random.randint(uint32_info.min, uint32_info.max + 1, size=(10,), dtype=np.uint32)
        uint32_arr[0] = uint32_info.min
        uint32_arr[-1] = uint32_info.max

        test_table = np.array(list(zip(float64_arr, uint32_arr)), dtype=dm.test_table.dtype)

        def assert_table_correct(model):
            for idx, (col_name, col_data) in enumerate([('float64_col', float64_arr), ('uint32_col', uint32_arr)]):
                # The table dtype and field dtype are stored separately, and may not
                # necessarily agree.
                assert np.can_cast(model.test_table.dtype[idx], col_data.dtype, 'equiv')
                assert np.can_cast(model.test_table.field(col_name).dtype, col_data.dtype, 'equiv')
                assert np.array_equal(model.test_table.field(col_name), col_data)

        # The datamodel casts our array to FITS_rec on assignment, so here we're
        # checking that the data survived the casting.
        dm.test_table = test_table
        assert_table_correct(dm)

        # Confirm that saving the table (and converting the uint32 values to signed int w/TZEROn)
        # doesn't mangle the data.
        dm.save(TMP_FITS)
        assert_table_correct(dm)

    # Confirm that the data loads from the file intact (converting the signed ints back to
    # the appropriate uint32 values).
    with DataModel(TMP_FITS, schema=schema) as dm2:
        assert_table_correct(dm2)


def test_metadata_from_fits():
    from astropy.io import fits

    mask = np.array([[0, 1], [2, 3]])
    fits.ImageHDU(data=mask, name='DQ').writeto(TMP_FITS, overwrite=True)
    with DataModel(init=TMP_FITS) as dm:
        dm.save(TMP_FITS2)

    with fits.open(TMP_FITS2, memmap=False) as hdulist:
        assert hdulist[2].name == 'ASDF'


# def test_float_as_int():
#     from astropy.io import fits

#     hdulist = fits.HDUList()
#     primary = fits.PrimaryHDU()
#     hdulist.append(primary)
#     hdulist[0].header['SUBSTRT1'] = 42.7
#     hdulist.writeto(TMP_FITS, overwrite=True)

#     with DataModel(TMP_FITS) as dm:
#         assert dm.meta.subarray.xstart == 42.7
