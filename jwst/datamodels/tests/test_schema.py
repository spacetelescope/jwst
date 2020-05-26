# Licensed under a 3-clause BSD style license - see LICENSE.rst

from datetime import datetime
import os
import shutil
import tempfile
import warnings

import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

import jsonschema
from astropy.io import fits
from astropy.modeling import models
from astropy import time

from .. import util, validate
from .. import _defined_models as defined_models
from .. import (DataModel, ImageModel, RampModel, MaskModel, MultiSlitModel,
    AsnModel, CollimatorModel, SourceModelContainer, MultiExposureModel,
    DrizProductModel, MultiProductModel, MIRIRampModel)
from ..schema import merge_property_trees, build_docstring

from ..extension import URL_PREFIX

import asdf
from asdf import schema as mschema

FITS_FILE = None
MASK_FILE = None
TMP_FITS = None
TMP_FITS2 = None
TMP_YAML = None
TMP_ASDF = None
TMP_DIR = None

def setup():
    global FITS_FILE, MASK_FILE, TMP_DIR, TMP_FITS, TMP_YAML, TMP_ASDF, TMP_FITS2
    ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')
    FITS_FILE = os.path.join(ROOT_DIR, 'test.fits')
    MASK_FILE = os.path.join(ROOT_DIR, 'mask.fits')

    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = os.path.join(TMP_DIR, 'tmp.fits')
    TMP_FITS2 = os.path.join(TMP_DIR, 'tmp2.fits')
    TMP_YAML = os.path.join(TMP_DIR, 'tmp.yaml')
    TMP_ASDF = os.path.join(TMP_DIR, 'tmp.asdf')


def teardown():
    shutil.rmtree(TMP_DIR)


def test_choice():
    with pytest.raises(jsonschema.ValidationError):
        with DataModel(FITS_FILE, strict_validation=True) as dm:
            assert dm.meta.instrument.name == 'MIRI'
            dm.meta.instrument.name = 'FOO'


def test_set_na_ra():
    with pytest.raises(jsonschema.ValidationError):
        with DataModel(FITS_FILE, strict_validation=True) as dm:
            # Setting an invalid value should raise a ValueError
            dm.meta.target.ra = "FOO"


def test_date2():
    with ImageModel((50, 50), strict_validation=True) as dm:
        time_obj = time.Time(dm.meta.date)
        assert isinstance(time_obj, time.Time)
        date_obj = datetime.strptime(dm.meta.date, '%Y-%m-%dT%H:%M:%S.%f')
        assert isinstance(date_obj, datetime)


TRANSFORMATION_SCHEMA = {
    "allOf": [
        mschema.load_schema(os.path.join(URL_PREFIX, "image.schema"),
            resolver=asdf.AsdfFile().resolver,
            resolve_references=True),
        {
            "type": "object",
            "properties": {
                "meta": {
                    "type": "object",
                    "properties": {
                        "transformations": {
                            "title": "A list of transformations",
                            "type": "array",
                            "items": {
                                "title": "A transformation",
                                "type": "object",
                                "properties": {
                                    "type": {
                                        "title": "Transformation type",
                                        "type": "string"
                                    },
                                    "coeff": {
                                        "title": "coefficients",
                                        "type": "number"
                                    }
                                },
                                "additionalProperties": False
                            }
                        }
                    }
                }
            }
        }
    ]
}


def test_list():
    with pytest.raises(jsonschema.ValidationError):
        with ImageModel((50, 50), schema=TRANSFORMATION_SCHEMA,
                        strict_validation=True) as dm:
            dm.meta.transformations = []
            dm.meta.transformations.item(transformation="SIN", coeff=2.0)


def test_list2():
    with pytest.raises(jsonschema.ValidationError):
        with ImageModel(
            (50, 50),
            schema=TRANSFORMATION_SCHEMA,
            strict_validation=True) as dm:
            dm.meta.transformations = []
            dm.meta.transformations.append({'transformation': 'FOO', 'coeff': 2.0})


def test_invalid_fits():
    hdulist = fits.open(FITS_FILE)
    header = hdulist[0].header
    header['INSTRUME'] = 'FOO'

    if os.path.exists(TMP_FITS):
        os.remove(TMP_FITS)

    hdulist.writeto(TMP_FITS)
    hdulist.close()

    with pytest.raises(validate.ValidationWarning):
        with warnings.catch_warnings():
            os.environ['PASS_INVALID_VALUES'] = '0'
            os.environ['STRICT_VALIDATION'] = '0'
            warnings.simplefilter('error')
            model = util.open(TMP_FITS)
            model.close()

    with pytest.raises(jsonschema.ValidationError):
        os.environ['STRICT_VALIDATION'] = '1'
        model = util.open(TMP_FITS)
        model.close()

    # Check that specifying an argument does not get
    # overridden by the environmental.
    with pytest.raises(jsonschema.ValidationError):
        os.environ['STRICT_VALIDATION'] = '0'
        model = util.open(TMP_FITS, strict_validation=True)
        model.close()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        os.environ['PASS_INVALID_VALUES'] = '0'
        os.environ['STRICT_VALIDATION'] = '0'
        model = util.open(TMP_FITS, pass_invalid_values=True)
        assert model.meta.instrument.name == 'FOO'
        model.close()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        os.environ['PASS_INVALID_VALUES'] = '0'
        os.environ['STRICT_VALIDATION'] = '0'
        model = util.open(TMP_FITS)
        assert model.meta.instrument.name is None
        model.close()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        os.environ['PASS_INVALID_VALUES'] = '1'
        model = util.open(TMP_FITS)
        assert model.meta.instrument.name == 'FOO'
        model.close()

    del os.environ['PASS_INVALID_VALUES']
    del os.environ['STRICT_VALIDATION']

    with pytest.raises(validate.ValidationWarning):
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            model = util.open(TMP_FITS,
                              pass_invalid_values=False,
                              strict_validation=False)
            model.close()

    with pytest.raises(jsonschema.ValidationError):
        model = util.open(TMP_FITS,
                          pass_invalid_values=False,
                          strict_validation=True)
        model.close()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = util.open(TMP_FITS,
                          pass_invalid_values=False,
                          strict_validation=False)
        assert model.meta.instrument.name is None
        model.close()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = util.open(TMP_FITS,
                          pass_invalid_values=True,
                          strict_validation=False)
        assert model.meta.instrument.name == 'FOO'
        model.close()

    if os.path.exists(TMP_FITS):
        os.remove(TMP_FITS)

def test_ad_hoc_json():
    with DataModel() as dm:
        dm.meta.foo = {'a': 42, 'b': ['a', 'b', 'c']}

        dm.save(TMP_ASDF)

    with DataModel(TMP_ASDF) as dm2:
        assert dm2.meta.foo == {'a': 42, 'b': ['a', 'b', 'c']}


def test_ad_hoc_fits():
    with DataModel() as dm:
        dm.meta.foo = {'a': 42, 'b': ['a', 'b', 'c']}

        dm.to_fits(TMP_FITS, overwrite=True)

    with DataModel.from_fits(TMP_FITS) as dm2:
        assert dm2.meta.foo == {'a': 42, 'b': ['a', 'b', 'c']}


def test_find_fits_keyword():
    with DataModel() as x:
        assert x.find_fits_keyword('DATE-OBS') == \
          ['meta.observation.date']


def test_search_schema():
    with DataModel() as x:
        results = x.search_schema('target')

    results = [x[0] for x in results]
    assert 'meta.target' in results
    assert 'meta.target.ra' in results


def test_dictionary_like():
    with DataModel(strict_validation=True) as x:
        x.meta.origin = 'FOO'
        assert x['meta.origin'] == 'FOO'

        with pytest.raises(jsonschema.ValidationError):
            x['meta.subarray.xsize'] = 'FOO'

        with pytest.raises(KeyError):
            x['meta.FOO.BAR.BAZ']


def test_to_flat_dict():
    with DataModel() as x:
        x.meta.origin = 'FOO'
        assert x['meta.origin'] == 'FOO'

        d = x.to_flat_dict()

        assert d['meta.origin'] == 'FOO'


def test_table_array_shape_ndim():
    table_schema = {
        "allOf": [
            mschema.load_schema(os.path.join(URL_PREFIX, "image.schema"),
                resolver=asdf.AsdfFile().resolver,
                resolve_references=True),
            {
                "type": "object",
                "properties": {
                    "table": {
                        'title': 'A structured table',
                        'fits_hdu': 'table',
                        'datatype': [
                            'bool8',
                            {'datatype': 'int16',
                             'name': 'my_int'},
                            {'datatype': 'float32',
                             'name': 'my_float1',
                             'shape': [3, 2]},
                            {'datatype': 'float32',
                             'name': 'my_float2',
                             'ndim': 2},
                            {'datatype': 'float32',
                             'name': 'my_float3'},
                            {'datatype': 'float32',
                             'name': 'my_float4'},
                            {'datatype': 'float32',
                             'name': 'my_float5'},
                            {'datatype': ['ascii', 64],
                             'name': 'my_string'}
                        ]
                    }
                }
            }
        ]
    }

    with DataModel(schema=table_schema) as x:
        x.table = [
            (
                True,
                42,
                [[37.5, 38.0], [39.0, 40.0], [41.0, 42.0]],
                [[37.5, 38.0], [39.0, 40.0], [41.0, 42.0]],
                [[37.5, 38.0], [39.0, 40.0], [41.0, 42.0]],
                [37.5, 38.0],
                37.5,
                'STRING'
            )
        ]
        assert x.table.dtype == [
            ('f0', '?'),
            ('my_int', '=i2'),
            ('my_float1', '=f4', (3, 2)),
            ('my_float2', '=f4', (3, 2)),
            ('my_float3', '=f4', (3, 2)),
            ('my_float4', '=f4', (2,)),
            ('my_float5', '=f4'),
            ('my_string', 'S64')
        ]

        x.to_fits(TMP_FITS, overwrite=True)

    with DataModel(TMP_FITS, schema=table_schema) as x:
        assert x.table.dtype == [
            ('f0', '?'),
            ('my_int', '=i2'),
            ('my_float1', '=f4', (3, 2)),
            ('my_float2', '=f4', (3, 2)),
            ('my_float3', '=f4', (3, 2)),
            ('my_float4', '=f4', (2,)),
            ('my_float5', '=f4'),
            ('my_string', 'S64')
        ]

    table_schema['allOf'][1]['properties']['table']['datatype'][3]['ndim'] = 3
    with DataModel(schema=table_schema) as x:
        with pytest.raises(ValueError) as e:
            x.table = [
                (
                    True,
                    42,
                    [[37.5, 38.0], [39.0, 40.0], [41.0, 42.0]],
                    [[37.5, 38.0], [39.0, 40.0], [41.0, 42.0]],
                    [[37.5, 38.0], [39.0, 40.0], [41.0, 42.0]],
                    [37.5, 38.0],
                    37.5,
                    'STRING'
                )
            ]

        assert str(e.value).startswith("Array has wrong number of dimensions.")


def test_table_array_convert():
    """
    Test that structured arrays are converted when necessary, and
    reused as views when not.
    """
    from jwst.datamodels import util

    table_schema = {
        "allOf": [
            mschema.load_schema(os.path.join(URL_PREFIX, "image.schema"),
                resolver=asdf.AsdfFile().resolver,
                resolve_references=True),
            {
                "type": "object",
                "properties": {
                    "table": {
                        'title': 'A structured table',
                        'fits_hdu': 'table',
                        'datatype': [
                            'bool8',
                            {'datatype': 'int16',
                             'name': 'my_int'},
                            {'datatype': ['ascii', 64],
                             'name': 'my_string'}
                        ]
                    }
                }
            }
        ]
    }

    table = np.array(
        [(42, 32000, 'foo')],
        dtype=[('f0', '?'), ('my_int', '=i2'), ('my_string', 'S64')]
    )

    x = util.gentle_asarray(
        table,
        dtype=[('f0', '?'), ('my_int', '=i2'), ('my_string', 'S64')]
    )

    assert x is table

    with DataModel(schema=table_schema) as x:
        x.table = table
        assert x.table is not table

    table = np.array(
        [(42, 32000, 'foo')],
        dtype=[('f0', '?'), ('my_int', '=i2'), ('my_string', 'S3')]
    )

    with DataModel(schema=table_schema) as x:
        x.table = table
        assert x.table is not table
        assert x.table['my_string'][0] != table['my_string'][0]


def test_mask_model():
    with MaskModel(MASK_FILE) as dm:
        assert dm.dq.dtype == np.uint32


def test_data_array():
    data_array_schema = {
        "allOf": [
            mschema.load_schema(os.path.join(URL_PREFIX, "core.schema"),
                resolver=asdf.AsdfFile().resolver,
                resolve_references=True),
            {
                "type": "object",
                "properties": {
                    "arr": {
                        'title': 'An array of data',
                        'type': 'array',
                        "fits_hdu": ["FOO", "DQ"],

                        "items": {
                            "title": "entry",
                            "type": "object",
                            "properties": {
                                "data": {
                                    "fits_hdu": "FOO",
                                    "default": 0.0,
                                    "max_ndim": 2,
                                    "datatype": "float64"
                                },
                                "dq": {
                                    "fits_hdu": "DQ",
                                    "default": 1,
                                    "datatype": "uint8"
                                },
                            }
                        }
                    }
                }
            }
        ]
    }

    array1 = np.random.rand(5, 5)
    array2 = np.random.rand(5, 5)
    array3 = np.random.rand(5, 5)

    with DataModel(schema=data_array_schema) as x:
        x.arr.append(x.arr.item())
        x.arr[0].data = array1
        assert len(x.arr) == 1
        x.arr.append(x.arr.item(data=array2))
        assert len(x.arr) == 2
        x.arr.append({})
        assert len(x.arr) == 3
        x.arr[2].data = array3
        del x.arr[1]
        assert len(x.arr) == 2
        x.to_fits(TMP_FITS, overwrite=True)

    with DataModel(TMP_FITS, schema=data_array_schema) as x:
        assert len(x.arr) == 2
        assert_array_almost_equal(x.arr[0].data, array1)
        assert_array_almost_equal(x.arr[1].data, array3)

        del x.arr[0]
        assert len(x.arr) == 1

        x.arr = []
        assert len(x.arr) == 0
        x.arr.append({'data': np.empty((5, 5))})
        assert len(x.arr) == 1
        x.arr.extend([
            x.arr.item(data=np.empty((5, 5))),
            x.arr.item(data=np.empty((5, 5)),
                       dq=np.empty((5, 5), dtype=np.uint8))])
        assert len(x.arr) == 3
        del x.arr[1]
        assert len(x.arr) == 2
        x.to_fits(TMP_FITS2, overwrite=True)

    from astropy.io import fits
    with fits.open(TMP_FITS2) as hdulist:
        x = set()
        for hdu in hdulist:
            x.add((hdu.header.get('EXTNAME'),
                   hdu.header.get('EXTVER')))

        assert x == set(
            [('FOO', 2), ('FOO', 1), ('ASDF', None), ('DQ', 2),
             (None, None)])


def test_multislit_model():
    array1 = np.asarray(np.random.rand(2, 2), dtype='float32')
    array2 = np.asarray(np.random.rand(2, 2), dtype='float32')

    with MultiSlitModel() as ms:
        assert len(ms.slits) == 0
        ms.slits.append(ms.slits.item())
        ms.slits[-1].data = array1
        assert len(ms.slits) == 1
        ms.slits.append(ms.slits.item())
        ms.slits[-1].data = array2
        assert len(ms.slits) == 2
        for i, slit in enumerate(ms.slits):
            assert slit == ms.slits[i]
        ms2 = ms.copy()
        assert len(ms2.slits) == 2
        assert_array_equal(ms.slits[-1].data, array2)
        del ms.slits[0]
        assert len(ms.slits) == 1
        assert_array_equal(ms.slits[0].data, array2)


def test_implicit_creation_lower_dimensionality():
    with RampModel(np.zeros((10, 20, 30, 40))) as rm:
        assert rm.pixeldq.shape == (30, 40)


def test_add_schema_entry():
    with DataModel(strict_validation=True) as dm:
        dm.add_schema_entry('meta.foo.bar', {'enum': ['foo', 'bar', 'baz']})
        dm.meta.foo.bar
        dm.meta.foo.bar = 'bar'
        try:
            dm.meta.foo.bar = 'what?'
        except jsonschema.ValidationError:
            pass
        else:
            assert False


def test_table_size_zero():
    with AsnModel() as dm:
        assert len(dm.asn_table) == 0


def test_copy_multslit():
    model1 = MultiSlitModel()
    model2 = MultiSlitModel()

    model1.slits.append(ImageModel(np.ones((1024, 1024))))
    model2.slits.append(ImageModel(np.ones((1024, 1024)) * 2))

    # Create the ouput model as a copy of the first input
    output = model1.copy()

    assert len(model1.slits) == 1
    assert len(model2.slits) == 1
    assert len(output.slits) == 1

    assert model1.slits[0].data[330, 330] == 1
    assert output.slits[0].data[330, 330] == 1
    assert id(model1.slits[0].data) != id(output.slits[0].data)

    output.slits[0].data = model1.slits[0].data - model2.slits[0].data

    assert model1.slits[0].data[330, 330] == 1
    assert output.slits[0].data[330, 330] == -1


def test_multislit_move_from_fits():
    from astropy.io import fits

    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    for i in range(5):
        hdu = fits.ImageHDU(data=np.zeros((64, 64)), name='SCI')
        hdu.ver = i + 1
        hdulist.append(hdu)

    hdulist.writeto(TMP_FITS, overwrite=True)

    n = MultiSlitModel()
    with MultiSlitModel(TMP_FITS) as m:
        n.slits.append(m.slits[2])

        assert len(n.slits) == 1


def test_validate_transform():
    """
    Tests that custom types, like transform, can be validated.
    """
    m = CollimatorModel(model=models.Shift(1) & models.Shift(2),
                        strict_validation=True)
    m.meta.description = "Test validate a WCS reference file."
    m.meta.author = "ND"
    m.meta.pedigree = "GROUND"
    m.meta.useafter = "2018/06/18"
    m.meta.reftype = "collimator"
    m.validate()


def test_validate_transform_from_file():
    """
    Tests that custom types, like transform, can be validated.
    """
    fname = os.path.join(os.path.dirname(__file__), 'data', 'collimator_fake.asdf')
    with CollimatorModel(fname, strict_validation=True) as m:
        m.validate()


def test_multislit_append_string():
    with pytest.raises(jsonschema.ValidationError):
        m = MultiSlitModel(strict_validation=True)
        m.slits.append('junk')


@pytest.mark.parametrize('combiner', ['anyOf', 'oneOf'])
def test_merge_property_trees(combiner):

    s = {
         'type': 'object',
         'properties': {
             'foobar': {
                 combiner: [
                     {
                         'type': 'array',
                         'items': [ {'type': 'string'}, {'type': 'number'} ],
                         'minItems': 2,
                         'maxItems': 2,
                     },
                     {
                         'type': 'array',
                         'items': [
                             {'type': 'number'},
                             {'type': 'string'},
                             {'type': 'number'}
                         ],
                         'minItems': 3,
                         'maxItems': 3,
                     }
                 ]
             }
         }
    }

    # Make sure that merge_property_trees does not destructively modify schemas
    f = merge_property_trees(s)
    assert f == s


def test_schema_docstring():
    template = "{fits_hdu} {title}"
    docstring = build_docstring(ImageModel, template).split("\n")
    for i, hdu in enumerate(('SCI', 'DQ', 'ERR', 'ZEROFRAME')):
        assert docstring[i].startswith(hdu)


@pytest.mark.parametrize("model", [v for v in defined_models.values()])
def test_all_datamodels_init(model):
    """
    Test that all current datamodels can be initialized.
    """
    if model is SourceModelContainer:
        # SourceModelContainer cannot have init=None
        model(MultiExposureModel())
    elif model in (DrizProductModel, MultiProductModel, MIRIRampModel):
        with pytest.warns(DeprecationWarning):
            model()
    else:
        model()


def test_datamodel_schema_entry_points():
    """Test that entry points for DataModelExtension works as expected"""
    resolver = asdf.AsdfFile().resolver
    mschema.load_schema('http://stsci.edu/schemas/jwst_datamodel/image.schema',
        resolver=resolver, resolve_references=True)
