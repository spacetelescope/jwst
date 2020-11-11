from datetime import datetime
import os

from asdf import schema as mschema
from astropy import time
from astropy.io import fits
from numpy.testing import assert_array_equal, assert_array_almost_equal
from stdatamodels import validate
import jsonschema
import numpy as np
import pytest

from jwst.datamodels import util
from jwst.datamodels import _defined_models as defined_models
from jwst.datamodels import (JwstDataModel, ImageModel, MaskModel, MultiSlitModel,
    AsnModel, SourceModelContainer, MultiExposureModel,
    DrizProductModel, MultiProductModel, MIRIRampModel)


def test_strict_validation_enum():
    with JwstDataModel(strict_validation=True) as dm:
        assert dm.meta.instrument.name is None
        with pytest.raises(jsonschema.ValidationError):
            # FOO is not in the allowed enumerated values
            dm.meta.instrument.name = 'FOO'


def test_strict_validation_type():
    with JwstDataModel(strict_validation=True) as dm:
        with pytest.raises(jsonschema.ValidationError):
            # Schema requires a float
            dm.meta.target.ra = "FOO"


def test_strict_validation_date():
    with JwstDataModel(strict_validation=True) as dm:
        time_obj = time.Time(dm.meta.date)
        assert isinstance(time_obj, time.Time)
        date_obj = datetime.strptime(dm.meta.date, '%Y-%m-%dT%H:%M:%S.%f')
        assert isinstance(date_obj, datetime)


TRANSFORMATION_SCHEMA = {
    "allOf": [
        mschema.load_schema("http://stsci.edu/schemas/jwst_datamodel/core.schema",
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
    with JwstDataModel(schema=TRANSFORMATION_SCHEMA,
                       strict_validation=True) as dm:
        dm.meta.transformations = []
        with pytest.raises(jsonschema.ValidationError):
            dm.meta.transformations.item(transformation="SIN", coeff=2.0)


def test_list2():
    with JwstDataModel(schema=TRANSFORMATION_SCHEMA,
                       strict_validation=True) as dm:
        dm.meta.transformations = []
        with pytest.raises(jsonschema.ValidationError):
            dm.meta.transformations.append({'transformation': 'FOO', 'coeff': 2.0})


def test_invalid_fits(tmp_path):
    path = str(tmp_path / "invalid.fits")

    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    header = hdulist[0].header
    header["DATAMODL"] = "JwstDataModel"
    # Add invalid keyword
    header["INSTRUME"] = "FOO"
    hdulist.writeto(path)

    os.environ['PASS_INVALID_VALUES'] = '0'
    os.environ['STRICT_VALIDATION'] = '0'
    with pytest.warns(validate.ValidationWarning):
        with util.open(path):
            pass

    os.environ['STRICT_VALIDATION'] = '1'
    with pytest.raises(jsonschema.ValidationError):
        with util.open(path):
            pass

    # Check that specifying an argument does not get
    # overridden by the environmental.
    os.environ['STRICT_VALIDATION'] = '0'
    with pytest.raises(jsonschema.ValidationError):
        with util.open(path, strict_validation=True):
            pass

    os.environ['PASS_INVALID_VALUES'] = '0'
    os.environ['STRICT_VALIDATION'] = '0'
    with pytest.warns(validate.ValidationWarning):
        with util.open(path, pass_invalid_values=True) as model:
            assert model.meta.instrument.name == 'FOO'

    os.environ['PASS_INVALID_VALUES'] = '0'
    os.environ['STRICT_VALIDATION'] = '0'
    with pytest.warns(validate.ValidationWarning):
        with util.open(path) as model:
            assert model.meta.instrument.name is None

    os.environ['PASS_INVALID_VALUES'] = '1'
    with pytest.warns(validate.ValidationWarning):
        with util.open(path) as model:
            assert model.meta.instrument.name == 'FOO'

    del os.environ['PASS_INVALID_VALUES']
    del os.environ['STRICT_VALIDATION']

    with pytest.warns(validate.ValidationWarning):
        with util.open(path, pass_invalid_values=False, strict_validation=False):
            pass

    with pytest.raises(jsonschema.ValidationError):
        with util.open(path, pass_invalid_values=False, strict_validation=True):
            pass

    with pytest.warns(validate.ValidationWarning):
        with util.open(path, pass_invalid_values=False, strict_validation=False) as model:
            assert model.meta.instrument.name is None

    with pytest.warns(validate.ValidationWarning):
        with util.open(path, pass_invalid_values=True, strict_validation=False) as model:
            pass
    assert model.meta.instrument.name == 'FOO'


def test_mask_model():
    with MaskModel((10, 10)) as dm:
        assert dm.dq.dtype == np.uint32


def test_data_array(tmp_path):
    path = str(tmp_path / "data_array.fits")
    data_array_schema = {
        "allOf": [
            mschema.load_schema("http://stsci.edu/schemas/jwst_datamodel/core.schema",
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

    with JwstDataModel(schema=data_array_schema) as x:
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
        x.save(path)

    with JwstDataModel(path, schema=data_array_schema) as x:
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
        x.save(path)

    from astropy.io import fits
    with fits.open(path) as hdulist:
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


def test_multislit_move_from_fits(tmp_path):
    path = tmp_path / "multislit_move.fits"

    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    for i in range(5):
        hdu = fits.ImageHDU(data=np.zeros((64, 64)), name='SCI')
        hdu.ver = i + 1
        hdulist.append(hdu)

    hdulist.writeto(path)

    n = MultiSlitModel()
    with MultiSlitModel(path) as m:
        n.slits.append(m.slits[2])

        assert len(n.slits) == 1


def test_multislit_append_string():
    with pytest.raises(jsonschema.ValidationError):
        m = MultiSlitModel(strict_validation=True)
        m.slits.append('junk')


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
    """Test that entry points for JwstDataModelExtension works as expected"""
    mschema.load_schema('http://stsci.edu/schemas/jwst_datamodel/image.schema',
        resolve_references=True)
