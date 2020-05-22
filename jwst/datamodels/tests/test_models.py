from copy import copy
import os
from os import path as op
import shutil
import tempfile
import warnings
import jsonschema

import pytest
from astropy.table import Table
from astropy.time import Time
import numpy as np
from numpy.testing import assert_allclose
from astropy.io import fits

from jwst.datamodels import (DataModel, ImageModel, MaskModel, QuadModel,
                MultiSlitModel, ModelContainer, SlitModel,
                SlitDataModel, IFUImageModel, ABVegaOffsetModel)
from jwst import datamodels
from jwst.lib.file_utils import pushdir


ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')
FITS_FILE = os.path.join(ROOT_DIR, 'test.fits')
ASDF_FILE = os.path.join(ROOT_DIR, 'collimator_fake.asdf')
ASN_FILE = os.path.join(ROOT_DIR, 'association.json')


def reset_group_id(container):
    """Remove group_id from all models in container"""
    for m in container:
        try:
            del m.meta.group_id
        except AttributeError:
            pass


def setup():
    global FITS_FILE, TMP_DIR, TMP_FITS, TMP_YAML, TMP_JSON, TMP_FITS2

    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = os.path.join(TMP_DIR, 'tmp.fits')
    TMP_YAML = os.path.join(TMP_DIR, 'tmp.yaml')
    TMP_JSON = os.path.join(TMP_DIR, 'tmp.json')
    TMP_FITS2 = os.path.join(TMP_DIR, 'tmp2.fits')


def teardown():
    shutil.rmtree(TMP_DIR)


def test_set_shape():
    with ImageModel((50, 50)) as dm:
        assert dm.shape == (50, 50)

        with pytest.raises(AttributeError):
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
    warnings.simplefilter("ignore")
    with fits.open(FITS_FILE, memmap=False) as hdulist:
        with datamodels.open(hdulist) as dm:
            dm.data
        assert hdulist.fileinfo(0)['file'].closed == False


def delete_array():
    with datamodels.open() as dm:
        del dm.data


def test_subarray():
    with DataModel(FITS_FILE) as dm:
        dm.meta.subarray.xstart


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
    dm.to_fits(TMP_FITS, overwrite=True)
    return DataModel.from_fits(TMP_FITS)


def test_delete():
    with DataModel() as dm:
        dm.meta.instrument.name = 'NIRCAM'
        assert dm.meta.instrument.name == 'NIRCAM'
        del dm.meta.instrument.name
        assert dm.meta.instrument.name is None


def test_open():
    with datamodels.open() as dm:
        pass

    with datamodels.open((50, 50)) as dm:
        pass

    warnings.simplefilter("ignore")
    with datamodels.open(FITS_FILE) as dm:
        assert isinstance(dm, QuadModel)

def test_open_warning():
    with warnings.catch_warnings(record=True) as warners:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        with datamodels.open(FITS_FILE) as model:
            pass

        class_name = model.__class__.__name__
        j = None
        for i, w in enumerate(warners):
            if class_name in str(w.message):
                j = i
        assert j is not None


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

def test_stringify():
    im = DataModel()
    assert str(im) == '<DataModel>'
    im.close()

    im = ImageModel((10,100))
    assert str(im) == '<ImageModel(10, 100)>'
    im.close()

    image = ROOT_DIR + "/nircam_mask.fits"
    im = MaskModel(image)
    assert str(im) ==  '<MaskModel(2048, 2048) from nircam_mask.fits>'
    im.close()

def test_section():
    with QuadModel((5, 35, 40, 32)) as dm:
        section = dm.get_section('data')[3:4, 1:3]
        assert section.shape == (1, 2, 40, 32)


def test_init_with_array():
    array = np.empty((50, 50))
    with datamodels.open(array) as dm:
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
            dm.data


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


def test_default_value_anyof_schema():
    """Make sure default values are set properly when anyOf in schema"""
    with ImageModel((100, 100)) as im:
        val = im.meta.instrument.channel
        assert val is None
        val = im.meta.observation.date
        assert val is None
        val = im.meta.guidestar
        assert val.instance == {}


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
        im.save(TMP_FITS, overwrite=True)

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
        im.name = "FOO"
        assert ms.slits[0].name == "FOO"


def test_multislit_metadata2():
    with MultiSlitModel() as ms:
        ms.slits.append(ms.slits.item())
        for key, val in ms.iteritems():
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
    with fits.open(TMP_FITS, memmap=False) as hdulist:
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

def test_initialize_arrays_with_arglist():
    shape = (10,10)
    thirteen = np.full(shape, 13.0)
    bitz = np.full(shape, 7)

    im = ImageModel(shape, zeroframe=thirteen, dq=bitz)
    assert np.array_equal(im.zeroframe, thirteen)
    assert np.array_equal(im.dq, bitz)

def test_open_asdf_model():
    # Open an empty asdf file, pass extra arguments

    model = DataModel(init=None,
                     ignore_version_mismatch=False,
                     ignore_unrecognized_tag=True)

    assert model._asdf._ignore_version_mismatch == False
    assert model._asdf._ignore_unrecognized_tag == True
    model.close()

    # Open an existing asdf file

    model = DataModel(ASDF_FILE,
                     ignore_version_mismatch=False,
                     ignore_unrecognized_tag=True)

    assert model._asdf._ignore_version_mismatch == False
    assert model._asdf._ignore_unrecognized_tag == True
    model.close()

def test_open_asdf_model_s3(s3_root_dir):
    path = str(s3_root_dir.join("test.asdf"))
    with DataModel() as dm:
        dm.save(path)

    model = DataModel("s3://test-s3-data/test.asdf")
    assert isinstance(model, DataModel)

def test_open_fits_model_s3(s3_root_dir):
    path = str(s3_root_dir.join("test.fits"))
    with DataModel() as dm:
        dm.save(path)

    model = DataModel("s3://test-s3-data/test.fits")
    assert isinstance(model, DataModel)

def test_image_with_extra_keyword_to_multislit():
    with ImageModel((32, 32)) as im:
        im.save(TMP_FITS, overwrite=True)

    from astropy.io import fits
    with fits.open(TMP_FITS, mode="update", memmap=False) as hdulist:
        hdulist[1].header['BUNIT'] = 'x'

    with ImageModel(TMP_FITS) as im:
        with MultiSlitModel() as ms:
            ms.update(im)
            for i in range(3):
                ms.slits.append(ImageModel((4, 4)))
            assert len(ms.slits) == 3

            ms.save(TMP_FITS2, overwrite=True)

    with MultiSlitModel(TMP_FITS2) as ms:
        assert len(ms.slits) == 3
        for slit in ms.slits:
            assert slit.data.shape == (4, 4)


@pytest.fixture
def datamodel_for_update(tmpdir):
    """Provide ImageModel with one keyword each defined in PRIMARY and SCI
    extensions from the schema, and one each not in the schema"""
    tmpfile = str(tmpdir.join("old.fits"))
    with ImageModel((5, 5)) as im:
        # Add schema keywords, one to each extension
        im.meta.telescope = "JWST"
        im.meta.wcsinfo.crval1 = 5

        im.save(tmpfile)
    # Add non-schema keywords that will get dumped in the extra_fits attribute
    with fits.open(tmpfile, mode="update") as hdulist:
        hdulist["PRIMARY"].header["FOO"] = "BAR"
        hdulist["SCI"].header["BAZ"] = "BUZ"

    return tmpfile


@pytest.mark.parametrize("extra_fits", [True, False])
@pytest.mark.parametrize("only", [None, "PRIMARY", "SCI"])
def test_update_from_datamodel(tmpdir, datamodel_for_update, only, extra_fits):
    """Test update method does not update from extra_fits unless asked"""
    tmpfile = datamodel_for_update

    newtmpfile = str(tmpdir.join("new.fits"))
    with ImageModel((5, 5)) as newim:
        with ImageModel(tmpfile) as oldim:

            # Verify the fixture returns keywords we expect
            assert oldim.meta.telescope == "JWST"
            assert oldim.meta.wcsinfo.crval1 == 5
            assert oldim.extra_fits.PRIMARY.header == [['FOO', 'BAR', '']]
            assert oldim.extra_fits.SCI.header == [['BAZ', 'BUZ', '']]

            newim.update(oldim, only=only, extra_fits=extra_fits)
        newim.save(newtmpfile)

    with fits.open(newtmpfile) as hdulist:
        if extra_fits:
            if only == "PRIMARY":
                assert "TELESCOP" in hdulist["PRIMARY"].header
                assert "CRVAL1" not in hdulist["SCI"].header
                assert "FOO" in hdulist["PRIMARY"].header
                assert "BAZ" not in hdulist["SCI"].header
            elif only == "SCI":
                assert "TELESCOP" not in hdulist["PRIMARY"].header
                assert "CRVAL1" in hdulist["SCI"].header
                assert "FOO" not in hdulist["PRIMARY"].header
                assert "BAZ" in hdulist["SCI"].header
            else:
                assert "TELESCOP" in hdulist["PRIMARY"].header
                assert "CRVAL1" in hdulist["SCI"].header
                assert "FOO" in hdulist["PRIMARY"].header
                assert "BAZ" in hdulist["SCI"].header

        else:
            assert "FOO" not in hdulist["PRIMARY"].header
            assert "BAZ" not in hdulist["SCI"].header

            if only == "PRIMARY":
                assert "TELESCOP" in hdulist["PRIMARY"].header
                assert "CRVAL1" not in hdulist["SCI"].header
            elif only == "SCI":
                assert "TELESCOP" not in hdulist["PRIMARY"].header
                assert "CRVAL1" in hdulist["SCI"].header
            else:
                assert "TELESCOP" in hdulist["PRIMARY"].header
                assert "CRVAL1" in hdulist["SCI"].header


def test_update_from_dict(tmpdir):
    """Test update method from a dictionary"""
    tmpfile = str(tmpdir.join("update.fits"))
    with ImageModel((5, 5)) as im:
        update_dict = dict(meta=dict(telescope="JWST", wcsinfo=dict(crval1=5)))
        im.update(update_dict)
        im.save(tmpfile)

    with fits.open(tmpfile) as hdulist:
        assert "TELESCOP" in hdulist[0].header
        assert "CRVAL1" in hdulist[1].header


@pytest.fixture
def container():
    warnings.simplefilter("ignore")
    asn_file_path, asn_file_name = op.split(ASN_FILE)
    with pushdir(asn_file_path):
        with ModelContainer(asn_file_name) as c:
            for m in c:
                m.meta.observation.program_number = '0001'
                m.meta.observation.observation_number = '1'
                m.meta.observation.visit_number = '1'
                m.meta.observation.visit_group = '1'
                m.meta.observation.sequence_id = '01'
                m.meta.observation.activity_id = '1'
                m.meta.observation.exposure_number = '1'
                m.meta.instrument.name = 'NIRCAM'
                m.meta.instrument.channel = 'SHORT'
        yield c


def test_modelcontainer_iteration(container):
    for model in container:
        assert model.meta.telescope == 'JWST'


def test_modelcontainer_indexing(container):
    assert isinstance(container[0], DataModel)


def test_modelcontainer_group1(container):
    for group in container.models_grouped:
        assert len(group) == 2
        for model in group:
            pass


def test_modelcontainer_group2(container):
    container[0].meta.observation.exposure_number = '2'
    for group in container.models_grouped:
        assert len(group) == 1
        for model in group:
            pass
    container[0].meta.observation.exposure_number = '1'


def test_modelcontainer_group_names(container):
    assert len(container.group_names) == 1
    reset_group_id(container)
    container[0].meta.observation.exposure_number = '2'
    assert len(container.group_names) == 2


def test_object_node_iterator():
    im = ImageModel()
    items = []
    for i in im.meta.items():
        items.append(i[0])

    assert 'date' in items
    assert 'model_type' in items


def test_hasattr():
    model = DataModel()
    assert model.meta.hasattr('date')
    assert not model.meta.hasattr('filename')


def test_info():
    warnings.simplefilter("ignore")
    with datamodels.open(FITS_FILE) as model:
        info = model.info()
    matches = 0
    for line in info.split("\n"):
        words = line.split()
        if len(words) > 0:
            if words[0] == "data":
                matches += 1
                assert words[1] == "(32,40,35,5)", "Correct size for data"
                assert words[2] == "float32", "Correct type for data"
            elif words[0] == "dq":
                matches += 1
                assert words[1] == "(32,40,35,5)", "Correct size for dq"
                assert words[2] == "uint32", "Correct type for dq"
            elif words[0] == "err":
                matches += 1
                assert words[1] == "(32,40,35,5)", "Correct size for err"
                assert words[2] == "float32", "Correct type for err"
    assert matches== 3, "Check all extensions are described"


def test_validate_on_read():
    schema = ImageModel((10, 10))._schema.copy()
    schema['properties']['meta']['properties']['calibration_software_version']['fits_required'] = True

    with pytest.raises(jsonschema.ValidationError):
        with ImageModel(FITS_FILE, schema=schema, strict_validation=True):
            pass


def test_validate_required_field():
    im = ImageModel((10, 10), strict_validation=True)
    schema = im.meta._schema
    schema['properties']['telescope']['fits_required'] = True

    with pytest.raises(jsonschema.ValidationError):
        im.validate_required_fields()

    im.meta.telescope = 'JWST'
    im.validate_required_fields()


def test_multislit_model():
    data = np.arange(24, dtype=np.float32).reshape((6, 4))
    err = np.arange(24, dtype=np.float32).reshape((6, 4)) + 2
    wav = np.arange(24, dtype=np.float32).reshape((6, 4)) + 3
    dq = np.arange(24,dtype=np.uint32).reshape((6, 4)) + 1
    s0 = SlitDataModel(data=data, err=err, dq=dq, wavelength=wav)
    s1 = SlitDataModel(data=data+1, err=err+1, dq=dq+1, wavelength=wav+1)

    ms = MultiSlitModel()
    ms.slits.append(s0)
    ms.slits.append(s1)
    ms.meta.instrument.name = 'NIRSPEC'
    ms.meta.exposure.type = 'NRS_IMAGE'
    slit1 = ms[1]
    assert isinstance(slit1, SlitModel)
    assert slit1.meta.instrument.name == 'NIRSPEC'
    assert slit1.meta.exposure.type == 'NRS_IMAGE'
    assert_allclose(slit1.data, data + 1)


def test_slit_from_image():
    data = np.arange(24, dtype=np.float32).reshape((6, 4))
    im = ImageModel(data=data, err=data/2, dq=data)
    im.meta.instrument.name = "MIRI"
    slit_dm = SlitDataModel(im)
    assert_allclose(im.data, slit_dm.data)
    assert hasattr(slit_dm, 'wavelength')
    # this should be enabled after gwcs starts using non-coordinate inputs
    #assert not hasattr(slit_dm, 'meta')

    slit = SlitModel(im)
    assert_allclose(im.data, slit.data)
    assert_allclose(im.err, slit.err)
    assert hasattr(slit, 'wavelength')
    assert slit.meta.instrument.name == "MIRI"

    im = ImageModel(slit)
    assert type(im) == ImageModel

    im = ImageModel(slit_dm)
    assert type(im) == ImageModel


def test_ifuimage():
    data = np.arange(24, dtype=np.float32).reshape((6, 4))
    im = ImageModel(data=data, err=data/2, dq=data)
    ifuimage = IFUImageModel(im)
    assert_allclose(im.data, ifuimage.data)
    assert_allclose(im.err, ifuimage.err)
    assert_allclose(im.dq, ifuimage.dq)

    im = ImageModel(ifuimage)
    assert type(im) == ImageModel


def test_datamodel_raises_filenotfound():
    with pytest.raises(FileNotFoundError):
        DataModel(init='file_does_not_exist.fits')

def test_abvega_offset_model():
    filepath = os.path.join(ROOT_DIR, 'nircam_abvega_offset.asdf')
    model = ABVegaOffsetModel(filepath)
    assert isinstance(model, ABVegaOffsetModel)
    assert hasattr(model, 'abvega_offset')
    assert isinstance(model.abvega_offset, Table)
    assert model.abvega_offset.colnames == ['filter', 'pupil', 'abvega_offset']
    model.validate()
    model.close()

@pytest.fixture(scope='module')
def empty_model():
    """Create an empy model"""
    return DataModel()


@pytest.fixture(scope='module')
def jail_environ():
    """Lock cheanges to the environment"""
    original = copy(os.environ)
    try:
        yield
    finally:
        os.environ = original


@pytest.mark.parametrize(
    'case, default, expected', [
        (None, False, False),
        (None, True, True),
        ('0',  False, False),
        ('0',  True, False),
        ('1', False, True),
        ('1', True, True),
        ('2', False, True),
        ('2', True, True),
        ('f', False, False),
        ('f', True, False),
        ('F', False, False),
        ('F', True, False),
        ('false', False, False),
        ('false', True, False),
        ('False', False, False),
        ('False', True, False),
        ('FALSE', False, False),
        ('FALSE', True, False),
        ('t', False, True),
        ('t', True, True),
        ('T', False, True),
        ('T', True, True),
        ('true', False, True),
        ('true', True, True),
        ('True', False, True),
        ('True', True, True),
        ('TRUE', False, True),
        ('TRUE', True, True),
        ('y', False, True),
        ('y', True, True),
        ('Y', False, True),
        ('Y', True, True),
        ('yes', False, True),
        ('yes', True, True),
        ('Yes', False, True),
        ('Yes', True, True),
        ('YES', False, True),
        ('YES', True, True),
    ]
)
def test_get_envar_as_boolean(case, default, expected, jail_environ, empty_model):
    """Test various options to a boolean environmental variable"""
    var = '__TEST_GET_ENVAR_AS_BOOLEAN'
    if case is None:
        try:
            del os.environ[var]
        except KeyError:
            pass
    else:
        os.environ[var] = case

    value = empty_model._get_envar_as_boolean(var, default=default)
    assert value == expected


def test_getarray_noinit_valid():
    """Test for valid value return"""
    arr = np.ones((5, 5))
    model = ImageModel(data=arr)
    fetched = model.getarray_noinit('data')
    assert (fetched == arr).all()


def test_getarray_noinit_raises():
    """Test for error when accessing non-existent array"""
    arr = np.ones((5, 5))
    model = ImageModel(data=arr)
    with pytest.raises(AttributeError):
        fetched = model.getarray_noinit('area')


def test_getarray_noinit_noinit():
    """Test that calling on a non-existant array does not initialize that array"""
    arr = np.ones((5, 5))
    model = ImageModel(data=arr)
    try:
        fetched = model.getarray_noinit('area')
    except AttributeError:
        pass
    assert 'area' not in model.instance
