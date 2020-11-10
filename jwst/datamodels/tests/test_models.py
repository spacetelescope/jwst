import os
import warnings

from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from numpy.testing import assert_allclose
import numpy as np
import pytest

from jwst.datamodels import (JwstDataModel, ImageModel, QuadModel,
                             MultiSlitModel, ModelContainer, SlitModel,
                             SlitDataModel, IFUImageModel, ABVegaOffsetModel)
from jwst import datamodels
from jwst.lib.file_utils import pushdir


ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')
FITS_FILE = os.path.join(ROOT_DIR, 'test.fits')
ASN_FILE = os.path.join(ROOT_DIR, 'association.json')


@pytest.fixture
def jail_environ():
    """Lock changes to the environment"""
    original = os.environ.copy()
    try:
        yield
    finally:
        os.environ = original


@pytest.fixture
def make_models(tmpdir_factory):
    """Create basic models

    Returns
    -------
    path_just_fits, path_model : (str, str)
        `path_just_fits` is a FITS file of `JwstDataModel` without the ASDF extension.
        `path_model` is a FITS file of `JwstDataModel` with the ASDF extension.
    """
    path = tmpdir_factory.mktemp('skip_fits_update')
    path_just_fits = str(path / 'just_fits.fits')
    path_model = str(path / 'model.fits')
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['EXP_TYPE'] = 'NRC_IMAGE'
    primary_hdu.header['DATAMODL'] = "JwstDataModel"
    hduls = fits.HDUList([primary_hdu])
    hduls.writeto(path_just_fits)
    model = JwstDataModel(hduls)
    model.save(path_model)
    return {
        'just_fits': path_just_fits,
        'model': path_model
    }


@pytest.mark.parametrize(
    'which_file, skip_fits_update, expected_exp_type',
    [
        ('just_fits', None,  'FGS_DARK'),
        ('just_fits', False, 'FGS_DARK'),
        ('just_fits', True,  'FGS_DARK'),
        ('model',     None,  'FGS_DARK'),
        ('model',     False, 'FGS_DARK'),
        ('model',     True,  'NRC_IMAGE')
    ]
)
@pytest.mark.parametrize(
    'open_func',
    [JwstDataModel, datamodels.open]
)
@pytest.mark.parametrize(
    'use_env',
    [False, True]
)
def test_skip_fits_update(jail_environ,
                          use_env,
                          make_models,
                          open_func,
                          which_file,
                          skip_fits_update,
                          expected_exp_type):
    """Test skip_fits_update setting"""
    # Setup the FITS file, modifying a header value
    path = make_models[which_file]
    hduls = fits.open(path)
    hduls[0].header['exp_type'] = 'FGS_DARK'

    # Decide how to skip. If using the environmental,
    # set that and pass None to the open function.
    try:
        del os.environ['SKIP_FITS_UPDATE']
    except KeyError:
        # No need to worry, environmental doesn't exist anyways
        pass
    if use_env:
        if skip_fits_update is not None:
            os.environ['SKIP_FITS_UPDATE'] = str(skip_fits_update)
            skip_fits_update = None

    model = open_func(hduls, skip_fits_update=skip_fits_update)
    assert model.meta.exposure.type == expected_exp_type


def test_subarray():
    with JwstDataModel(FITS_FILE) as dm:
        dm.meta.subarray.xstart



def test_open():
    with datamodels.open():
        pass

    with datamodels.open((50, 50)):
        pass

    with pytest.warns(datamodels.util.NoTypeWarning):
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


def test_imagemodel():
    dims = (10, 10)
    with ImageModel(dims) as dm:
        assert dm.data.shape == dims
        assert dm.err.shape == dims
        assert dm.dq.shape == dims
        assert dm.data.mean() == 0.0
        assert dm.err.mean() == 0.0
        assert dm.dq.mean() == 0.0
        assert dm.zeroframe.shape == dims
        assert dm.area.shape == dims
        assert dm.pathloss_point.shape == dims
        assert dm.pathloss_uniform.shape == dims


def test_multislit():
    with MultiSlitModel() as dm:
        dm.slits.append(dm.slits.item())
        slit = dm.slits[-1]
        slit.data = np.random.rand(5, 5)
        slit.dm = np.random.rand(5, 5)
        slit.err = np.random.rand(5, 5)
        assert slit.wavelength.shape == (0, 0)
        assert slit.pathloss_point.shape == (0, 0)
        assert slit.pathloss_uniform.shape == (0, 0)
        assert slit.barshadow.shape == (0, 0)


def test_multislit_from_image():
    with ImageModel((64, 64)) as im:
        with MultiSlitModel(im) as ms:
            assert len(ms.slits) == 1
            assert ms.slits[0].data.shape == (64, 64)


def test_multislit_from_saved_image(tmp_path):
    path = tmp_path / "multislit_from_image.fits"
    with ImageModel((64, 64)) as im:
        im.save(path)

    with MultiSlitModel(path) as ms:
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
        for key, val in ms.items():
            assert isinstance(val, (bytes, str, int, float, bool, Time))


def test_multislit_copy(tmp_path):
    path = tmp_path / "multislit.fits"
    with MultiSlitModel() as input:
        for i in range(4):
            input.slits.append(input.slits.item(data=np.empty((50, 50))))

        assert len(input.slits) == 4
        input.save(path)

        output = input.copy()
        assert len(output.slits) == 4

    from astropy.io import fits
    with fits.open(path, memmap=False) as hdulist:
        assert len(hdulist) == 6

    with MultiSlitModel(path) as model:
        for i, slit in enumerate(model.slits):
            pass
        assert i+1 == 4

        output = model.copy()
        assert len(output.slits) == 4


def test_model_with_nonstandard_primary_array():
    class NonstandardPrimaryArrayModel(JwstDataModel):
        schema_url = os.path.join(ROOT_DIR, "nonstandard_primary_array.schema.yaml")

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


def test_image_with_extra_keyword_to_multislit(tmp_path):
    path = tmp_path / "extra_keyword_image.fits"
    path2 = tmp_path / "extra_keyword_multislit.fits"
    with ImageModel((32, 32)) as im:
        im.save(path)

    from astropy.io import fits
    with fits.open(path, mode="update", memmap=False) as hdulist:
        hdulist[1].header['BUNIT'] = 'x'

    with ImageModel(path) as im:
        with MultiSlitModel() as ms:
            ms.update(im)
            for i in range(3):
                ms.slits.append(ImageModel((4, 4)))
            assert len(ms.slits) == 3

            ms.save(path2)

    with MultiSlitModel(path2) as ms:
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
    asn_file_path, asn_file_name = os.path.split(ASN_FILE)
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


def reset_group_id(container):
    """Remove group_id from all models in container"""
    for m in container:
        try:
            del m.meta.group_id
        except AttributeError:
            pass


def test_modelcontainer_iteration(container):
    for model in container:
        assert model.meta.telescope == 'JWST'


def test_modelcontainer_indexing(container):
    assert isinstance(container[0], JwstDataModel)


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


def test_modelcontainer_error_from_asn(tmpdir):
    from jwst.associations.asn_from_list import asn_from_list

    asn = asn_from_list(["foo.fits"], product_name="foo_out")
    name, serialized = asn.dump(format="json")
    path = str(tmpdir.join(name))
    with open(path, "w") as f:
        f.write(serialized)

    # The foo.fits file doesn't exist
    with pytest.raises(FileNotFoundError):
        datamodels.open(path)


def test_multislit_model():
    data = np.arange(24, dtype=np.float32).reshape((6, 4))
    err = np.arange(24, dtype=np.float32).reshape((6, 4)) + 2
    wav = np.arange(24, dtype=np.float32).reshape((6, 4)) + 3
    dq = np.arange(24, dtype=np.uint32).reshape((6, 4)) + 1
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
    # assert not hasattr(slit_dm, 'meta')

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


def test_abvega_offset_model():
    path = os.path.join(ROOT_DIR, 'nircam_abvega_offset.asdf')
    with ABVegaOffsetModel(path) as model:
        assert isinstance(model, ABVegaOffsetModel)
        assert hasattr(model, 'abvega_offset')
        assert isinstance(model.abvega_offset, Table)
        assert model.abvega_offset.colnames == ['filter', 'pupil', 'abvega_offset']
        model.validate()
