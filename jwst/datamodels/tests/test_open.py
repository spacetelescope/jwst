"""
Test datamodel.open
"""

import os
import os.path
from pathlib import Path, PurePath
import warnings

import pytest
import numpy as np
from astropy.io import fits
from stdatamodels import DataModel
from stdatamodels.validate import ValidationError, ValidationWarning

from jwst.datamodels import (JwstDataModel, ModelContainer, ImageModel,
                             RampModel, CubeModel, ReferenceFileModel, ReferenceImageModel,
                             ReferenceCubeModel, ReferenceQuadModel)
from jwst import datamodels
from jwst.datamodels import util

import asdf

# Define artificial memory size
MEMORY = 100  # 100 bytes


@pytest.mark.parametrize('guess', [True, False])
def test_guess(guess):
    """Test the guess parameter to the open func"""
    path = Path(t_path('test.fits'))

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        if guess is None or guess:
            # Default is to guess at the model type.
            with datamodels.open(path, guess=guess) as model:
                assert isinstance(model, JwstDataModel)
        else:
            # Without guessing, the call should fail.
            with pytest.raises(TypeError):
                with datamodels.open(path, guess=guess) as model:
                    pass


def test_mirirampmodel_deprecation(tmp_path):
    """Test that a deprecated MIRIRampModel can be opened"""
    path = str(tmp_path / "ramp.fits")
    # Create a MIRIRampModel, working around the deprecation.
    model = datamodels.RampModel((1, 1, 10, 10))
    model.save(path)
    hduls = fits.open(path, mode='update')
    hduls[0].header['datamodl'] = 'MIRIRampModel'
    hduls.close()

    # Test it.
    with pytest.warns(DeprecationWarning):
        miri_ramp = datamodels.open(path)
    assert isinstance(miri_ramp, datamodels.RampModel)


@pytest.fixture
def mock_get_available_memory(monkeypatch):
    def mock(include_swap=True):
        available = MEMORY
        if include_swap:
            available *= 2
        return available
    monkeypatch.setattr(util, 'get_available_memory', mock)


@pytest.mark.parametrize(
    'allowed_env, allowed_explicit, result',
    [
        (None, None, True),  # Perform no check.
        (0.1, None, False),  # Force too little memory.
        (0.1, 1.0, True),    # Explicit overrides environment.
        (1.0, 0.1, False),   # Explicit overrides environment.
        (None, 0.1, False),  # Explicit overrides environment.
    ]
)
def test_check_memory_allocation_env(monkeypatch, mock_get_available_memory,
                                     allowed_env, allowed_explicit, result):
    """Check environmental control over memory check"""
    if allowed_env is None:
        monkeypatch.delenv('DMODEL_ALLOWED_MEMORY', raising=False)
    else:
        monkeypatch.setenv('DMODEL_ALLOWED_MEMORY', str(allowed_env))

    # Allocate amount that would fit at 100% + swap.
    can_allocate, required = util.check_memory_allocation(
        (MEMORY // 2, 1), allowed=allowed_explicit,
    )
    assert can_allocate is result


@pytest.mark.parametrize(
    'dim, allowed, include_swap, result',
    [
        (MEMORY // 2, 1.0, True, True),    # Fit within memory and swap
        (MEMORY // 2, 1.0, False, False),  # Does not fit without swap
        (MEMORY, 1.0, True, False),        # Does not fit at all
        (MEMORY, None, True, True),        # Check disabled
        (MEMORY // 2, 0.1, True, False),   # Does not fit in restricted memory
    ]
)
def test_check_memory_allocation(mock_get_available_memory, dim, allowed, include_swap, result):
    """Check general operation of check_memory_allocation"""
    can_allocate, required = util.check_memory_allocation(
        (dim, 1), allowed=allowed, include_swap=include_swap
    )
    assert can_allocate is result


def test_open_from_pathlib():
    """Test opening a PurePath object"""
    path = Path(t_path('test.fits'))
    assert isinstance(path, PurePath)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        with datamodels.open(path) as model:
            assert isinstance(model, JwstDataModel)


def test_open_fits():
    """Test opening a model from a FITS file"""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        fits_file = t_path('test.fits')
        with datamodels.open(fits_file) as model:
            assert isinstance(model, JwstDataModel)


def test_open_fits_s3(s3_root_dir):
    """Test opening a model from a FITS file on S3"""
    path = str(s3_root_dir.join("test.fits"))
    with JwstDataModel() as dm:
        dm.save(path)

    with datamodels.open("s3://test-s3-data/test.fits") as m:
        assert isinstance(m, JwstDataModel)


def test_open_asdf_s3(s3_root_dir):
    """Test opening a model from an ASDF file on S3"""
    path = str(s3_root_dir.join("test.asdf"))
    with JwstDataModel() as dm:
        dm.save(path)

    with datamodels.open("s3://test-s3-data/test.asdf") as m:
        assert isinstance(m, JwstDataModel)


def test_open_association():
    """Test for opening an association"""

    asn_file = t_path('association.json')
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        with datamodels.open(asn_file) as c:
            assert isinstance(c, ModelContainer)
            for model in c:
                assert model.meta.asn.table_name == "association.json"
                assert model.meta.asn.pool_name == "pool"


def test_container_open_asn_with_sourcecat():
    path = t_path("association_w_cat.json")
    with datamodels.open(path, asn_exptypes="science") as c:
        for model in c:
            assert model.meta.asn.table_name == "association_w_cat.json"


def test_open_none():
    with datamodels.open() as model:
        assert isinstance(model, JwstDataModel)


def test_open_shape():
    shape = (50, 20)
    with datamodels.open(shape) as model:
        assert isinstance(model, ImageModel)
        assert model.shape == shape


def test_open_illegal():
    with pytest.raises(ValueError):
        init = 5
        datamodels.open(init)


def test_open_hdulist(tmp_path):
    hdulist = fits.HDUList()
    primary = fits.PrimaryHDU()
    data = np.empty((50, 50), dtype=np.float32)
    science = fits.ImageHDU(data=data, name='SCI', ver=1)
    hdulist.append(primary)
    hdulist.append(science)

    # datamodels.open() can't open pathlib objects
    path = str(tmp_path / "jwst_image.fits")
    hdulist.writeto(path)

    with datamodels.open(hdulist) as model:
        assert isinstance(model, ImageModel)

    with pytest.warns(datamodels.util.NoTypeWarning) as record:
        with datamodels.open(path) as model:
            assert isinstance(model, ImageModel)
            assert len(record) == 1
            assert "model_type not found" in record[0].message.args[0]


def test_open_ramp(tmp_path):
    """Open 4D data without a DQ as RampModel"""
    path = str(tmp_path / "ramp.fits")
    shape = (2, 3, 4, 5)
    with fits.HDUList(fits.PrimaryHDU()) as hdulist:
        hdulist.append(fits.ImageHDU(data=np.zeros(shape), name="SCI", ver=1))
        hdulist.writeto(path)

    with pytest.warns(datamodels.util.NoTypeWarning):
        with datamodels.open(path) as model:
            assert isinstance(model, RampModel)


def test_open_cube(tmp_path):
    """Open 3D data as CubeModel"""
    path = str(tmp_path / "ramp.fits")
    shape = (2, 3, 4)
    with fits.HDUList(fits.PrimaryHDU()) as hdulist:
        hdulist.append(fits.ImageHDU(data=np.zeros(shape), name="SCI", ver=1))
        hdulist.writeto(path)

    with pytest.warns(datamodels.util.NoTypeWarning):
        with datamodels.open(path) as model:
            assert isinstance(model, CubeModel)


@pytest.mark.parametrize("model_class, shape", [
    (ReferenceFileModel, None),
    (ReferenceImageModel, (10, 10)),
    (ReferenceCubeModel, (3, 3, 3)),
    (ReferenceQuadModel, (2, 2, 2, 2)),
])
def test_open_reffiles(tmp_path, model_class, shape):
    """Try opening files with a REFTYPE keyword and different data/dq shapes"""
    path = str(tmp_path / "reffile.fits")
    with fits.HDUList(fits.PrimaryHDU()) as hdulist:
        hdulist["PRIMARY"].header.append(("REFTYPE", "foo"))
        if shape is not None:
            hdulist.append(fits.ImageHDU(data=np.zeros(shape), name="SCI", ver=1))
            hdulist.append(fits.ImageHDU(data=np.zeros(shape, dtype=np.uint), name="DQ", ver=1))
        hdulist.writeto(path)

    with pytest.warns(datamodels.util.NoTypeWarning):
        with datamodels.open(path) as model:
            assert isinstance(model, model_class)


@pytest.mark.parametrize("suffix", [".asdf", ".fits"])
def test_open_readonly(tmp_path, suffix):
    """Test opening a FITS-format datamodel that is read-only on disk"""
    path = str(tmp_path / f"readonly{suffix}")

    with ImageModel(data=np.zeros((10, 10))) as model:
        model.meta.telescope = 'JWST'
        model.meta.instrument.name = 'NIRCAM'
        model.meta.instrument.detector = 'NRCA4'
        model.meta.instrument.channel = 'SHORT'
        model.save(path)

    os.chmod(path, 0o440)
    assert os.access(path, os.W_OK) is False

    with datamodels.open(path) as model:
        assert model.meta.telescope == 'JWST'
        assert isinstance(model, ImageModel)


# Utilities
def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(test_dir, partial_path)


@pytest.mark.parametrize("suffix", ["asdf", "fits"])
def test_open_asdf_no_datamodel_class(tmp_path, suffix):
    path = str(tmp_path / f"no_model.{suffix}")
    model = DataModel()
    model.save(path)

    # Note: only the fits open emits a "model_type not found" warning.  Both
    # fits and asdf should behave the same
    with datamodels.open(path) as m:
        assert isinstance(m, DataModel)


def test_open_asdf(tmp_path):
    path = str(tmp_path / "straight_asdf.asdf")
    tree = {"foo": 42, "bar": 13, "seq": np.arange(100)}
    with asdf.AsdfFile(tree) as af:
        af.write_to(path)

    with datamodels.open(path) as m:
        assert isinstance(m, DataModel)


def test_open_kwargs_asdf(tmp_path):
    """
    Test that unrecognized kwargs to the datamodels.open function
    are passed on to the model class constructor.
    """
    file_path = tmp_path / "test.asdf"

    with pytest.warns(ValidationWarning):
        model = ImageModel((4, 4), pass_invalid_values=True)
        model.meta.instrument.name = "CELESTRON"
        model.save(file_path)

    with pytest.raises(ValidationError):
        with datamodels.open(file_path, strict_validation=True) as model:
            model.validate()


def test_open_kwargs_fits(tmp_path):
    """
    Test that unrecognized kwargs to the datamodels.open function
    are passed on to the model class constructor.  Similar to the
    above, except the invalid file must be created differently
    because DataModel can't save an invalid .fits file.
    """
    file_path = tmp_path / "test.fits"

    model = ImageModel((4, 4))
    model.save(file_path)

    with fits.open(file_path, mode="update") as hdul:
        hdul[0].header["INSTRUME"] = "CELESTRON"
        hdul.flush()

    with pytest.raises(ValidationError):
        with datamodels.open(file_path, strict_validation=True) as model:
            model.validate()
