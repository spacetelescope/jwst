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

from jwst.datamodels import (DataModel, ModelContainer, ImageModel,
    DistortionModel, RampModel, CubeModel, ReferenceFileModel, ReferenceImageModel,
    ReferenceCubeModel, ReferenceQuadModel)
from jwst import datamodels
from jwst.datamodels import util

# Define artificial memory size
MEMORY = 100  # 100 bytes


@pytest.fixture
def mock_get_available_memory(monkeypatch):
    def mock(include_swap=True):
        avaliable = MEMORY
        if include_swap:
            avaliable *= 2
        return avaliable
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
        monkeypatch.setenv('DMODEL_ALLOWED_MEMORY', allowed_env)

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
            assert isinstance(model, DataModel)


def test_open_fits():
    """Test opening a model from a FITS file"""

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        fits_file = t_path('test.fits')
        with datamodels.open(fits_file) as model:
            assert isinstance(model, DataModel)


def test_open_fits_s3(s3_root_dir):
    """Test opening a model from a FITS file on S3"""
    path = str(s3_root_dir.join("test.fits"))
    with DataModel() as dm:
        dm.save(path)

    with datamodels.open("s3://test-s3-data/test.fits") as m:
        assert isinstance(m, DataModel)


def test_open_asdf_s3(s3_root_dir):
    """Test opening a model from an ASDF file on S3"""
    path = str(s3_root_dir.join("test.asdf"))
    with DataModel() as dm:
        dm.save(path)

    with datamodels.open("s3://test-s3-data/test.asdf") as m:
        assert isinstance(m, DataModel)


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


def test_open_shape():
    init = (200, 200)
    with datamodels.open(init) as model:
        assert type(model) == ImageModel


def test_open_illegal():
    with pytest.raises(ValueError):
        init = 5
        datamodels.open(init)


def test_open_hdulist():
    hdulist = fits.HDUList()
    data = np.empty((50, 50), dtype=np.float32)
    primary = fits.PrimaryHDU()
    hdulist.append(primary)
    science = fits.ImageHDU(data=data, name='SCI')
    hdulist.append(science)

    with datamodels.open(hdulist) as model:
        assert type(model) == ImageModel


def test_open_image():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        image_name = t_path('jwst_image.fits')
        with datamodels.open(image_name) as model:
            assert type(model) == ImageModel


def test_open_ramp(tmpdir):
    """Open 4D data without a DQ as RampModel"""
    path = str(tmpdir.join("ramp.fits"))
    shape = (2, 3, 4, 5)
    with fits.HDUList(fits.PrimaryHDU()) as hdulist:
        hdulist.append(fits.ImageHDU(data=np.zeros(shape), name="SCI", ver=1))
        hdulist.writeto(path)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        with datamodels.open(path) as model:
            assert isinstance(model, RampModel)


def test_open_cube(tmpdir):
    """Open 3D data as CubeModel"""
    path = str(tmpdir.join("ramp.fits"))
    shape = (2, 3, 4)
    with fits.HDUList(fits.PrimaryHDU()) as hdulist:
        hdulist.append(fits.ImageHDU(data=np.zeros(shape), name="SCI", ver=1))
        hdulist.writeto(path)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        with datamodels.open(path) as model:
            assert isinstance(model, CubeModel)


@pytest.mark.parametrize("model_class, shape", [
    (ReferenceFileModel, None),
    (ReferenceImageModel, (10, 10)),
    (ReferenceCubeModel, (3, 3, 3)),
    (ReferenceQuadModel, (2, 2, 2, 2)),
])
def test_open_reffiles(tmpdir, model_class, shape):
    """Try opening files with a REFTYPE keyword and different data/dq shapes"""
    path = str(tmpdir.join("reffile.fits"))
    with fits.HDUList(fits.PrimaryHDU()) as hdulist:
        hdulist["PRIMARY"].header.append(("REFTYPE", "foo"))
        if shape is not None:
            hdulist.append(fits.ImageHDU(data=np.zeros(shape), name="SCI", ver=1))
            hdulist.append(fits.ImageHDU(data=np.zeros(shape, dtype=np.uint), name="DQ", ver=1))
        hdulist.writeto(path)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        with datamodels.open(path) as model:
            assert isinstance(model, model_class)


def test_open_fits_readonly(tmpdir):
    """Test opening a FITS-format datamodel that is read-only on disk"""
    tmpfile = str(tmpdir.join('readonly.fits'))
    data = np.arange(100, dtype=np.float).reshape(10, 10)

    with ImageModel(data=data) as model:
        model.meta.telescope = 'JWST'
        model.meta.instrument.name = 'NIRCAM'
        model.meta.instrument.detector = 'NRCA4'
        model.meta.instrument.channel = 'SHORT'
        model.save(tmpfile)

    os.chmod(tmpfile, 0o440)
    assert os.access(tmpfile, os.W_OK) == False

    with datamodels.open(tmpfile) as model:
        assert model.meta.telescope == 'JWST'


def test_open_asdf_readonly(tmpdir):
    tmpfile = str(tmpdir.join('readonly.asdf'))

    with DistortionModel() as model:
        model.meta.telescope = 'JWST'
        model.meta.instrument.name = 'NIRCAM'
        model.meta.instrument.detector = 'NRCA4'
        model.meta.instrument.channel = 'SHORT'
        model.save(tmpfile)

    os.chmod(tmpfile, 0o440)
    assert os.access(tmpfile, os.W_OK) == False

    with datamodels.open(tmpfile) as model:
        assert model.meta.telescope == 'JWST'

# Utilities
def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(test_dir, partial_path)
