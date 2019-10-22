"""
Test datamodel.open
"""

import os
import os.path
import warnings

import pytest
import numpy as np
from astropy.io import fits

from ..util import open
from .. import (DataModel, ModelContainer, ImageModel, ReferenceFileModel,
                ReferenceImageModel, ReferenceCubeModel, ReferenceQuadModel,
                FlatModel, MaskModel, NircamPhotomModel, GainModel,
                ReadnoiseModel, DistortionModel)
from jwst import datamodels


def test_open_fits():
    """Test opening a model from a FITS file"""

    warnings.simplefilter("ignore")
    fits_file = t_path('test.fits')
    m = open(fits_file)
    assert isinstance(m, DataModel)

def test_open_fits_s3():
    """Test opening a model from a FITS file on S3"""
    m = open("s3://test-s3-data/data_model.fits")
    assert isinstance(m, DataModel)

def test_open_asdf_s3():
    """Test opening a model from an ASDF file on S3"""
    m = open("s3://test-s3-data/data_model.asdf")
    assert isinstance(m, DataModel)

def test_open_association():
    """Test for opening an association"""

    asn_file = t_path('association.json')
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "model_type not found")
        m = open(asn_file)
    assert isinstance(m, ModelContainer)

def test_open_shape():
    init = (200, 200)
    model = open(init)
    assert type(model) == ImageModel
    model.close()

def test_open_illegal():
    with pytest.raises(ValueError):
        init = 5
        open(init)

def test_open_hdulist():
    hdulist = fits.HDUList()
    data = np.empty((50, 50), dtype=np.float32)
    primary = fits.PrimaryHDU()
    hdulist.append(primary)
    science = fits.ImageHDU(data=data, name='SCI')
    hdulist.append(science)

    model = open(hdulist)
    assert type(model) == ImageModel
    model.close()

def test_open_image():
    warnings.simplefilter("ignore")
    image_name = t_path('jwst_image.fits')
    model = open(image_name)
    assert type(model) == ImageModel
    model.close()


def test_open_reference_files():
    files = {'nircam_flat.fits' : FlatModel,
             'nircam_mask.fits' : MaskModel,
             'nircam_photom.fits' : NircamPhotomModel,
             'nircam_gain.fits' : GainModel,
             'nircam_readnoise.fits' : ReadnoiseModel}

    warnings.simplefilter("ignore")
    for base_name, klass in files.items():
        file = t_path(base_name)
        model = open(file)
        if model.shape:
            ndim = len(model.shape)
        else:
            ndim = 0

        if ndim == 0:
            my_klass = ReferenceFileModel
        elif ndim == 2:
            my_klass = ReferenceImageModel
        elif ndim == 3:
            my_klass = ReferenceCubeModel
        elif ndim == 4:
            my_klass = ReferenceQuadModel
        else:
            my_klass = None

        assert isinstance(model, my_klass)
        model.close()

        model = klass(file)
        assert isinstance(model, klass)
        model.close()

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
