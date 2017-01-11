from __future__ import absolute_import, unicode_literals, division, print_function

"""
Test datamodel.open
"""

import os
import os.path
import numpy as np
from astropy.io import fits

from .. import (DataModel, ModelContainer)
from ..util import open

from jwst import datamodels
import jwst.datamodels.image
import jwst.datamodels.cube
import jwst.datamodels.reference

from jwst.datamodels.reference import ReferenceFileModel, ReferenceImageModel, \
                                ReferenceCubeModel, ReferenceQuadModel 
from jwst.datamodels.flat import FlatModel
from jwst.datamodels.mask import MaskModel
from jwst.datamodels.photom import NircamPhotomModel
from jwst.datamodels.gain import GainModel
from jwst.datamodels.readnoise import ReadnoiseModel

ROOT_DIR = os.path.join(os.path.dirname(__file__), 'data')

def test_open_fits():
    """Test opening a model from a FITS file"""

    fits_file = t_path('test.fits')
    m = open(fits_file)
    assert isinstance(m, DataModel)


def test_open_association():
    """Test for opening an association"""

    asn_file = t_path('association.json')
    m = open(asn_file)
    assert isinstance(m, ModelContainer)

def test_open_shape():
    init = (200, 200)
    model = jwst.datamodels.open(init)
    assert type(model) == jwst.datamodels.image.ImageModel
    model.close()

def test_open_illegal():
    init = 5
    try:
        model = jwst.datamodels.open(init)
    except ValueError:
        fail = 1
    else:
        fail = 0
    assert fail == 1

def test_open_hdulist():
    hdulist = fits.HDUList()
    data = np.empty((50, 50), dtype=np.float32)
    primary = fits.PrimaryHDU()
    hdulist.append(primary)
    science = fits.ImageHDU(data=data, name='SCI')
    hdulist.append(science)

    model = jwst.datamodels.open(hdulist)
    assert type(model) == jwst.datamodels.image.ImageModel
    model.close()

def test_open_image():
    image_name = t_path('jwst_image.fits')
    model = jwst.datamodels.open(image_name)
    assert type(model) == jwst.datamodels.image.ImageModel
    model.close()


def test_open_reference_files():
    files = {'nircam_flat.fits' : FlatModel,
             'nircam_mask.fits' : MaskModel,
             'nircam_photom.fits' : NircamPhotomModel,
             'nircam_gain.fits' : GainModel,
             'nircam_readnoise.fits' : ReadnoiseModel}
    
    for base_name, klass in files.items():
        file = t_path(base_name)
        model = jwst.datamodels.open(file)
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
            my_class = None
            
        assert isinstance(model, my_klass)
        model.close()
        
        model = klass(file)
        assert isinstance(model, klass)
        model.close()

# Utilities
def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(test_dir, partial_path)
