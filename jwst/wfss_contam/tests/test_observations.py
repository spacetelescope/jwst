import pytest
import numpy as np
import asdf
import os

from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from astropy.stats import sigma_clipped_stats
from photutils.datasets import make_100gaussians_image
from photutils.segmentation import SourceFinder

from jwst.wfss_contam.observations import background_subtract
from jwst.wfss_contam.disperse import dispersed_pixel
from jwst.wfss_contam.tests import data
from jwst.datamodels import SegmentationMapModel, ImageModel

data_path = os.path.split(os.path.abspath(data.__file__))[0]
DIR_IMAGE = "direct_image.fits"


@pytest.fixture(scope='module')
def direct_image():
    data = make_100gaussians_image()
    kernel = make_2dgaussian_kernel(3, size=5)
    data = convolve(data, kernel)
    return data


@pytest.fixture(scope='module')
def direct_image_with_gradient(tmp_cwd_module, direct_image):
    ny, nx = direct_image.shape
    y, x = np.mgrid[:ny, :nx]
    gradient = x * y / 5000.0
    data = direct_image + gradient

    # obs expects input list of direct image filenames
    model = ImageModel(data=data)
    model.save(DIR_IMAGE)

    return model


@pytest.fixture(scope='module')
def segmentation_map(direct_image):
    mean, median, stddev = sigma_clipped_stats(direct_image, sigma=3.0)
    threshold = median+3*stddev
    finder = SourceFinder(npixels=10)
    segm = finder(direct_image, threshold)

    # turn this into a jwst datamodel
    model = SegmentationMapModel(data=segm.data)
    asdf_file = asdf.open(os.path.join(data_path, "segmentation_wcs.asdf"))
    wcsobj = asdf_file.tree['wcs']
    model.meta.wcs = wcsobj

    return model


@pytest.fixture(scope='module')
def grism_wcs():
    asdf_file = asdf.open(os.path.join(data_path, "grism_wcs.asdf"))
    wcsobj = asdf_file.tree['wcs']
    return wcsobj


def test_background_subtract(direct_image_with_gradient):
    data = direct_image_with_gradient.data
    subtracted_data = background_subtract(data)
    mean, median, stddev = sigma_clipped_stats(subtracted_data, sigma=3.0)
    assert np.isclose(mean, 0.0, atol=0.2*stddev)


def test_disperse_oversample_same_result(grism_wcs, segmentation_map):
    '''
    Coverage for bug where wavelength oversampling led to double-counted fluxes

    note: segmentation_map fixture needs to be able to find module-scoped direct_image
    fixture, so it must be imported here
    '''

    # manual input of input params set the same as test_observations.py
    x0 = 300.5
    y0 = 300.5
    order = 1
    width = 1.0
    height = 1.0
    lams = [2.0]
    flxs = [1.0]
    ID = 0
    naxis = (300, 500)
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)
    seg_wcs = segmentation_map.meta.wcs
    0, (300, 500), 2, False, 
    xoffset = 2200
    yoffset = 1000


    xs, ys, areas, lams_out, counts_1, ID = dispersed_pixel(
                    x0, y0, width, height, lams, flxs, order, wmin, wmax,
                    sens_waves, sens_resp, seg_wcs, grism_wcs, ID, naxis,
                    oversample_factor=1, extrapolate_sed=False, xoffset=xoffset,
                    yoffset=yoffset)
    
    xs, ys, areas, lams_out, counts_3, ID = dispersed_pixel(
                x0, y0, width, height, lams, flxs, order, wmin, wmax,
                sens_waves, sens_resp, seg_wcs, grism_wcs, ID, naxis,
                oversample_factor=3, extrapolate_sed=False, xoffset=xoffset,
                yoffset=yoffset)

    assert np.isclose(np.sum(counts_1), np.sum(counts_3), rtol=1/sens_waves.size)