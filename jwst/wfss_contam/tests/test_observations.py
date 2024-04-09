import pytest
import numpy as np
import asdf
import os

from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from astropy.stats import sigma_clipped_stats
from photutils.datasets import make_100gaussians_image
from photutils.segmentation import SourceFinder

from jwst.wfss_contam.observations import background_subtract, _select_ids, Observation
from jwst.wfss_contam.tests import data
from jwst.datamodels import SegmentationMapModel, ImageModel

data_path = os.path.split(os.path.abspath(data.__file__))[0]
DIR_IMAGE = "direct_image.fits"

#@pytest.fixture(scope='module')
#def create_source_catalog():
#    '''Mock source catalog'''
#    pass

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


@pytest.fixture(scope='module')
def observation(direct_image_with_gradient, segmentation_map, grism_wcs):
    '''
    set up observation object with mock data.
    direct_image_with_gradient still needs to be run to produce the file,
    even though it is not called directly
    '''
    filter_name = "F200W"
    seg = segmentation_map.data
    all_IDs = np.array(list(set(np.ravel(seg))))
    IDs = all_IDs[50:52]
    obs = Observation([DIR_IMAGE], segmentation_map, grism_wcs, filter_name, ID=IDs,
                 sed_file=None, extrapolate_sed=False,
                 boundaries=[], offsets=[0, 0], renormalize=True, max_cpu=1)
    return obs


@pytest.mark.parametrize("ID, expected", [(None, [1,2,3]), (2, [2]), (np.array([1,3]), [1,3])])
def test_select_ids(ID, expected):
    all_ids = [1, 2, 3]
    assert _select_ids(ID, all_ids) == expected
    assert isinstance(_select_ids(ID, all_ids), list)


def test_background_subtract(direct_image_with_gradient):

    data = direct_image_with_gradient.data
    subtracted_data = background_subtract(data)
    mean, median, stddev = sigma_clipped_stats(subtracted_data, sigma=3.0)
    assert np.isclose(mean, 0.0, atol=0.2*stddev)


def test_create_pixel_list(observation, segmentation_map):
    '''
    create_pixel_list is called on initialization
    compare the pixel lists to values determined directly
    from segmentation map

    Note: still need test coverage for flux dictionary
    '''
    seg = segmentation_map.data
    all_IDs = np.array(list(set(np.ravel(seg))))
    IDs = all_IDs[50:52]
    for i, ID in enumerate(IDs):
        pixels_y, pixels_x = np.where(seg == ID)
        assert np.all(np.isin(observation.xs[i], pixels_x))
        assert np.all(np.isin(observation.ys[i], pixels_y))
        assert len(observation.fluxes[2.0][i]) == pixels_x.size


def test_disperse_chunk(observation):
    '''
    disperse_chunk is a static method so need to give it lots of observation attributes as input

    Note: it's not obvious how to get a trivial flux example from first principles
    even setting all input fluxes in dict to 1, because transforms change
    pixel areas in nontrivial ways. seems a bad idea to write a test that 
    asserts the answer as it is currently, in case step gets updated slightly
    in the future
    '''
    obs = observation
    i = 1
    order = 1
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)

    # manually change x,y offset because took transform from a real direct image, with different
    # pixel 0,0 than the mock data. This puts i=1, order 1 onto the real grism image
    obs.xoffset = 2200
    obs.yoffset = 1000

    # set all fluxes to unity to try to make a trivial example
    obs.fluxes[2.0][i] = np.ones(obs.fluxes[2.0][i].shape)

    disperse_chunk_args = [i, order, wmin, wmax, sens_waves, sens_resp,
                       obs.IDs[i], obs.xs[i], obs.ys[i], 
                       obs.fluxes, 
                       obs.seg_wcs, obs.grism_wcs, obs.dims, 
                       obs.extrapolate_sed, obs.xoffset, obs.yoffset]

    (chunk, chunk_bounds, sid, order_out) = obs.disperse_chunk(*disperse_chunk_args)

    #trivial bookkeeping
    assert sid == obs.IDs[i] 
    assert order == order_out

    # check size of object is same as input dims
    assert chunk.shape == obs.dims

    #check that the chunk is zero outside the bounds
    assert np.all(chunk[:chunk_bounds[2]-1,:] == 0)
    assert np.all(chunk[chunk_bounds[3]+1:,] == 0)
    assert np.all(chunk[:,:chunk_bounds[0]-1] == 0)
    assert np.all(chunk[:,chunk_bounds[1]+1:] == 0)


def test_disperse_chunk_null(observation):
    '''
    ensure bounds return None when all dispersion is off image
    '''
    obs = observation
    i = 0 #i==0 source happens to be far left on image
    order = 3
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)

    # manually change x,y offset because took transform from a real direct image, with different
    # pixel 0,0 than the mock data. This puts i=1, order 1 onto the real grism image
    obs.xoffset = 2200
    obs.yoffset = 1000

    disperse_chunk_args = [i, order, wmin, wmax, sens_waves, sens_resp,
                       obs.IDs[i], obs.xs[i], obs.ys[i], 
                       obs.fluxes, 
                       obs.seg_wcs, obs.grism_wcs, obs.dims, 
                       obs.extrapolate_sed, obs.xoffset, obs.yoffset]

    (chunk, chunk_bounds, sid, order_out) = obs.disperse_chunk(*disperse_chunk_args)

    assert chunk_bounds is None
    assert np.all(chunk == 0)
