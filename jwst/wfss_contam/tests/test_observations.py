import warnings

import pytest
import numpy as np
import asdf
from astropy.convolution import convolve
from astropy.stats import sigma_clipped_stats
from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_allclose
from photutils.datasets import make_100gaussians_image
from photutils.segmentation import make_2dgaussian_kernel, SourceFinder

from jwst.wfss_contam.observations import background_subtract, _select_ids, Observation
from jwst.wfss_contam.disperse import disperse
from stdatamodels.jwst.datamodels import SegmentationMapModel, ImageModel, MultiSlitModel

DIR_IMAGE = "direct_image.fits"


# TODO: mock the WCS objects so the data/ directory can be removed


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
    with warnings.catch_warnings():
        # asdf.exceptions.AsdfPackageVersionWarning in oldestdeps job
        warnings.filterwarnings("ignore", message="File .* was created with extension URI .* which is not currently installed")
        with asdf.open(get_pkg_data_filename("data/segmentation_wcs.asdf", package="jwst.wfss_contam.tests")) as asdf_file:
            wcsobj = asdf_file.tree['wcs']
            model.meta.wcs = wcsobj

    return model


@pytest.fixture(scope='module')
def grism_wcs():
    with warnings.catch_warnings():
        # asdf.exceptions.AsdfPackageVersionWarning in oldestdeps job
        warnings.filterwarnings("ignore", message="File .* was created with extension URI .* which is not currently installed")
        with asdf.open(get_pkg_data_filename("data/grism_wcs.asdf", package="jwst.wfss_contam.tests")) as asdf_file:
            return asdf_file.tree['wcs']


@pytest.fixture(scope='module')
def observation(direct_image_with_gradient, segmentation_map, grism_wcs):
    '''
    set up observation object with mock data.
    direct_image_with_gradient still needs to be run to produce the file,
    even though it is not called directly
    '''
    filter_name = "F200W"
    seg = segmentation_map.data
    all_ids = np.array(list(set(np.ravel(seg))))
    source_ids = all_ids[50:52]
    obs = Observation(
        direct_image_with_gradient.data,
        segmentation_map,
        grism_wcs,
        filter_name,
        source_id=source_ids,
    )
    return obs


@pytest.mark.parametrize("source_id, expected", [(None, [1,2,3]), (2, [2]), (np.array([1,3]), [1,3])])
def test_select_ids(source_id, expected):
    all_ids = [1, 2, 3]
    assert _select_ids(source_id, all_ids) == expected
    assert isinstance(_select_ids(source_id, all_ids), list)


def test_select_ids_expected_raises():
    with pytest.raises(ValueError):
        _select_ids("all", [1, 2, 3])


def test_background_subtract(direct_image_with_gradient):

    data = direct_image_with_gradient.data
    subtracted_data = background_subtract(data)
    mean, median, stddev = sigma_clipped_stats(subtracted_data, sigma=3.0)
    assert_allclose(mean, 0.0, atol=0.2*stddev)


def test_create_pixel_list(observation, segmentation_map):
    '''
    create_pixel_list is called on initialization
    compare the pixel lists to values determined directly
    from segmentation map

    Note: still need test coverage for flux dictionary
    '''
    seg = segmentation_map.data
    all_ids = np.array(list(set(np.ravel(seg))))
    source_ids = all_ids[50:52]
    for i, source_id in enumerate(source_ids):
        # TODO: add here
        pass


def test_disperse_order(observation):

    obs = observation
    order = 1
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)

    # manually change x,y offset because took transform from a real direct image, with different
    # pixel 0,0 than the mock data. This puts i=1, order 1 onto the real grism image
    obs.xoffset = 2200
    obs.yoffset = 1000

    # shorten pixel list to make this test take less time
    obs.disperse_order(order, wmin, wmax, sens_waves, sens_resp)

    # test simulated image. should be mostly but not all zeros
    assert obs.simulated_image.shape == obs.dims
    assert not np.allclose(obs.simulated_image, 0.0)
    assert np.median(obs.simulated_image) == 0.0

    # test simulated slits and their associated metadata
    # only the second of the two obs ids is in the simulated image
    assert type(obs.simulated_slits) == MultiSlitModel
    # TODO: check more here


def test_disperse_oversample_same_result(grism_wcs, segmentation_map):
    '''
    Coverage for bug where wavelength oversampling led to double-counted fluxes

    note: segmentation_map fixture needs to be able to find module-scoped direct_image
    fixture, so it must be imported here
    '''

    # manual input of input params set the same as test_observations.py
    x0 = np.array([300.5])
    y0 = np.array([300.5])
    order = 1
    flxs = np.array([1.0])
    source_id = np.array([0])
    naxis = (300, 500)
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)
    seg_wcs = segmentation_map.meta.wcs
    0, (300, 500), 2, False,
    xoffset = 2200
    yoffset = 1000

    output_images = []
    for os in [2, 3]:
        src = disperse(
            x0,
            y0,
            flxs,
            source_id,
            order,
            wmin,
            wmax,
            sens_waves,
            sens_resp,
            seg_wcs,
            grism_wcs,
            naxis,
            oversample_factor=os,
            xoffset=xoffset,
            yoffset=yoffset,
        )
        output_images.append(src[source_id[0]]["image"])

    # different oversampling gives different effects at the ends
    # unsure if this is a bug or not, but the middle should definitely be the same
    assert_allclose(output_images[0][2:-2,:], output_images[1][2:-2], rtol=1e-5)


def test_construct_slitmodel(observation):
    '''
    test that the chunk is constructed correctly
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

    # disperse_chunk_args = [i, order, wmin, wmax, sens_waves, sens_resp]
    # (chunk, chunk_bounds, sid, order_out) = obs.disperse_chunk(*disperse_chunk_args)
    # TODO: replace this call with obs.disperse_one_order constrained to one source_id

    slit = obs.construct_slitmodel(chunk, chunk_bounds, sid, order_out)

    # check that the metadata is correct
    assert slit.xstart == chunk_bounds[0]
    assert slit.xsize == chunk_bounds[1] - chunk_bounds[0] + 1
    assert slit.ystart == chunk_bounds[2]
    assert slit.ysize == chunk_bounds[3] - chunk_bounds[2] + 1
    assert slit.source_id == sid
    assert slit.meta.wcsinfo.spectral_order == order_out
    assert np.allclose(slit.data, chunk[chunk_bounds[2]:chunk_bounds[3]+1, chunk_bounds[0]:chunk_bounds[1]+1])
