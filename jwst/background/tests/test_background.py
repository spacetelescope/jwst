"""
Unit tests for background subtraction
"""
import os

from astropy.stats import sigma_clipped_stats
import pytest
import numpy as np
from numpy.testing import assert_allclose

from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.background import BackgroundStep
from jwst.background.tests import data as data_directory
from jwst.stpipe import Step
from jwst.background.background_sub import robust_mean, mask_from_source_cat, no_NaN


data_path = os.path.split(os.path.abspath(data_directory.__file__))[0]


def get_file_path(filename):
    """Construct an absolute path."""
    return os.path.join(data_path, filename)


@pytest.fixture(scope='module')
def background(tmpdir_factory):
    """Generate a  background image to feed to background step"""

    filename = tmpdir_factory.mktemp('background_input')
    filename = str(filename.join('background.fits'))
    with datamodels.IFUImageModel((10, 10)) as image:
        image.data[:, :] = 10
        image.meta.instrument.name = 'NIRSPEC'
        image.meta.instrument.detector = 'NRS1'
        image.meta.instrument.filter = 'CLEAR'
        image.meta.instrument.grating = 'PRISM'
        image.meta.exposure.type = 'NRS_IFU'
        image.meta.observation.date = '2019-02-27'
        image.meta.observation.time = '13:37:18.548'
        image.meta.date = '2019-02-27T13:37:18.548'

        image.meta.subarray.xstart = 1
        image.meta.subarray.ystart = 1

        image.meta.subarray.xsize = image.data.shape[-1]
        image.meta.subarray.ysize = image.data.shape[-2]

        image.meta.instrument.gwa_xtilt = 0.0001
        image.meta.instrument.gwa_ytilt = 0.0001
        image.meta.instrument.gwa_tilt = 37.0610

        image.save(filename)

    return filename


@pytest.fixture(scope='function')
def science_image():
    """Generate science image"""

    image = datamodels.IFUImageModel((10, 10))
    image.data[:, :] = 100
    image.meta.instrument.name = 'NIRSPEC'
    image.meta.instrument.detector = 'NRS1'
    image.meta.instrument.filter = 'CLEAR'
    image.meta.instrument.grating = 'PRISM'
    image.meta.exposure.type = 'NRS_IFU'
    image.meta.observation.date = '2019-02-27'
    image.meta.observation.time = '13:37:18.548'
    image.meta.date = '2019-02-27T13:37:18.548'
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1

    image.meta.subarray.xsize = image.data.shape[-1]
    image.meta.subarray.ysize = image.data.shape[-2]

    image.meta.instrument.gwa_xtilt = 0.0001
    image.meta.instrument.gwa_ytilt = 0.0001
    image.meta.instrument.gwa_tilt = 37.0610

    return image


def test_nirspec_gwa(_jail, background, science_image):
    """Verify NIRSPEC GWA logic for in the science and background"""

    # open the background to read in the GWA values
    back_image = datamodels.open(background)
    science_image.meta.instrument.gwa_xtilt = back_image.meta.instrument.gwa_xtilt
    science_image.meta.instrument.gwa_ytilt = back_image.meta.instrument.gwa_ytilt
    science_image.meta.instrument.gwa_tilt = back_image.meta.instrument.gwa_tilt

    bkg = [background]
    # Test Run with GWA values the same - confirm it runs
    # And gives the predicted result

    result = BackgroundStep.call(science_image, bkg)

    test = science_image.data - back_image.data
    assert_allclose(result.data, test)
    assert type(result) is type(science_image)
    assert result.meta.cal_step.back_sub == 'COMPLETE'
    back_image.close()


def test_nirspec_gwa_xtilt(_jail, background, science_image):
    """Verify NIRSPEC GWA Xtilt must be the same in the science and background image"""

    # open the background to read in the GWA values
    back_image = datamodels.open(background)
    science_image.meta.instrument.gwa_xtilt = back_image.meta.instrument.gwa_xtilt
    science_image.meta.instrument.gwa_ytilt = back_image.meta.instrument.gwa_ytilt
    science_image.meta.instrument.gwa_tilt = back_image.meta.instrument.gwa_tilt

    bkg = [background]

    # Test change xtilt
    science_image.meta.instrument.gwa_xtilt = \
        science_image.meta.instrument.gwa_xtilt + 0.00001

    result = BackgroundStep.call(science_image, bkg)

    assert type(result) is type(science_image)
    assert result.meta.cal_step.back_sub == 'SKIPPED'
    back_image.close()


def test_nirspec_gwa_ytilt(_jail, background, science_image):
    """Verify NIRSPEC GWA Ytilt must be the same in the science and background image"""

    # open the background to read in the GWA values
    back_image = datamodels.open(background)
    science_image.meta.instrument.gwa_xtilt = back_image.meta.instrument.gwa_xtilt
    science_image.meta.instrument.gwa_ytilt = back_image.meta.instrument.gwa_ytilt
    science_image.meta.instrument.gwa_tilt = back_image.meta.instrument.gwa_tilt

    bkg = [background]

    # Test different ytilt
    science_image.meta.instrument.gwa_ytilt = \
        science_image.meta.instrument.gwa_ytilt + 0.00001

    result = BackgroundStep.call(science_image, bkg)

    assert type(result) is type(science_image)
    assert result.meta.cal_step.back_sub == 'SKIPPED'

    back_image.close()


@pytest.fixture(scope='module')
def make_wfss_datamodel():
    """Generate WFSS Observation"""
    wcsinfo = {
        'dec_ref': -27.79156387419731,
        'ra_ref': 53.16247756038121,
        'roll_ref': 0.04254766236781744,
        'v2_ref': -290.1,
        'v3_ref': -697.5,
        'v3yangle': 0.56987,
        'vparity': -1}

    observation = {
        'date': '2023-01-05',
        'time': '8:59:37'}

    exposure = {
        'duration': 11.805952,
        'end_time': 58119.85416,
        'exposure_time': 11.776,
        'frame_time': 0.11776,
        'group_time': 0.11776,
        'groupgap': 0,
        'integration_time': 11.776,
        'nframes': 1,
        'ngroups': 8,
        'nints': 1,
        'nresets_between_ints': 0,
        'nsamples': 1,
        'sample_time': 10.0,
        'start_time': 58668.72509857639,
        'zero_frame': False}

    subarray = {'xsize': 2048,
                'ysize': 2048,
                'xstart': 1,
                'ystart': 1}

    instrument = {
        'filter_position': 1,
        'pupil_position': 1}

    image = datamodels.ImageModel((2048, 2048))

    image.meta.wcsinfo._instance.update(wcsinfo)
    image.meta.instrument._instance.update(instrument)
    image.meta.observation._instance.update(observation)
    image.meta.subarray._instance.update(subarray)
    image.meta.exposure._instance.update(exposure)
    image.data = np.random.rand(2048, 2048)
    image.meta.source_catalog = get_file_path('test_cat.ecsv')

    return image


filter_list = ['F250M', 'F277W', 'F335M', 'F356W', 'F460M',
               'F356W', 'F410M', 'F430M', 'F444W']  # + ['F480M', 'F322W2', 'F300M']


@pytest.mark.parametrize("pupils", ['GRISMC', 'GRISMR'])
@pytest.mark.parametrize("filters", filter_list)
@pytest.mark.parametrize("detectors", ['NRCALONG', 'NRCBLONG'])
def test_nrc_wfss_background(filters, pupils, detectors, make_wfss_datamodel):
    """Test background subtraction for NIRCAM WFSS modes."""
    data = make_wfss_datamodel

    data.meta.instrument.filter = filters
    data.meta.instrument.pupil = pupils
    data.meta.instrument.detector = detectors
    data.meta.instrument.channel = 'LONG'
    data.meta.instrument.name = 'NIRCAM'
    data.meta.exposure.type = 'NRC_WFSS'

    if data.meta.instrument.detector == 'NRCALONG':
        data.meta.instrument.module = 'A'
    elif data.meta.instrument.detector == 'NRCBLONG':
        data.meta.instrument.module = 'B'

    wcs_corrected = AssignWcsStep.call(data)

    # Get References
    wavelenrange = Step().get_reference_file(wcs_corrected, "wavelengthrange")
    bkg_file = Step().get_reference_file(wcs_corrected, 'wfssbkg')

    mask = mask_from_source_cat(wcs_corrected, wavelenrange)

    with datamodels.open(bkg_file) as bkg_ref:
        bkg_ref = no_NaN(bkg_ref)

        # calculate backgrounds
        pipeline_data_mean = robust_mean(wcs_corrected.data[mask])
        test_data_mean, _, _ = sigma_clipped_stats(wcs_corrected.data, sigma=2)

        pipeline_reference_mean = robust_mean(bkg_ref.data[mask])
        test_reference_mean, _, _ = sigma_clipped_stats(bkg_ref.data, sigma=2)

        assert np.isclose([pipeline_data_mean], [test_data_mean], rtol=1e-3)
        assert np.isclose([pipeline_reference_mean], [test_reference_mean], rtol=1e-1)


@pytest.mark.parametrize("filters", ['GR150C', 'GR150R'])
@pytest.mark.parametrize("pupils", ['F090W', 'F115W', 'F140M', 'F150W', 'F158M', 'F200W'])
def test_nis_wfss_background(filters, pupils, make_wfss_datamodel):
    """Test background subtraction for NIRISS WFSS modes."""
    data = make_wfss_datamodel

    data.meta.instrument.filter = filters
    data.meta.instrument.pupil = pupils
    data.meta.instrument.detector = 'NIS'
    data.meta.instrument.name = 'NIRISS'
    data.meta.exposure.type = 'NIS_WFSS'

    wcs_corrected = AssignWcsStep.call(data)

    # Get References
    wavelenrange = Step().get_reference_file(wcs_corrected, "wavelengthrange")
    bkg_file = Step().get_reference_file(wcs_corrected, 'wfssbkg')

    mask = mask_from_source_cat(wcs_corrected, wavelenrange)

    with datamodels.open(bkg_file) as bkg_ref:
        bkg_ref = no_NaN(bkg_ref)

        # calculate backgrounds
        pipeline_data_mean = robust_mean(wcs_corrected.data[mask])
        test_data_mean, _, _ = sigma_clipped_stats(wcs_corrected.data, sigma=2)

        pipeline_reference_mean = robust_mean(bkg_ref.data[mask])
        test_reference_mean, _, _ = sigma_clipped_stats(bkg_ref.data, sigma=2)

        assert np.isclose([pipeline_data_mean], [test_data_mean], rtol=1e-3)
        assert np.isclose([pipeline_reference_mean], [test_reference_mean], rtol=1e-1)


def test_robust_mean():
    """Test robust mean calculation"""
    data = np.random.rand(2048, 2048)
    result = robust_mean(data)
    test = np.mean(data)

    assert np.isclose([test], [result], rtol=1e-3)


def test_no_Nan():
    """Make sure that nan values are filled with fill value"""
    # Make data model
    model = datamodels.ImageModel()
    data = np.random.rand(10, 10)

    # Randomly insert NaNs
    data.ravel()[np.random.choice(data.size, 10, replace=False)] = np.nan
    model.data = data

    # Randomly select fill value
    fill_val = np.random.randint(0, 20)

    # Call no_NaN
    result = no_NaN(model, fill_value=fill_val)

    # Use np.nan to find NaNs.
    test_result = np.isnan(model.data)
    # Assign fill values to NaN indices
    model.data[test_result] = fill_val

    # Make sure arrays are equal.
    assert np.array_equal(model.data, result.data)
