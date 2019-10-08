"""
Unit tests for background subtraction
"""
import pytest
import numpy as np
from numpy.testing.utils import assert_allclose

from jwst import datamodels
from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.util import create_grism_bbox
from jwst.background import BackgroundStep
from jwst.background.background_sub import robust_mean, subtract_wfss_bkg, mask_from_source_cat, no_NaN
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs


@pytest.fixture(scope='module')
def background(tmpdir_factory):
    """Generate a  background image to feed to background step"""

    filename = tmpdir_factory.mktemp('background_input')
    filename = str(filename.join('background.fits'))
    image = datamodels.IFUImageModel((10, 10))
    image.data[:,:] = 10
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

    image.meta.instrument.gwa_xtilt = 0.0001
    image.meta.instrument.gwa_ytilt = 0.0001
    image.meta.instrument.gwa_tilt = 37.0610

    image.save(filename)
    return filename


@pytest.fixture(scope='function')
def science_image():
    """Generate science image"""

    image = datamodels.IFUImageModel((10, 10))
    image.data[:,:] = 100
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

    collect_pipeline_cfgs('./config')
    result = BackgroundStep.call(
        science_image, bkg,
        config_file='config/background.cfg',
        )

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
    collect_pipeline_cfgs('./config')
    result = BackgroundStep.call(
        science_image, bkg,
        config_file='config/background.cfg',
        )
    assert type(result) is type(science_image)
    assert result.meta.cal_step.back_sub == 'SKIPPED'
    back_image.close()


def test_nirspec_gwa_ytitl(_jail, background, science_image):
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
    collect_pipeline_cfgs('./config')
    result = BackgroundStep.call(
        science_image, bkg,
        config_file='config/background.cfg',
        )
    assert type(result) is type(science_image)
    assert result.meta.cal_step.back_sub == 'SKIPPED'

    back_image.close()


def make_wfss_datamodel():
    """Generate NIRISS WFSS Observation"""

    wcsinfo = {
        'dec_ref': -27.79156387419731,
        'ra_ref': 53.16247756038121,
        'roll_ref': 0.04254766236781744,
        'v2_ref': -290.1,
        'v3_ref': -697.5,
        'v3yangle': 0.56987,
        'vparity': -1}

    observation = {
        'date': '2016-09-05',
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
        'type': 'NIS_WFSS',
        'zero_frame': False}

    subarray = {'xsize':2048,
                'ysize':2048}

    instrument = {
        'detector': 'NIS',
        'filter': 'GR150C',
        'pupil': 'F090W',
        'name': 'NIRISS',
        'filter_position':1,
        'pupil_position':1}

    image = datamodels.ImageModel((2048, 2048))

    image.meta.wcsinfo._instance.update(wcsinfo)
    image.meta.instrument._instance.update(instrument)
    image.meta.observation._instance.update(observation)
    image.meta.subarray._instance.update(subarray)
    image.meta.exposure._instance.update(exposure)
    image.data = np.random.rand(2048, 2048)
    image.meta.source_catalog.filename = '/Users/mfix/Desktop/swara_nis_cat.ecsv'

    return image

def test_wfss_background():
    """Test background subtraction for NIRISS WFSS modes."""
    data = make_wfss_datamodel()
    bkg = datamodels.ImageModel()
    
    wcs_corrected = AssignWcsStep.call(data)

    step_result = BackgroundStep.call(wcs_corrected, [bkg])
