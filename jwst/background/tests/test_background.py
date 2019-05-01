"""
Unit tests for background subtraction
"""
import pytest
from numpy.testing.utils import assert_allclose
from jwst.background import BackgroundStep

from jwst import datamodels
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
