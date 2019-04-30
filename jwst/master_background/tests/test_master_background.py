"""
Unit tests for master background subtraction
"""
import pytest
import numpy as np

from jwst.master_background import MasterBackgroundStep
from jwst.assign_wcs import AssignWcsStep
# from jwst.assign_wcs.tests.test_nirspec import (
#     create_nirspec_fs_file,
#     create_nirspec_ifu_file,
#     )
from jwst.master_background.create_master_bkg import create_background
from jwst import datamodels
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
# from jwst.extract_2d import Extract2dStep


@pytest.fixture(scope='module')
def user_background(tmpdir_factory):
    """Generate a user background spectrum"""

    filename = tmpdir_factory.mktemp('master_background_user_input')
    filename = str(filename.join('user_background.fits'))
    wavelength = np.linspace(0.5, 25, num=100)
    flux = np.linspace(2.0, 2.2, num=100)
    data = create_background(wavelength, flux)
    data.save(filename)
    return filename


@pytest.fixture(scope='module')
def science_image(tmpdir_factory):
    """Generate science image """

    filename = tmpdir_factory.mktemp('master_background_image')
    filename = str(filename.join('science_image.fits'))

    image = datamodels.ImageModel((10, 10))
    image.meta.instrument.name = 'MIRI'
    image.meta.instrument.detector = 'MIRIMAGE'
    image.meta.exposure.type = 'MIR_LRS-FIXEDSLIT'
    image.meta.observation.date = '2018-01-01'
    image.meta.observation.time = '00:00:00'
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1
    image.meta.wcsinfo.v2_ref = 0
    image.meta.wcsinfo.v3_ref = 0
    image.meta.wcsinfo.roll_ref = 0
    image.meta.wcsinfo.ra_ref = 0
    image.meta.wcsinfo.dec_ref = 0
    image = AssignWcsStep.call(image)
    image.save(filename)
    return filename


def test_master_background_userbg(_jail, user_background, science_image):
    """Verify data can run through the step with a user-supplied background"""

    # Run with a user-supplied background and verify this is recorded in header
    collect_pipeline_cfgs('./config')
    result = MasterBackgroundStep.call(
        science_image,
        config_file='config/master_background.cfg',
        user_background=user_background,
        )

    # with the fixture science_image is a filename now - should
    # I use datamodels open to  compre
    # For inputs that are not files, the following should be true
    # assert type(science_image) is type(result)
    # assert result is not science_image
    assert result.meta.cal_step.master_background == 'COMPLETE'
    assert result.meta.background.master_background_file == 'user_background.fits'


def test_master_background_logic(_jail, user_background, science_image):
    """Verify if calspec 2 background step was run the master background step will be skipped"""

    image = datamodels.open(science_image)
    image.meta.cal_step.back_sub = 'COMPLETE'

    # Run with a user-supplied background
    collect_pipeline_cfgs('./config')
    result = MasterBackgroundStep.call(
        image,
        config_file='config/master_background.cfg',
        user_background=user_background,
        )

    assert result.meta.cal_step.master_background == 'SKIPPED'
    assert type(image) is type(result)

    # Now force it
    result = MasterBackgroundStep.call(
        image,
        config_file='config/master_background.cfg',
        user_background=user_background,
        force_subtract=True
        )

    assert result.meta.cal_step.master_background == 'COMPLETE'
    assert type(image) is type(result)
    image.close()
