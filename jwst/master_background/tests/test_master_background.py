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
    data.meta.model_type = 'MultiSpecModel'
    data.save(filename)

    return filename


def _generate_data():
    """Generate data of each type of input for master background step"""

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

    # hdulist = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    # hdulist[1].data = np.zeros((2048, 2048))
    # fs = datamodels.ImageModel(hdulist)
    # fs = AssignWcsStep.call(fs)
    # multislit = Extract2dStep.call(fs)

    # hdulist = create_nirspec_ifu_file("OPAQUE", "G140M")
    # hdulist[1].data = np.zeros((2048, 2048))
    # ifu_image = datamodels.IFUImageModel(hdulist)
    # ifu_image = AssignWcsStep.call(ifu_image)
    
    container = datamodels.ModelContainer([image])
    
    cube = datamodels.CubeModel((2, 10, 10))

    # Return the data and a status dependent on whether the step can process it
    return [(image, 'COMPLETE'),
            # (multislit, 'COMPLETE'),
            # (ifu_image, 'COMPLETE'),
            (container, 'COMPLETE'),
            (cube, 'SKIPPED'),
            ]

@pytest.mark.parametrize('input_data, status', _generate_data())
def test_master_background_init(input_data, status, _jail, user_background):
    """Verify data can run through the step"""

    result = MasterBackgroundStep.call(input_data)

    collect_pipeline_cfgs('./config')
    result = MasterBackgroundStep.call(
        input_data,
        config_file='config/master_background.cfg'
        )

    assert type(input_data) is type(result)

    try:
        assert result.meta.cal_step.master_back_sub == status
    except AttributeError:
        for model in result:
            assert model.meta.cal_step.master_back_sub == status

    # Run with a user-supplied background and verify this is recorded in header
    result = MasterBackgroundStep.call(input_data, user_background=user_background)

    try:
        if result.meta.cal_step.master_back_sub == 'COMPLETE':
            assert result.meta.master_background == 'user_background.fits'
    except AttributeError:
        for model in result:
            if model.meta.cal_step.master_back_sub == 'COMPLETE':
                assert model.meta.master_background == 'user_background.fits'

    # Make sure saving the computed background works
    result = MasterBackgroundStep.call(input_data, save_background=True)
