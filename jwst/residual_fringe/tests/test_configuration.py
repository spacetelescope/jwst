"""
Unit test for Residual Fringe Correction for testing interface
"""

import pytest
import numpy as np
from jwst import datamodels
from jwst.residual_fringe import ResidualFringeStep
from jwst.residual_fringe import residual_fringe


@pytest.fixture(scope='function')
def miri_image():

    image = datamodels.IFUImageModel((20, 20))
    image.data = np.random.random((20, 20))
    image.meta.instrument.name = 'MIRI'
    image.meta.instrument.detector = 'MIRIFULONG'
    image.meta.exposure.type = 'MIR_MRS'
    image.meta.instrument.channel = '12'
    image.meta.instrument.band = 'SHORT'
    image.meta.filename = 'test_miri.fits'
    return image


@pytest.mark.xfail(resason='bug in parsing in parameters')
def test_call_residual_fringe(_jail,  miri_image):
    """ test defaults of step are set up and user input are defined correctly """

    # testing the ignore_regions_min
    # There has to be an equal number of min and max ignore region values
    # --ignore_region_min="4.9,"  --ignore_region_max='5.5,"

    step = ResidualFringeStep()
    step.ignore_region_min = [4.9, 5.7]
    step.ignore_region_max = [5.6, 6.5, 9.0]

    # If the number ignore min and max regions is not the same a value error is returned
    with pytest.raises(ValueError):
        step.run(miri_image)


def test_fringe_flat_applied(_jail, miri_image):

    miri_image.meta.cal_step.fringe = 'SKIP'
    residual_fringe_reference_file = None
    regions_reference_file = None
    save_intermediate_results = False
    transmission_level = 2
    ignore_regions = {}
    pars = {'save_intermediate_results': save_intermediate_results,
            'transmission_level': transmission_level}

    rfc = residual_fringe.ResidualFringeCorrection(miri_image,
                                                   residual_fringe_reference_file,
                                                   regions_reference_file,
                                                   ignore_regions,
                                                   **pars)
    # test that the fringe flat step has to be already run on the data before running residual fringe step

    with pytest.raises(residual_fringe.ErrorNoFringeFlat):
        rfc.do_correction()
