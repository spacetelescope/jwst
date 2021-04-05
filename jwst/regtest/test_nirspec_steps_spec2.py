"""Regression tests for NIRSpec IFU"""
import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

import jwst.datamodels as dm
from jwst.flatfield import FlatFieldStep
from jwst.flatfield.flat_field import nirspec_ifu
from jwst.pathloss import PathLossStep

# Define artifactory source and truth
INPUT_PATH = 'nirspec/ifu'
TRUTH_PATH = 'truth/test_nirspec_ifu'


@pytest.mark.bigdata
def test_nirspec_ifu_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test using predefined interpolated flat"""
    with dm.open(rtdata.get_data('nirspec/ifu/nrs_ifu_nrs1_assign_wcs.fits')) as data:
        with dm.open(rtdata.get_data('nirspec/ifu/nrs_ifu_nrs1_interpolated_flat.fits')) as user_supplied_flat:
            nirspec_ifu(data, None, None, None, None, user_supplied_flat=user_supplied_flat)
            rtdata.output = 'ff_using_interpolated.fits'
            data.write(rtdata.output)

    rtdata.get_truth(TRUTH_PATH + '/' + 'ff_using_interpolated.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_flat_field_step_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    data = rtdata.get_data('nirspec/ifu/nrs_ifu_nrs1_assign_wcs.fits')
    user_supplied_flat = rtdata.get_data('nirspec/ifu/nrs_ifu_nrs1_interpolated_flat.fits')

    data_flat_fielded = FlatFieldStep.call(data, user_supplied_flat=user_supplied_flat)
    rtdata.output = 'flat_fielded_step_user_supplied.fits'
    data_flat_fielded.write(rtdata.output)
    del data_flat_fielded

    rtdata.get_truth(TRUTH_PATH + '/' + 'flat_fielded_step_user_supplied.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.slow
@pytest.mark.bigdata
def test_ff_inv(rtdata, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    with dm.open(rtdata.get_data('nirspec/ifu/nrs_ifu_nrs1_assign_wcs.fits')) as data:
        flatted = FlatFieldStep.call(data)
        unflatted = FlatFieldStep.call(flatted, inverse=True)

    assert np.allclose(data.data, unflatted.data), 'Inversion failed'


@pytest.mark.slow
@pytest.mark.bigdata
def test_pathloss_corrpars(rtdata):
    """Test PathLossStep using correction_pars"""
    with dm.open(rtdata.get_data('nirspec/ifu/nrs1_flat_field.fits')) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.use_correction_pars = True
        corrected_corrpars = pls.run(data)

    assert np.allclose(corrected.data, corrected_corrpars.data, equal_nan=True)


@pytest.mark.slow
@pytest.mark.bigdata
def test_pathloss_inverse(rtdata):
    """Test PathLossStep using correction_pars"""
    with dm.open(rtdata.get_data('nirspec/ifu/nrs1_flat_field.fits')) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.inverse = True
        corrected_inverse = pls.run(corrected)
        non_nan = ~np.isnan(corrected_inverse.data)

    assert np.allclose(corrected.data[non_nan], corrected_inverse.data[non_nan])


@pytest.mark.slow
@pytest.mark.bigdata
def test_pathloss_source_type(rtdata):
    """Test PathLossStep forcing source type"""
    with dm.open(rtdata.get_data('nirspec/ifu/nrs1_flat_field.fits')) as data:
        pls = PathLossStep()
        pls.source_type = 'extended'
        pls.run(data)

    assert np.allclose(pls.correction_pars.data,
                       pls.correction_pars.pathloss_uniform,
                       equal_nan=True)
