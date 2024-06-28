import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

import stdatamodels.jwst.datamodels as dm

from jwst.flatfield import FlatFieldStep
from jwst.stpipe import Step


@pytest.mark.bigdata
def test_nirspec_image2(rtdata, fitsdiff_default_kwargs):
    rtdata.get_data("nirspec/imaging/jw84600010001_02102_00001_nrs2_rate.fits")

    args = ["calwebb_image2", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw84600010001_02102_00001_nrs2_cal.fits"

    rtdata.get_truth("truth/test_nirspec_image2/jw84600010001_02102_00001_nrs2_cal.fits")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_flat_field_step_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    data = rtdata.get_data('nirspec/imaging/usf_assign_wcs.fits')
    user_supplied_flat = rtdata.get_data('nirspec/imaging/usf_flat.fits')

    data_flat_fielded = FlatFieldStep.call(data, user_supplied_flat=user_supplied_flat)
    rtdata.output = 'flat_fielded_step_user_supplied.fits'
    data_flat_fielded.write(rtdata.output)

    rtdata.get_truth('truth/test_nirspec_image2/flat_fielded_step_user_supplied.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_ff_inv(rtdata, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    with dm.open(rtdata.get_data('nirspec/imaging/usf_assign_wcs.fits')) as data:
        flatted = FlatFieldStep.call(data)
        unflatted = FlatFieldStep.call(flatted, inverse=True)

        # flat fielding may set some new NaN values - ignore these in test
        is_nan = np.isnan(unflatted.data)
        assert np.allclose(data.data[~is_nan], unflatted.data[~is_nan]), 'Inversion failed'

        # make sure NaNs are only at do_not_use pixels
        assert np.all(unflatted.dq[is_nan] & dm.dqflags.pixel['DO_NOT_USE'])


@pytest.mark.bigdata
def test_correction_pars(rtdata, fitsdiff_default_kwargs):
    """Test use of correction parameters"""
    with dm.open(rtdata.get_data('nirspec/imaging/usf_assign_wcs.fits')) as data:

        # First use of FlatFieldStep will store the correction.
        # The next use will use that correction
        step = FlatFieldStep()
        flatted = step.run(data)
        assert step.correction_pars['flat'] is not None

        step.use_correction_pars = True
        reflatted = step.run(data)

        # flat fielding may set some new NaN values - ignore these in test
        is_nan = np.isnan(reflatted.data)
        assert np.allclose(flatted.data[~is_nan], reflatted.data[~is_nan]), 'Re-run with correction parameters failed'

        # make sure NaNs are only at do_not_use pixels
        assert np.all(reflatted.dq[is_nan] & dm.dqflags.pixel['DO_NOT_USE'])
