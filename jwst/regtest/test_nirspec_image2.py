import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
import numpy as np

import stdatamodels.jwst.datamodels as dm

from jwst.flatfield import FlatFieldStep
from jwst.stpipe import Step


@pytest.mark.bigdata
def test_nirspec_image2(rtdata, fitsdiff_default_kwargs):
    basename = "jw03290001001_03101_00001_nrs2"
    rtdata.get_data(f"nirspec/imaging/{basename}_rate.fits")

    args = ["calwebb_image2", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = f"{basename}_cal.fits"

    rtdata.get_truth(f"truth/test_nirspec_image2/{basename}_cal.fits")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_flat_field_step_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    basename = "jw03290001001_03101_00001_nrs2"
    output_file = f"{basename}_flat_from_user_file.fits"
    data = rtdata.get_data(f"nirspec/imaging/{basename}_assign_wcs.fits")
    user_supplied_flat = rtdata.get_data(f"nirspec/imaging/{basename}_interpolatedflat.fits")

    data_flat_fielded = FlatFieldStep.call(data, user_supplied_flat=user_supplied_flat)
    rtdata.output = output_file
    data_flat_fielded.save(rtdata.output)

    rtdata.get_truth(f"truth/test_nirspec_image2/{output_file}")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_ff_inv(rtdata, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    input_file = "nirspec/imaging/jw03290001001_03101_00001_nrs2_assign_wcs.fits"
    with dm.open(rtdata.get_data(input_file)) as data:
        flatted = FlatFieldStep.call(data)
        unflatted = FlatFieldStep.call(flatted, inverse=True)

        # flat fielding may set some new NaN values - ignore these in test
        is_nan = np.isnan(unflatted.data)
        assert np.allclose(data.data[~is_nan], unflatted.data[~is_nan]), "Inversion failed"

        # make sure NaNs are only at do_not_use pixels
        assert np.all(unflatted.dq[is_nan] & dm.dqflags.pixel["DO_NOT_USE"])


@pytest.mark.bigdata
def test_correction_pars(rtdata, fitsdiff_default_kwargs):
    """Test use of correction parameters"""
    input_file = "nirspec/imaging/jw03290001001_03101_00001_nrs2_assign_wcs.fits"
    with dm.open(rtdata.get_data(input_file)) as data:
        # First use of FlatFieldStep will store the correction.
        # The next use will use that correction
        step = FlatFieldStep()
        flatted = step.run(data)
        assert step.correction_pars["flat"] is not None

        step.use_correction_pars = True
        reflatted = step.run(data)

        # flat fielding may set some new NaN values - ignore these in test
        is_nan = np.isnan(reflatted.data)
        assert np.allclose(flatted.data[~is_nan], reflatted.data[~is_nan]), (
            "Re-run with correction parameters failed"
        )

        # make sure NaNs are only at do_not_use pixels
        assert np.all(reflatted.dq[is_nan] & dm.dqflags.pixel["DO_NOT_USE"])
