import os
import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.master_background import MasterBackgroundStep


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    rtdata = rtdata_module

    # This is the user-supplied background file.
    rtdata.get_data("miri/lrs/jw01529-o004_t002_miri_p750l_x1d.fits")
    user_bkg = rtdata.input

    # This is the input file for the master_background step.
    rtdata.get_data("miri/lrs/jw01529003001_03103_00011_mirimage_cal.fits")

    MasterBackgroundStep.call(
        rtdata.input,
        user_background=user_bkg,
        save_results=True,
        suffix="user_mbsub",
        save_background=True,
    )

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "output",
    [
        "jw01529003001_03103_00011_mirimage_masterbg2d.fits",
        "jw01529003001_03103_00011_mirimage_user_mbsub.fits",
    ],
)
def test_miri_lrs_masterbg_user(run_pipeline, fitsdiff_default_kwargs, output):
    """Run a test for MIRI LRS data with a user-supplied background file."""

    rtdata = run_pipeline
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_miri_lrs_masterbg_user", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
