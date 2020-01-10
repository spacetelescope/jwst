import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):

    rtdata = rtdata_module

    collect_pipeline_cfgs("config")

    # This is the user-supplied background file.
    rtdata.get_data("miri/lrs/miri_lrs_bkg_x1d.fits")
    user_bkg = rtdata.input

    # This is the input file for the master_background step.
    rtdata.get_data("miri/lrs/miri_lrs_sci+bkg_cal.fits")

    args = ["config/master_background.cfg", rtdata.input,
            "--user_background=" + user_bkg,
            "--save_results=True"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
def test_miri_lrs_masterbg_user(run_pipeline, fitsdiff_default_kwargs):
    """Run a test for MIRI LRS data with a user-supplied background file."""

    rtdata = run_pipeline
    rtdata.output = "miri_lrs_sci+bkg_master_background.fits"

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_miri_lrs_masterbg_user",
                                  "miri_lrs_sci+bkg_master_background.fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
