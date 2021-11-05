import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.master_background import MasterBackgroundStep


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):

    rtdata = rtdata_module

    rtdata.get_asn("miri/mrs/miri_mrs_mbkg_0304_spec3_asn.json")

    MasterBackgroundStep.call(rtdata.input, save_results=True, suffix='master_background')

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("exp_seq", ["exp1", "exp2", "exp3"])
def test_miri_mrs_dedicated_mbkg(run_pipeline, fitsdiff_default_kwargs,
                                 exp_seq):
    """Run a test for MIRI MRS data with dedicated background exposures."""

    rtdata = run_pipeline
    rtdata.output = "miri_mrs_seq1_long_34_" + exp_seq + \
                    "_master_background.fits"

    # Get the truth file
    rtdata.get_truth(os.path.join(
        "truth/test_miri_mrs_dedicated_mbkg",
        "miri_mrs_seq1_long_34_" + exp_seq + "_master_background.fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
