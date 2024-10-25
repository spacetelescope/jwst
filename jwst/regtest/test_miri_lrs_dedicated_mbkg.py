import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.master_background import MasterBackgroundStep


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):

    rtdata = rtdata_module

    rtdata.get_asn("miri/lrs/jw01529-o003_spec3_with_bg_00001_asn.json")

    MasterBackgroundStep.call(rtdata.input, save_results=True, suffix='mbsub')

    return rtdata


@pytest.mark.bigdata
def test_miri_lrs_dedicated_mbkg(run_pipeline, fitsdiff_default_kwargs):
    """Run a test for MIRI LRS data with dedicated background exposures."""

    rtdata = run_pipeline
    rtdata.output = "jw01529003001_03103_00011_mirimage_mbsub.fits"

    # Get the truth file
    rtdata.get_truth(os.path.join(
        "truth/test_miri_lrs_dedicated_mbkg",
        "jw01529003001_03103_00011_mirimage_mbsub.fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
