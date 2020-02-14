import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):

    rtdata = rtdata_module

    collect_pipeline_cfgs("config")

    rtdata.get_asn("miri/lrs/miri_lrs_mbkg_dedicated_spec3_asn.json")

    args = ["config/master_background.cfg", rtdata.input,
            "--save_results=True"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
def test_miri_lrs_dedicated_mbkg(run_pipeline, fitsdiff_default_kwargs):
    """Run a test for MIRI LRS data with dedicated background exposures."""

    rtdata = run_pipeline
    rtdata.output = "miri_lrs_seq2_exp2_master_background.fits"

    # Get the truth file
    rtdata.get_truth(os.path.join(
                "truth/test_miri_lrs_dedicated_mbkg",
                "miri_lrs_seq2_exp2_master_background.fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
