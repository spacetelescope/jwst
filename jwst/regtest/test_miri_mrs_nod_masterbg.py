import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):

    rtdata = rtdata_module

    collect_pipeline_cfgs("config")

    rtdata.get_asn("miri/mrs/miri_mrs_mbkg_nodded_spec3_asn.json")

    args = ["config/master_background.cfg", rtdata.input,
            "--save_results=True"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("nod_seq", ["seq1", "seq2"])
def test_miri_mrs_nod_masterbg(run_pipeline, fitsdiff_default_kwargs, nod_seq):
    """Run a test for MIRI MRS data with nodded background exposures."""

    rtdata = run_pipeline
    rtdata.output = "miri_mrs_nod_" + nod_seq + \
                    "_short12_exp1_master_background.fits"

    # Get the truth file
    rtdata.get_truth(os.path.join(
        "truth/test_miri_mrs_nod_masterbg",
        "miri_mrs_nod_" + nod_seq + "_short12_exp1_master_background.fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
