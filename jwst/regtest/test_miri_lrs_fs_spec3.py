""" Test for the spec3 pipeline using MIR data in LRS mode. This takes
    the an association and generates the level 3 products."""
import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run the calwebb_spec3 pipeline on an ASN of nodded MIRI LRS
       fixedslit exposures."""

    rtdata = rtdata_module

    # Get the cfg files
    collect_pipeline_cfgs("config")

    # Get the spec3 ASN and its members
    rtdata.get_asn("miri/lrs/jw00623-o044_20191106t132852_spec3_001_asn.json")

    # Run the calwebb_spec3 pipeline; save results from intermediate steps
    args = ["config/calwebb_spec3.cfg", rtdata.input]

    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("output",["s2d", "x1d"])
def test_miri_lrs_fs_spec3(run_pipeline, fitsdiff_default_kwargs, output):
    """Regression test of the calwebb_spec3 pipeline on MIRI
       LRS fixedslit data using along-slit-nod pattern for
       background subtraction."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    rtdata.output = "jw00623-o044_t007_miri_p750l_" + output + ".fits"

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_miri_lrs_fs_spec3",
                                  "jw00623-o044_t007_miri_p750l_" + output + ".fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
