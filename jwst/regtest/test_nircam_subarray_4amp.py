import os

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

"""
jw02459005001_03103_00001-seg001_nrca3_uncal.fits       the input (uncal) file
jw02459005001_03103_00001-seg001_nrca3_rate.fits        rate file
jw02459005001_03103_00001-seg001_nrca3_rateints.fits    rateints file
"""


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run calwebb_detector1 pipeline on NIRCAM subarray data."""
    rtdata = rtdata_module
    rtdata.get_data("nircam/subarray/jw02459005001_03103_00001-seg001_nrca3_uncal.fits")

    args = ["jwst.pipeline.Detector1Pipeline", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "output",
    [
        "jw02459005001_03103_00001-seg001_nrca3_rate.fits",
        "jw02459005001_03103_00001-seg001_nrca3_rateints.fits",
    ],
    ids=["rate", "rateints"],
)
def test_nircam_detector1_subarray(run_pipeline, fitsdiff_default_kwargs, output):
    """
    Regression test of calwebb_detector1 pipeline performed on NIRSpec data.
    """
    rtdata = run_pipeline
    rtdata.output = output
    rtdata.get_truth(os.path.join("truth/test_nircam_subarray_4amp", output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
