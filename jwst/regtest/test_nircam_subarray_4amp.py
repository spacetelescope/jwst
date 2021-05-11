import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step

"""
jw00617196001_02102_00001_nrca4_uncal.fits       the input (uncal) file
jw00617196001_02102_00001_nrca4_rate.fits        rate file
jw00617196001_02102_00001_nrca4_rateints.fits        rateints file
jw00617196001_02102_00001_nrca4_trapsfilled.fits        trapsfilled file
"""


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_detector1 pipeline on NIRCAM subarray data."""
    rtdata = rtdata_module
    rtdata.get_data("nircam/subarray/jw00617196001_02102_00001_nrca4_uncal.fits")

    args = ["jwst.pipeline.Detector1Pipeline", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("output", [
    'jw00617196001_02102_00001_nrca4_rate.fits',
    'jw00617196001_02102_00001_nrca4_rateints.fits',
    'jw00617196001_02102_00001_nrca4_trapsfilled.fits', ],
    ids=['rate', 'rateints', 'trapsfilled'])
def test_nircam_detector1_subarray(run_pipeline, fitsdiff_default_kwargs, output):
    """
    Regression test of calwebb_detector1 pipeline performed on NIRSpec data.
    """
    rtdata = run_pipeline
    rtdata.output = output
    rtdata.get_truth(os.path.join("truth/test_nircam_subarray_4amp", output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
