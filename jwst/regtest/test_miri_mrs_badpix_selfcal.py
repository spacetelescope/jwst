import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.badpix_selfcal import BadpixSelfcalStep


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):

    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01204-o021_20240127t024203_spec2_00010_asn.json")
    BadpixSelfcalStep.call(rtdata.input, skip=False, save_results=True, flagfrac=0.0005)
    return rtdata


@pytest.mark.bigdata
def test_miri_mrs_badpix_selfcal(run_pipeline, fitsdiff_default_kwargs,):
    """Run a test for MIRI MRS data with dedicated background exposures."""

    rtdata = run_pipeline
    rtdata.output = "jw01204-o021_20240127t024203_spec2_00010_asn_0_badpixselfcalstep.fits"

    # Get the truth file
    rtdata.get_truth("truth/test_miri_mrs_badpix_selfcal/jw01204-o021_20240127t024203_spec2_00010_asn_0_badpixselfcalstep.fits")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()