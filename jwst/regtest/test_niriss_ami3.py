import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step

from jwst.ami import AmiAnalyzeStep


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run calwebb_ami3 on NIRISS AMI data."""
    rtdata = rtdata_module
    rtdata.get_asn("niriss/ami/ami3_test_asn.json")

    # Run the calwebb_ami3 pipeline on the association
    args = ["calwebb_ami3", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("obs", ["012", "015"])
@pytest.mark.parametrize("suffix", ["ami-oi", "psf-ami-oi"])
def test_niriss_ami3_exp(run_pipeline, obs, fitsdiff_default_kwargs):
    """Check exposure-level results of calwebb_ami3"""
    rtdata = run_pipeline

    output = "jw01093" + obs + "001_03102_" + "00001_nis_a3002_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_niriss_ami3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_niriss_ami3_product(run_pipeline, fitsdiff_default_kwargs):
    """Check final products of calwebb_ami3"""
    rtdata = run_pipeline

    output = "jw01093012001_aminorm-oi.fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_niriss_ami3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

