import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_coron3 on coronographic data."""
    rtdata = rtdata_module
    rtdata.get_asn("nircam/coron/jw01386-c1020_20220909t073458_coron3_002a_asn.json")

    # Run the calwebb_coron3 pipeline on the association
    args = ["calwebb_coron3", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["psfalign", "psfsub", "crfints"])
@pytest.mark.parametrize("obs", ["002", "003"])
def test_nircam_coron3_sci_exp(run_pipeline, suffix, obs, fitsdiff_default_kwargs):
    """Check intermediate results of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw01386" + obs + "001_03108_00001_nrcalong_c1020_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nircam_coron3/" + output)

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["crfints"])
@pytest.mark.parametrize("exposure", ["00001", "00002", "00003", "00004"])
def test_nircam_coron3_psf_exp(run_pipeline, suffix, exposure, fitsdiff_default_kwargs):
    """Check intermediate results of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw01386001001_0310a_" + exposure + "_nrcalong_c1020_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nircam_coron3/" + output)

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["psfstack", "i2d"])
def test_nircam_coron3_product(run_pipeline, suffix, fitsdiff_default_kwargs):
    """Check final products of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw01386-c1020_t001_nircam_f410m-maskrnd-sub320a335r_mini_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nircam_coron3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
