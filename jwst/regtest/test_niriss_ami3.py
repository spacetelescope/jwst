import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run calwebb_ami3 on NIRISS AMI data."""
    rtdata = rtdata_module
    rtdata.get_asn("niriss/ami/ami3_test_asn.json")

    # Run the calwebb_ami3 pipeline on the association
    args = ["calwebb_ami3", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def run_step_with_cal(rtdata_module):
    rtdata = rtdata_module
    # run step?
    rtdata.get_data("niriss/ami/jw04478001001_03102_00001_nis_cal.fits")
    args = [
        'ami_analyze',
        rtdata.input,
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("obs, suffix", [("012", "ami-oi"), ("015", "psf-ami-oi")])
def test_niriss_ami3_exp(run_pipeline, obs, suffix, fitsdiff_default_kwargs):
    """Check exposure-level results of calwebb_ami3"""
    rtdata = run_pipeline

    output = f"jw01093{obs}001_03102_00001_nis_a3002_{suffix}.fits"
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


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ("ami-oi", "amimulti-oi", "amilg"))
def test_niriss_ami3_cal(run_step_with_cal, suffix, fitsdiff_default_kwargs):
    rtdata = run_step_with_cal

    output = f"jw04478001001_03102_00001_nis_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_niriss_ami3/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
