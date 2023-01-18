import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_coron3 on MIRI 4QPM coronographic data."""
    rtdata = rtdata_module
    rtdata.get_asn("miri/coron/jw01386-c1002_20230109t015044_coron3_00001_asn.json")

    # Run the calwebb_coron3 pipeline on the association
    args = [
        "calwebb_coron3", rtdata.input,
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["crfints", "psfalign", "psfsub"])
@pytest.mark.parametrize("exposure", ["4", "5"])
def test_miri_coron3_sci_exp(run_pipeline, suffix, exposure, fitsdiff_default_kwargs):
    """Check intermediate results of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw0138600" + exposure + "001_04101_00001_mirimage_c1002_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_coron3/" + output)

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["psfstack", "i2d"])
def test_miri_coron3_product(run_pipeline, suffix, fitsdiff_default_kwargs):
    """Check final products of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw01386-c1002_t001_miri_f1140c-mask1140_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_coron3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
