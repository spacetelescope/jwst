import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_coron3 on MIRI 4QPM coronographic data."""
    rtdata = rtdata_module
    rtdata.get_asn("miri/coron/miri_1065_coron3_asn.json")

    # Run the calwebb_coron3 pipeline on the association
    args = [
        "calwebb_coron3", rtdata.input,
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["psfalign", "psfsub"])
@pytest.mark.parametrize("exposure", ["roll1", "roll2"])
def test_miri_coron3_sci_exp(run_pipeline, suffix, exposure, fitsdiff_default_kwargs):
    """Check intermediate results of calwebb_coron3"""
    rtdata = run_pipeline

    output = "MIRIM_CORON1065-F1065C-" + exposure + "_a3001_" + suffix + ".fits"
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

    output = "coron1065_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_coron3/" + output)

    fitsdiff_default_kwargs['atol'] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
