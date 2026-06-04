import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_miri_wfss_spec3(rtdata_module, resource_tracker):
    """Run the calwebb_spec3 pipeline on MIRI WFSS"""
    rtdata = rtdata_module

    # Get the level3 association file and run the spec3 pipeline on it.
    rtdata.get_asn("miri/wfss/jw09505-o001_spec3_00001_asn.json")
    args = ["calwebb_spec3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.mark.parametrize("suffix", ["x1d", "c1d"])
def test_miri_wfss_spec3(run_miri_wfss_spec3, rtdata_module, suffix, fitsdiff_default_kwargs):
    """Regression test of the calwebb_spec3 pipeline applied to MIRI WFSS data"""
    rtdata = rtdata_module
    rtdata.input = "jw09505-o001_spec3_00001_asn.json"
    output = "jw09505-o001_t001_miri_p750l_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_miri_wfss_spec3/{output}")

    # Compare the results
    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
