import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_tso_spec2(rtdata_module):
    """Run stage 2 pipeline on NIRISS SOSS TSO data with FULL subarray."""
    rtdata = rtdata_module

    # Run tso-spec2 pipeline on a _rateints file
    rtdata.get_data("niriss/soss/jw02113004001_02101_00001-seg001_nis_rateints.fits")
    args = ["calwebb_spec2", rtdata.input]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso_tso3(rtdata_module, run_tso_spec2, resource_tracker):
    """Run stage 3 pipeline on NIRISS SOSS TSO data with FULL subarray."""
    rtdata = rtdata_module

    # Get the level3 association json file (though not its members) and run
    # the tso3 pipeline on all the _calints file listed in association
    rtdata.get_data("niriss/soss/jw02113-o004_20250910t061917_tso3_00001_asn.json")
    args = ["calwebb_tso3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_niriss_soss_tso_stage2(rtdata_module, run_tso_spec2, fitsdiff_default_kwargs):
    """Test spec2 pipeline on NIRISS SOSS TSO FULL frame data."""
    rtdata = rtdata_module

    # only 1 file expected
    output = "jw02113004001_02101_00001-seg001_nis_calints.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_full/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_niriss_soss_tso_stage3(rtdata_module, run_tso_tso3, fitsdiff_default_kwargs):
    """Test tso3 pipeline on NIRISS SOSS TSO FULL frame data."""
    rtdata = rtdata_module

    # only 1 file expected
    output = "jw02113004001_02101_00001-seg001_nis_o004_crfints.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_full/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
