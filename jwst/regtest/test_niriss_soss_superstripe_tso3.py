import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_tso3(rtdata_module, resource_tracker):
    """Run stage 3 pipeline on NIRISS SOSS superstripe data."""
    rtdata = rtdata_module
    rtdata.get_data("niriss/soss/jwst_niriss_specprofile_0022.fits")
    rtdata.get_asn("niriss/soss/jw07073-o021_20251221t152654_tso3_00001_asn.json")
    args = [
        "calwebb_tso3",
        rtdata.input,
        "--steps.extract_1d.override_specprofile=jwst_niriss_specprofile_0022.fits",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_tso3(log_tracked_resources, run_tso3):
    log_tracked_resources()


def test_niriss_soss_superstripe_tso3_crfints(rtdata_module, run_tso3, fitsdiff_default_kwargs):
    """Regression test for tso3 crfints with NIRISS SOSS superstripe data."""
    rtdata = rtdata_module

    output = "jw07073021001_04102_00001-seg001_nis_SUB680STRIPE_SOSS_o021_crfints.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_superstripe_tso3/{output}")

    # Ignore the custom specprofile reference file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_SPPROF")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_niriss_soss_superstripe_tso3_x1dints(rtdata_module, run_tso3, fitsdiff_default_kwargs):
    """Regression test for tso3 x1dints with NIRISS SOSS superstripe data."""
    rtdata = rtdata_module

    output = "jw07073-o021_t002_niriss_clear-gr700xd-sub680stripe_x1dints.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_superstripe_tso3/{output}")

    # Ignore the custom specprofile reference file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_SPPROF")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_niriss_soss_superstripe_tso3_whtlt(rtdata_module, run_tso3, diff_astropy_tables):
    """Regression test for tso3 whitelight table with NIRISS SOSS superstripe data."""
    rtdata = rtdata_module

    output = "jw07073-o021_t002_niriss_clear-gr700xd-sub680stripe_whtlt.ecsv"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_superstripe_tso3/{output}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)
