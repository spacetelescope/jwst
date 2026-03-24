import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_spec2(rtdata_module, resource_tracker):
    """Run stage 2 pipeline on NIRISS SOSS superstripe data."""
    rtdata = rtdata_module
    rtdata.get_data("niriss/soss/jwst_niriss_specprofile_0022.fits")
    rtdata.get_data(
        "niriss/soss/jw07073021001_04102_00001-seg001_nis_SUB680STRIPE_SOSS_rateints.fits"
    )

    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=True",
        "--steps.bkg_subtract.save_results=True",
        "--steps.flat_field.save_results=True",
        "--steps.extract_1d.override_specprofile=jwst_niriss_specprofile_0022.fits",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_spec2(log_tracked_resources, run_spec2):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["assign_wcs", "bsubints", "flat_field", "calints", "x1dints"])
def test_niriss_soss_superstripe_spec2(rtdata_module, run_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test for spec2 with NIRISS SOSS superstripe data."""
    rtdata = rtdata_module
    basename = "jw07073021001_04102_00001-seg001_nis_SUB680STRIPE_SOSS"
    output = f"{basename}_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_superstripe_spec2/{output}")

    # Ignore the custom specprofile reference file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_SPPROF")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
