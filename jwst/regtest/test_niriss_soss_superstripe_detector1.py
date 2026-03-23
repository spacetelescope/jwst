import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_detector1(rtdata_module, resource_tracker):
    """Run calwebb_detector1 on NIRISS SOSS superstripe data."""
    rtdata = rtdata_module
    rtdata.get_data("niriss/soss/jw07073021001_04102_00001-seg001_nis_SUB680STRIPE_SOSS_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.superbias.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.jump.save_results=True",
        "--save_calibrated_ramp=True",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_det1(log_tracked_resources, run_detector1):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    [
        "dq_init",
        "saturation",
        "superbias",
        "refpix",
        "linearity",
        "dark_current",
        "jump",
        "ramp",
        "rate",
        "rateints",
    ],
)
def test_niriss_soss_superstripe_detector1(
    rtdata_module, run_detector1, fitsdiff_default_kwargs, suffix
):
    """Regression test for detector1 with NIRISS superstripe data."""
    rtdata = rtdata_module
    basename = "jw07073021001_04102_00001-seg001_nis_SUB680STRIPE_SOSS"
    rtdata.input = f"{basename}_uncal.fits"
    output = f"{basename}_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_superstripe_det1/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
