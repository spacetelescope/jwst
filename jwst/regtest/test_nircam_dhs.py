import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_detector1pipeline(rtdata_module, resource_tracker):
    """Run calwebb_detector1 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.get_data("nircam/dhs/sub164stripe4_dhs_mock_dark.fits")
    rtdata.get_data("nircam/dhs/jw04453010001_02106_00001_nrca1_genheader_uncal.fits")

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
        "--steps.dark_current.override_dark=sub164stripe4_dhs_mock_dark.fits",
        "--steps.jump.save_results=True",
        "--steps.jump.rejection_threshold=50.0",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_det1(log_tracked_resources, run_detector1pipeline):
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
        "rate",
        "rateints",
    ],
)
def test_nircam_dhs_detector1(
    run_detector1pipeline, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Regression test of detector1 pipeline performed on NIRCam DHS mock data."""
    rtdata = rtdata_module
    rtdata.input = "jw04453010001_02106_00001_nrca1_genheader_uncal.fits"
    output = "jw04453010001_02106_00001_nrca1_genheader_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_dhs/{output}")

    # Ignore the custom dark because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_DARK")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
