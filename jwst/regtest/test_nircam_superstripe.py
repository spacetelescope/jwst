from pathlib import Path

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]

# Note: these files are from an unrealistic simulation, extracted from an existing full
# frame exposure. The data shape and headers should be correct for the expected
# superstripe parameters but the rate values are not.
INPUT_FILES = [
    "jw01714006006_02101_00001_nrcb1_FULL_SUPSTP120_uncal.fits",  # FULL frame superstripe + refpix
    "jw01714006006_02101_00001_nrcb1_SUB64P_SUPSTP008_uncal.fits",  # SUB64P pure superstripe
    "jw01714006006_02101_00001_nrcb1_SUB64P_SUPSTP002_uncal.fits",  # SUB64P superstripe with repeat
]


@pytest.fixture(scope="module", params=INPUT_FILES)
def run_det1_nrc_superstripe(rtdata_module, request, resource_tracker):
    """Run calwebb_detector1 on NIRCam imaging in various superstripe modes."""
    rtdata = rtdata_module

    # Get the input exposure
    rtdata.get_data("nircam/superstripe/" + request.param)

    # Run detector1 pipeline
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.superbias.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.skip=True",  # no dark file available for this mode yet
        # "--steps.dark_current.save_results=True",
        "--steps.jump.save_results=True",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_nrc_superstripe_det1(
    log_tracked_resources, run_det1_nrc_superstripe
):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    [
        "dq_init",
        "saturation",
        "superbias",
        "refpix",
        "linearity",
        # "dark_current",
        "jump",
        "rate",
        "rateints",
    ],
)
def test_nircam_superstripe_det1(run_det1_nrc_superstripe, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 pipeline performed on NIRCam superstripe data."""
    rtdata = run_det1_nrc_superstripe
    input_name = Path(rtdata.input).name
    output = input_name.replace("uncal", suffix)
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_superstripe/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
