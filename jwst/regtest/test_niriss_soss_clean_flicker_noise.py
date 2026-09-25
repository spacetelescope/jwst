"""Test for the detector1 pipeline with clean_flicker_noise, using NIRISS SOSS data."""

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_detector1_with_clean_flicker_noise(rtdata_module):
    """Run calwebb_detector1 pipeline on NIRISS SOSS data with clean_flicker_noise."""
    rtdata = rtdata_module
    rtdata.get_data("niriss/soss/jw02734002001_04101_00001-seg003_nis_uncal.fits")

    Step.from_cmdline(
        [
            "calwebb_detector1",
            rtdata.input,
            "--steps.clean_flicker_noise.skip=False",
            "--steps.clean_flicker_noise.save_results=True",
            "--steps.clean_flicker_noise.save_background=True",
            "--steps.clean_flicker_noise.save_mask=True",
            "--steps.clean_flicker_noise.save_noise=True",
        ]
    )


@pytest.mark.parametrize(
    "suffix",
    [
        "clean_flicker_noise",
        "mask",
        "flicker_bkg",
        "flicker_noise",
        "rateints",
    ],
)
def test_niriss_soss_detector1_with_clean_flicker_noise(
    run_detector1_with_clean_flicker_noise, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Test detector1 pipeline for NIRISS SOSS with 1/f noise cleaning."""
    rtdata = rtdata_module

    output_filename = f"jw02734002001_04101_00001-seg003_nis_{suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_niriss_soss_clean_flicker_noise/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
