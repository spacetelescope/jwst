"""Test for the detector1 pipeline using NIRSpec data in IRS2 mode. This takes
an uncal file and generates the stage 1 FITS files (rate) along with the
intermediate products."""

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_detector1pipeline(rtdata_module):
    """Run calwebb_detector1 pipeline on NIRSpec data with IRS2 readout mode."""
    rtdata = rtdata_module
    rtdata.get_data("nirspec/irs2/jw01335004001_03101_00002_nrs2_uncal.fits")

    Step.from_cmdline(
        [
            "calwebb_detector1",
            rtdata.input,
            "--steps.group_scale.save_results=True",
            "--steps.dq_init.save_results=True",
            "--steps.saturation.save_results=True",
            "--steps.superbias.save_results=True",
            "--steps.refpix.save_results=True",
            "--steps.rscd.save_results=True",
            "--steps.linearity.save_results=True",
            "--steps.dark_current.save_results=True",
            "--steps.jump.save_results=True",
            "--steps.ramp_fit.save_results=True",
            "--steps.gain_scale.save_results=True",
        ]
    )


@pytest.fixture(scope="module")
def run_detector1_with_clean_flicker_noise(rtdata_module, resource_tracker):
    """Run detector1 pipeline on NIRSpec IRS2 data with noise cleaning."""
    rtdata_module.get_data("nirspec/irs2/jw01335004001_03101_00002_nrs2_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files.
    # Run optional clean_flicker_noise step,
    # saving extra outputs and masking to science regions
    args = [
        "jwst.pipeline.Detector1Pipeline",
        rtdata_module.input,
        "--output_file=jw01335004001_03101_00002_nrs2_cfn",
        "--save_calibrated_ramp=True",
        "--steps.clean_flicker_noise.skip=False",
        "--steps.clean_flicker_noise.mask_science_regions=True",
        "--steps.clean_flicker_noise.save_results=True",
        "--steps.clean_flicker_noise.save_mask=True",
        "--steps.clean_flicker_noise.save_background=True",
        "--steps.clean_flicker_noise.save_noise=True",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_det1(log_tracked_resources, run_detector1_with_clean_flicker_noise):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    [
        "group_scale",
        "dq_init",
        "saturation",
        "superbias",
        "refpix",
        "linearity",
        "dark_current",
        "jump",
        "0_ramp_fit",
        "gain_scale",
        "rate",
    ],
)
def test_nirspec_irs2_detector1(
    run_detector1pipeline, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """
    Regression test of calwebb_detector1 pipeline performed on NIRSpec IRS2 data.
    """
    rtdata = rtdata_module

    output_filename = f"jw01335004001_03101_00002_nrs2_{suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_nirspec_irs2_detector1/{output_filename}")

    # Set tolerances so comparisons work across architectures
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "suffix",
    [
        "cfn_clean_flicker_noise",
        "mask",
        "flicker_bkg",
        "flicker_noise",
        "cfn_ramp",
        "cfn_rate",
        "cfn_rateints",
    ],
)
def test_nirspec_irs2_detector1_with_clean_flicker_noise(
    run_detector1_with_clean_flicker_noise, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Test detector1 pipeline for NIRSpec IRS2 data with noise cleaning."""
    rtdata = rtdata_module

    output_filename = f"jw01335004001_03101_00002_nrs2_{suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_nirspec_irs2_clean_flicker_noise/{output_filename}")

    # Set tolerances so comparisons work across architectures
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
