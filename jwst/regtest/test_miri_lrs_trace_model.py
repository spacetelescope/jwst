"""Regression tests for MIRI LRS trace modeling."""

import os

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = "truth/test_miri_lrs_trace_model"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_spec2_trace_model_fs(rtdata_module):
    """Run the Spec2Pipeline and create a trace model with oversampling."""
    rtdata = rtdata_module

    # Get the input rate file
    rtdata.get_asn("miri/lrs/jw01530-o005_20221202t204827_spec2_00001_asn.json")

    # Run the pipeline
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.adaptive_trace_model.skip=false",
        "--steps.adaptive_trace_model.save_results=true",
        "--steps.adaptive_trace_model.oversample=2.0",
        "--steps.resample_spec.pixel_scale_ratio=0.5",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope="module")
def run_spec2_trace_model_slitless(rtdata_module):
    """Run the Spec2Pipeline and create a trace model with oversampling."""
    rtdata = rtdata_module

    # Get the input rate file
    rtdata.get_data("miri/lrs/jw04496004001_03103_00001-seg001_mirimage_truncated_rateints.fits")

    # Run the pipeline
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.adaptive_trace_model.skip=false",
        "--steps.adaptive_trace_model.save_results=true",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.parametrize(
    "suffix",
    ["cal", "adaptive_trace_model", "s2d", "x1d"],
)
def test_miri_lrs_fs_spec2_trace_model(
    run_spec2_trace_model_fs, fitsdiff_default_kwargs, suffix, rtdata_module
):
    """Test trace model and oversample on a MIRI LRS FS exposure."""

    rtdata = run_spec2_trace_model_fs
    output = f"jw01530005001_03103_00001_mirimage_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "suffix",
    ["calints", "adaptive_trace_model", "x1dints"],
)
def test_miri_lrs_slitless_spec2_trace_model(
    run_spec2_trace_model_slitless, fitsdiff_default_kwargs, suffix, rtdata_module
):
    """Test trace model and oversample on a MIRI LRS slitless exposure."""

    rtdata = run_spec2_trace_model_slitless
    output = f"jw04496004001_03103_00001-seg001_mirimage_truncated_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
